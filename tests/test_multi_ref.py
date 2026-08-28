"""
Tests for the multi-reference dropdown feature.

Rewritten to pytest style (functions + fixtures + plain assert). Fixtures and
synthetic-data builders live in ``conftest.py``; these tests exercise the CLI
end-to-end via ``bamdash.command.main`` against a synthetic 2-reference BAM,
multi-reference VCF, and multi-reference BED, all generated programmatically
so the suite is self-contained (no dependency on the shipped ``test/`` data).
"""
from __future__ import annotations

import logging
from pathlib import Path

import pytest

from bamdash.command import main
from tests.conftest import make_gb_record


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _run(bam: Path, prefix: Path, extra_args: list[str]) -> Path:
    args = ["-b", str(bam), "-p", str(prefix)] + list(extra_args)
    main(args)
    return prefix


def _capture_logs() -> tuple[list[str], logging.Handler]:
    """Attach a recording handler to the bamdash logger; return (records, handler)."""
    records: list[str] = []
    handler = logging.Handler()
    handler.emit = lambda rec: records.append(rec.getMessage())
    logger = logging.getLogger("bamdash")
    logger.addHandler(handler)
    return records, handler


def _detach(handler: logging.Handler) -> None:
    logging.getLogger("bamdash").removeHandler(handler)


# --------------------------------------------------------------------------- #
# -r / --ref-id: zero, one, many
# --------------------------------------------------------------------------- #
class TestRefIdArg:
    def test_no_ref_id_uses_all_references(self, two_ref_bam, two_ref_vcf, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out, ["-t", str(two_ref_vcf), str(two_ref_bed)])
        html = Path(f"{out}.html").read_text()
        assert "bamdash-ref-select" in html
        assert 'value="refA"' in html
        assert 'value="refB"' in html

    def test_single_ref_id_produces_master_html(self, two_ref_bam, two_ref_vcf, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out, ["-r", "refA", "-t", str(two_ref_vcf), str(two_ref_bed)])
        html = Path(f"{out}.html").read_text()
        # single reference also uses the master HTML wrapper (dropdown lists
        # the one available reference) for consistency with multi-ref
        assert "bamdash-ref-select" in html
        assert 'value="refA"' in html
        assert 'value="refB"' not in html

    def test_multiple_ref_ids_preserve_user_order(self, two_ref_bam, two_ref_vcf, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out, ["-r", "refB", "refA", "-t", str(two_ref_vcf), str(two_ref_bed)])
        html = Path(f"{out}.html").read_text()
        assert "bamdash-ref-select" in html
        # order follows the user-supplied order
        assert html.index('value="refB"') < html.index('value="refA"')


# --------------------------------------------------------------------------- #
# GenBank track handling across references
# --------------------------------------------------------------------------- #
class TestGenbankTracks:
    def test_one_gb_per_reference_all_matched(self, two_ref_bam, make_gb_path):
        gbA = make_gb_path([make_gb_record("refA", 50, features=[{"type": "CDS", "location": [(1, 10, 1)]}])], "refA.gb")
        gbB = make_gb_path([make_gb_record("refB", 30, features=[{"type": "CDS", "location": [(1, 10, 1)]}])], "refB.gb")
        out = _run(two_ref_bam, gbA.parent / "out", ["-t", str(gbA), str(gbB)])
        html = Path(f"{out}.html").read_text()
        assert 'id="plot-refA"' in html
        assert 'id="plot-refB"' in html
        assert "CDS" in html

    def test_single_multi_record_gb_all_matched(self, two_ref_bam, make_gb_path):
        """A single .gb file with both references: both refs get features.

        Regression for the ``break`` bug in ``genbank_to_dict`` that stopped
        scanning at the first non-matching record, so only the first reference
        in a multi-record file was ever found.
        """
        multi_gb = make_gb_path(
            [
                make_gb_record("refA", 50, features=[{"type": "CDS", "location": [(1, 10, 1)]}]),
                make_gb_record("refB", 30, features=[{"type": "CDS", "location": [(1, 10, 1)]}]),
            ],
            "multi.gb",
        )
        out = _run(two_ref_bam, multi_gb.parent / "out", ["-t", str(multi_gb)])
        html = Path(f"{out}.html").read_text()
        assert 'id="plot-refA"' in html
        assert 'id="plot-refB"' in html
        assert "CDS" in html

    def test_non_matching_gb_warns_once(self, two_ref_bam, make_gb_path, tmp_out):
        gbX = make_gb_path([make_gb_record("refX", 40)], "refX.gb")
        records, handler = _capture_logs()
        try:
            _run(two_ref_bam, tmp_out, ["-t", str(gbX)])
        finally:
            _detach(handler)
        gb_warnings = [m for m in records if "refX.gb" in m]
        assert len(gb_warnings) == 1

    def test_matching_gb_no_warning(self, two_ref_bam, make_gb_path, tmp_out):
        gbA = make_gb_path([make_gb_record("refA", 50)], "refA.gb")
        records, handler = _capture_logs()
        try:
            _run(two_ref_bam, tmp_out, ["-r", "refA", "-t", str(gbA)])
        finally:
            _detach(handler)
        assert [m for m in records if "does not contain" in m] == []


# --------------------------------------------------------------------------- #
# Invalid reference id
# --------------------------------------------------------------------------- #
class TestInvalidRef:
    def test_invalid_ref_id_exits_nonzero(self, two_ref_bam, tmp_out):
        with pytest.raises(SystemExit) as exc:
            main(["-b", str(two_ref_bam), "-r", "NOPE", "-p", str(tmp_out)])
        assert exc.value.code == 1


# --------------------------------------------------------------------------- #
# Static image export
# --------------------------------------------------------------------------- #
class TestStaticExport:
    def test_multi_ref_static_one_image_per_reference(self, two_ref_bam, two_ref_vcf, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out, ["-t", str(two_ref_vcf), str(two_ref_bed), "--suffix", "png"])
        assert Path(f"{out}_refA.png").exists()
        assert Path(f"{out}_refB.png").exists()

    def test_single_ref_static_uses_plain_name(self, two_ref_bam, two_ref_vcf, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out, ["-r", "refA", "-t", str(two_ref_vcf), str(two_ref_bed), "--suffix", "png"])
        # single ref => no ref token in the filename (backward compatible)
        assert Path(f"{out}.png").exists()
        assert not Path(f"{out}_refA.png").exists()


# --------------------------------------------------------------------------- #
# Per-reference record split in dump output
# --------------------------------------------------------------------------- #
class TestMultiRefRecordSplit:
    def test_vcf_records_split_per_reference(self, two_ref_bam, two_ref_vcf, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out, ["-t", str(two_ref_vcf), str(two_ref_bed), "--dump"])
        refA_vcf = Path(f"{out}_refA_vcf_data_0.tabular").read_text()
        refB_vcf = Path(f"{out}_refB_vcf_data_0.tabular").read_text()
        # refA has two variants (positions 5 and 10); refB has one (position 7)
        assert "5\tA\tT" in refA_vcf
        assert "10\tC\tG" in refA_vcf
        assert "7\tG\tA" not in refA_vcf
        assert "7\tG\tA" in refB_vcf
        assert "5\tA\tT" not in refB_vcf

    def test_bed_records_split_per_reference(self, two_ref_bam, two_ref_vcf, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out, ["-t", str(two_ref_vcf), str(two_ref_bed), "--dump"])
        refA_bed = Path(f"{out}_refA_bed_data_0.tabular").read_text()
        refB_bed = Path(f"{out}_refB_bed_data_0.tabular").read_text()
        assert "regionA1" in refA_bed
        assert "regionB1" not in refA_bed
        assert "regionB1" in refB_bed
        assert "regionA1" not in refB_bed

    def test_single_ref_dump_uses_plain_names(self, two_ref_bam, two_ref_vcf, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out, ["-r", "refA", "-t", str(two_ref_vcf), str(two_ref_bed), "--dump"])
        # single ref => no ref token in sidecar names (backward compatible)
        assert Path(f"{out}_bam_stats.tabular").exists()
        assert Path(f"{out}_vcf_data_0.tabular").exists()
        assert not Path(f"{out}_refA_bam_stats.tabular").exists()


# --------------------------------------------------------------------------- #
# Master HTML structure
# --------------------------------------------------------------------------- #
class TestMasterHtmlStructure:
    def test_master_html_has_panels_and_inlined_plotly(self, two_ref_bam, two_ref_vcf, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out, ["-t", str(two_ref_vcf), str(two_ref_bed)])
        html = Path(f"{out}.html").read_text()
        # two ref-panel divs, first one active
        assert 'id="plot-refA" class="ref-panel active"' in html
        assert 'id="plot-refB" class="ref-panel"' in html
        # plotly.js is inlined (offline-capable)
        assert "Plotly.newPlot" in html
        # no CDN script src (offline)
        assert 'src="https://cdn' not in html

    def test_master_html_has_no_dark_mode(self, two_ref_bam, two_ref_vcf, two_ref_bed, tmp_out):
        """Dark mode and its button were removed from the master HTML."""
        import re
        out = _run(two_ref_bam, tmp_out, ["-t", str(two_ref_vcf), str(two_ref_bed)])
        html = Path(f"{out}.html").read_text()
        # no dark template registration
        assert "plotly_dark_custom" not in html
        # no per-figure updatemenu buttons in the figure layout JSON. The
        # plotly.js bundle itself contains the string "updatemenus" (it is
        # a plotly component name), so check the per-figure Plotly.newPlot
        # layout payload instead, which is where our buttons would live.
        assert re.findall(r'"updatemenus":\s*\[', html) == []
        assert ">dark<" not in html

    def test_master_html_has_global_log_linear_toggle(self, two_ref_bam, two_ref_vcf, two_ref_bed, tmp_out):
        """The master HTML has a global log/linear toggle next to the dropdown."""
        out = _run(two_ref_bam, tmp_out, ["-t", str(two_ref_vcf), str(two_ref_bed)])
        html = Path(f"{out}.html").read_text()
        assert 'id="bamdash-scale-toggle"' in html
        assert "bamdashAxisRanges" in html
        assert '"linear"' in html
        assert '"log"' in html
        assert "applyScale" in html

    def test_single_ref_html_uses_master_with_toggle(self, two_ref_bam, two_ref_vcf, two_ref_bed, tmp_out):
        """Single-reference output is the master HTML with the global toggle."""
        out = _run(two_ref_bam, tmp_out, ["-r", "refA", "-t", str(two_ref_vcf), str(two_ref_bed)])
        html = Path(f"{out}.html").read_text()
        assert "bamdash-ref-select" in html
        assert 'id="bamdash-scale-toggle"' in html


# --------------------------------------------------------------------------- #
# --offline / --no-offline: plotly.js inlined vs CDN
# --------------------------------------------------------------------------- #
class TestOfflineOnlineArg:
    def test_offline_inlines_plotly_bundle(self, two_ref_bam, two_ref_vcf, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out, ["-t", str(two_ref_vcf), str(two_ref_bed), "--offline"])
        html = Path(f"{out}.html").read_text()
        # the plotly.js bundle is inlined in a <script> block (no src)
        assert "Plotly.newPlot" in html
        assert 'src="https://cdn' not in html
        # the inlined bundle is large (~3MB); a CDN-free offline file is big
        assert len(html) > 1_000_000

    def test_no_offline_loads_plotly_from_cdn(self, two_ref_bam, two_ref_vcf, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out, ["-t", str(two_ref_vcf), str(two_ref_bed), "--no-offline"])
        html = Path(f"{out}.html").read_text()
        # plotly.js is loaded from the CDN, not inlined
        import plotly
        assert f"cdn.plot.ly/plotly-{plotly.__version__}.min.js" in html
        # the CDN file is much smaller than the inlined bundle
        assert len(html) < 500_000

    def test_offline_is_default(self, two_ref_bam, tmp_out):
        """Without --offline/--no-offline, the plotly bundle is inlined."""
        out = _run(two_ref_bam, tmp_out, [])
        html = Path(f"{out}.html").read_text()
        assert 'src="https://cdn' not in html
        assert len(html) > 1_000_000
