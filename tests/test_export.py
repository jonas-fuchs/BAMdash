"""
Export and CLI-behavior tests for ``bamdash``.

Covers the surface not already exercised by ``test_multi_ref.py``:

* ``html._sanitize_ref_id`` — table-driven sanitization
* master-HTML ``bamdashAxisRanges`` JSON correctness (per-ref linear/log ranges
  match the independently-computed coverage bounds)
* static-image export: log y-axis applied via ``prepare_static`` while the
  interactive HTML keeps linear; ``--dimensions`` honored; deepcopy preserves
  the subplot grid (a track subplot survives into the static svg)
* ``--dump`` sidecar contents (bam_stats / vcf / bed / gb) and multi- vs
  single-reference naming
* CLI error/edge behaviour: empty ``--suffix``, unsupported suffix value,
  unsupported track extension, ``--binsize 0``, invalid ``--ref-id``,
  ``--dimensions`` warning when only html is requested, no-args prints help
"""
from __future__ import annotations

import json
import logging
import math
import re
from pathlib import Path

import pytest

from bamdash.command import main
from bamdash.scripts import html as html_mod
from tests.conftest import (
    make_bam,
    make_bed,
    make_gb,
    make_gb_record,
    make_vcf,
    run_cli,
)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _run(bam: Path, prefix: Path, extra_args: list[str]) -> Path:
    return run_cli(bam, prefix, extra_args)


def _capture_logs() -> tuple[list[str], logging.Handler]:
    records: list[str] = []
    handler = logging.Handler()
    handler.emit = lambda rec: records.append(rec.getMessage())
    logger = logging.getLogger("bamdash")
    logger.addHandler(handler)
    return records, handler


def _detach(handler: logging.Handler) -> None:
    logging.getLogger("bamdash").removeHandler(handler)


def _axis_ranges(html_text: str) -> dict:
    m = re.search(r"var bamdashAxisRanges = (\{.*?\});", html_text, re.DOTALL)
    assert m, "bamdashAxisRanges not found in master HTML"
    return json.loads(m.group(1))


# --------------------------------------------------------------------------- #
# _sanitize_ref_id (table-driven)
# --------------------------------------------------------------------------- #
class TestSanitizeRefId:
    @pytest.mark.parametrize("raw, expected", [
        ("refA", "refA"),
        ("ref A", "ref_A"),          # spaces -> single underscore
        ("ref/A", "ref_A"),          # slash -> underscore
        ("ref-1.2", "ref_1_2"),      # runs of non-alnum collapsed
        ("a b/c", "a_b_c"),          # mixed runs collapsed to one _
        ("__x__", "x"),              # leading/trailing stripped
        ("ref(1)", "ref_1"),         # parens -> underscore
        ("NC_001", "NC_001"),        # already-safe id unchanged
    ])
    def test_sanitize(self, raw, expected):
        assert html_mod._sanitize_ref_id(raw) == expected


# --------------------------------------------------------------------------- #
# bamdashAxisRanges JSON correctness
# --------------------------------------------------------------------------- #
class TestAxisRangesJson:
    def test_ranges_match_independent_coverage_bounds(self, two_ref_bam, tmp_out):
        import pysam
        # compute per-ref coverage maxima independently from the bam fixture
        bam = pysam.AlignmentFile(str(two_ref_bam), "rb")
        maxima = {}
        for ref, length in [("refA", 50), ("refB", 30)]:
            cov = bam.count_coverage(ref, quality_threshold=0)
            total = [sum(c[i] for c in cov) for i in range(length)]
            maxima[ref] = max(total)
        bam.close()
        out = _run(two_ref_bam, tmp_out, [])
        ranges = _axis_ranges(Path(f"{out}.html").read_text())
        for ref, cov_max in maxima.items():
            upper = cov_max + cov_max / 10
            log_upper = max(1, math.ceil(math.log10(upper)))
            assert ranges[ref]["linear"] == [0, upper]
            assert ranges[ref]["log"] == [0, log_upper]

    def test_ranges_keys_are_sanitized_ref_ids(self, two_ref_bam, tmp_out):
        out = _run(two_ref_bam, tmp_out, [])
        ranges = _axis_ranges(Path(f"{out}.html").read_text())
        # refA/refB are already safe => keys are the raw ids
        assert set(ranges) == {"refA", "refB"}

    def test_ranges_for_unsafe_ref_id_uses_sanitized_key(self, tmp_path, tmp_out):
        # a reference id with a space is sanitized in both the panel id and the
        # axis_ranges key
        bam = tmp_path / "b.bam"
        make_bam(bam, refs=[("ref A", 20)],
                 reads=[{"ref": "ref A", "start": 0, "length": 10, "mapq": 60}])
        out = _run(bam, tmp_out, [])
        ranges = _axis_ranges(Path(f"{out}.html").read_text())
        assert "ref_A" in ranges
        assert "ref A" not in ranges


# --------------------------------------------------------------------------- #
# Static-image export
# --------------------------------------------------------------------------- #
class TestStaticExport:
    def test_dimensions_honored_in_svg(self, two_ref_bam, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out,
                   ["-t", str(two_ref_bed), "--suffix", "svg",
                    "--dimensions", "800", "600"])
        svg = Path(f"{out}_refA.svg").read_text()
        w = re.search(r'width="(\d+)"', svg)
        h = re.search(r'height="(\d+)"', svg)
        assert w is not None and h is not None
        assert (w.group(1), h.group(1)) == ("800", "600")

    def test_default_dimensions_when_omitted(self, two_ref_bam, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out, ["-t", str(two_ref_bed), "--suffix", "svg"])
        svg = Path(f"{out}_refA.svg").read_text()
        w = re.search(r'width="(\d+)"', svg)
        h = re.search(r'height="(\d+)"', svg)
        # default is 1920x1080
        assert (w.group(1), h.group(1)) == ("1920", "1080")

    def test_static_log_axis_applied_while_html_stays_linear(
            self, two_ref_bam, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out,
                   ["-t", str(two_ref_bed), "--suffix", "html", "svg"])
        svg = Path(f"{out}_refA.svg").read_text()
        # a log axis with dtick=1 renders tick labels at powers of ten; for a
        # coverage max around 8.8 the log upper is 1, so "10" appears as a tick
        assert "10" in svg
        # the interactive HTML figure must NOT force log on its yaxis
        from tests.conftest import parse_plotly_payload
        _data, layout = parse_plotly_payload(
            Path(f"{out}.html").read_text(), "plotly-refA")
        assert layout.get("yaxis", {}).get("type") != "log"

    def test_static_svg_hides_stats_annotation(self, two_ref_bam, tmp_out):
        out = _run(two_ref_bam, tmp_out, ["--suffix", "svg"])
        svg = Path(f"{out}_refA.svg").read_text()
        # prepare_static sets annotations visible=False; the stats title text
        # (which contains "<b>reference:</b>") must not render in the static image
        assert "<b>reference:</b>" not in svg

    def test_deepcopy_preserves_track_subplot_in_static(
            self, two_ref_bam, two_ref_bed, tmp_out):
        # with a bed track the figure has a second subplot (row 2). deepcopy in
        # command.main preserves the subplot grid, so the static svg must still
        # render the second axis (anchored subplot) — checked via the presence
        # of subplot anchor metadata and the coverage x-axis title.
        out = _run(two_ref_bam, tmp_out, ["-t", str(two_ref_bed), "--suffix", "svg"])
        svg = Path(f"{out}_refA.svg").read_text()
        assert "genome position" in svg
        assert re.search(r"anchor", svg) is not None
        # contrast: a coverage-only static svg has only one subplot, so it has
        # fewer rendered path elements than the tracked one
        out_cov = _run(two_ref_bam, tmp_out, ["--suffix", "svg"])
        svg_cov = Path(f"{out_cov}_refA.svg").read_text()
        assert svg.count("<path") > svg_cov.count("<path")

    def test_multiple_static_suffixes_each_written(self, two_ref_bam, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out,
                   ["-t", str(two_ref_bed), "--suffix", "svg", "png"])
        assert Path(f"{out}_refA.svg").exists()
        assert Path(f"{out}_refA.png").exists()
        assert Path(f"{out}_refB.svg").exists()
        assert Path(f"{out}_refB.png").exists()


# --------------------------------------------------------------------------- #
# --dump sidecar contents
# --------------------------------------------------------------------------- #
class TestDumpContents:
    def test_bam_stats_tabular_has_expected_rows(self, two_ref_bam, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out, ["-t", str(two_ref_bed), "--dump"])
        stats = Path(f"{out}_refA_bam_stats.tabular").read_text()
        # the stat dict keys written as the first column
        for key in ["reference", "reference length (bp)", "mapped", "unmapped",
                    "total", "mean coverage", "% recovery >= 5x", "gc content (%)"]:
            assert key in stats

    def test_vcf_dump_drops_position_jittered_column(self, two_ref_bam, two_ref_vcf, tmp_out):
        out = _run(two_ref_bam, tmp_out, ["-t", str(two_ref_vcf), "--dump"])
        header = Path(f"{out}_refA_vcf_data_0.tabular").read_text().splitlines()[0]
        cols = header.split("\t")
        assert "position_jittered" not in cols
        assert "position" in cols

    def test_bed_dump_columns_and_values(self, two_ref_bam, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out, ["-t", str(two_ref_bed), "--dump"])
        text = Path(f"{out}_refA_bed_data_0.tabular").read_text()
        lines = text.splitlines()
        header = lines[0].split("\t")
        # the 'track' column is dropped before writing
        assert "track" not in header
        for col in ["start", "stop", "strand", "name", "score",
                    "mean coverage", "% recovery >= 5x"]:
            assert col in header
        # the data row for regionA1 (start 1 -> 1-based 2, stop 15)
        assert "regionA1" in text

    def test_gb_dump_is_valid_json_with_feature_structure(self, two_ref_bam, make_gb_path, tmp_out):
        gb = make_gb_path(
            [make_gb_record("refA", 50,
                            features=[{"type": "CDS", "location": [(1, 10, 1)],
                                       "qualifiers": {"product": ["protein"]}}])],
            "refA.gb")
        out = _run(two_ref_bam, tmp_out, ["-t", str(gb), "--dump"])
        raw = Path(f"{out}_refA_gb_data_0.json").read_text()
        obj = json.loads(raw)
        assert "CDS" in obj
        # the feature key is "<start> <stop>" and contains strand/product
        feat = obj["CDS"]
        key = next(iter(feat))
        assert "strand" in feat[key]
        assert feat[key]["product"] == "protein"

    def test_dump_naming_multi_vs_single_ref(self, two_ref_bam, two_ref_vcf, two_ref_bed, tmp_out):
        # multi-ref: ref token in sidecar names
        out_multi = _run(two_ref_bam, tmp_out, ["-t", str(two_ref_vcf), str(two_ref_bed), "--dump"])
        assert Path(f"{out_multi}_refA_bam_stats.tabular").exists()
        assert Path(f"{out_multi}_refB_bam_stats.tabular").exists()
        assert not Path(f"{out_multi}_bam_stats.tabular").exists()

    def test_dump_single_ref_uses_plain_names(self, two_ref_bam, two_ref_vcf, two_ref_bed, tmp_out):
        out = _run(two_ref_bam, tmp_out,
                   ["-r", "refA", "-t", str(two_ref_vcf), str(two_ref_bed), "--dump"])
        assert Path(f"{out}_bam_stats.tabular").exists()
        assert Path(f"{out}_vcf_data_0.tabular").exists()
        assert Path(f"{out}_bed_data_0.tabular").exists()
        assert not Path(f"{out}_refA_bam_stats.tabular").exists()


# --------------------------------------------------------------------------- #
# CLI error / edge behaviour
# --------------------------------------------------------------------------- #
class TestCliErrors:
    def test_empty_suffix_exits(self, two_ref_bam, tmp_out):
        with pytest.raises(SystemExit) as exc:
            main(["-b", str(two_ref_bam), "-p", str(tmp_out), "-s"])
        # argparse error => exit code 2
        assert exc.value.code == 2

    def test_unsupported_suffix_value_exits(self, two_ref_bam, tmp_out):
        with pytest.raises(SystemExit) as exc:
            main(["-b", str(two_ref_bam), "-p", str(tmp_out), "-s", "tiff"])
        assert exc.value.code == 2

    def test_unsupported_track_extension_exits(self, two_ref_bam, tmp_path, tmp_out):
        bad_track = tmp_path / "weird.txt"
        bad_track.write_text("not a supported track")
        with pytest.raises(SystemExit) as exc:
            main(["-b", str(two_ref_bam), "-p", str(tmp_out), "-t", str(bad_track)])
        assert exc.value.code == 1

    def test_binsize_zero_exits(self, two_ref_bam, tmp_out):
        with pytest.raises(SystemExit) as exc:
            main(["-b", str(two_ref_bam), "-p", str(tmp_out), "--binsize", "0"])
        assert exc.value.code == 1

    def test_invalid_ref_id_exits(self, two_ref_bam, tmp_out):
        with pytest.raises(SystemExit) as exc:
            main(["-b", str(two_ref_bam), "-p", str(tmp_out), "-r", "NOPE"])
        assert exc.value.code == 1

    def test_dimensions_warning_when_html_only(self, two_ref_bam, two_ref_bed, tmp_out):
        records, handler = _capture_logs()
        try:
            _run(two_ref_bam, tmp_out,
                 ["-t", str(two_ref_bed), "--dimensions", "800", "600"])
            assert any("--dimensions has no effect" in r for r in records)
        finally:
            _detach(handler)

    def test_no_args_prints_help_and_exits(self, capsys):
        with pytest.raises(SystemExit):
            main([])
        captured = capsys.readouterr()
        # argparse prints usage/help to stdout on missing required args
        assert "usage" in captured.out.lower() or "usage" in captured.err.lower()

    def test_version_flag_exits(self, capsys):
        with pytest.raises(SystemExit) as exc:
            main(["--version"])
        # argparse version action exits 0
        assert exc.value.code == 0
        captured = capsys.readouterr()
        assert "bamdash" in captured.out
