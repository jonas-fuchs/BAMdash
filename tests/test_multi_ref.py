"""
Tests for the multi-reference dropdown feature.

Fixtures (a synthetic 2-reference BAM, a multi-reference VCF, and a
multi-reference BED) are generated programmatically in ``setUpClass`` so the
suite is self-contained and does not depend on the single-reference HEV
example data shipped in ``test/``.
"""
import unittest
from pathlib import Path

import pysam

from bamdash.command import main


def _make_multi_ref_bam(path: Path):
    """Create a tiny 2-reference BAM (refA=50bp, refB=30bp) with reads on both."""
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "refA", "LN": 50}, {"SN": "refB", "LN": 30}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for i in range(8):
            r = pysam.AlignedSegment()
            r.query_name = f"readA{i}"
            r.query_sequence = "ACGT" * 5
            r.flag = 0
            r.reference_id = 0
            r.reference_start = i * 2
            r.mapping_quality = 60
            r.cigar = [(pysam.CMATCH, 20)]
            r.query_qualities = pysam.qualitystring_to_array("I" * 20)
            bam.write(r)
        for i in range(4):
            r = pysam.AlignedSegment()
            r.query_name = f"readB{i}"
            r.query_sequence = "ACGTACGTACGTACGTACGTACGTACGTAC"
            r.flag = 0
            r.reference_id = 1
            r.reference_start = i * 3
            r.mapping_quality = 60
            r.cigar = [(pysam.CMATCH, 30)]
            r.query_qualities = pysam.qualitystring_to_array("I" * 30)
            bam.write(r)
    pysam.index(str(path))


def _make_multi_ref_vcf(path: Path):
    """VCF with records on both refA and refB."""
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=refA,length=50>\n"
        "##contig=<ID=refB,length=30>\n"
        '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "refA\t5\t.\tA\tT\t.\tPASS\tAF=0.5\n"
        "refA\t10\t.\tC\tG\t.\tPASS\tAF=0.3\n"
        "refB\t7\t.\tG\tA\t.\tPASS\tAF=0.8\n"
    )


def _make_multi_ref_bed(path: Path):
    """BED with records on both refA and refB."""
    path.write_text("refA\t1\t15\tregionA1\t0\t+\nrefB\t5\t25\tregionB1\t0\t-\n")


class MultiRefFixture(unittest.TestCase):
    """Shared setUp/tearDown for the multi-ref test classes."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(__file__).parent / "_fixtures_tmp"
        cls.tmp.mkdir(exist_ok=True)
        cls.bam = cls.tmp / "multi.bam"
        cls.vcf = cls.tmp / "multi.vcf"
        cls.bed = cls.tmp / "multi.bed"
        _make_multi_ref_bam(cls.bam)
        _make_multi_ref_vcf(cls.vcf)
        _make_multi_ref_bed(cls.bed)

    @classmethod
    def tearDownClass(cls):
        for p in cls.tmp.glob("*"):
            p.unlink()
        cls.tmp.rmdir()

    def _run(self, prefix: str, extra_args):
        """Invoke the bamdash CLI in-process and return the output prefix path."""
        out = self.tmp / prefix
        args = ["-b", str(self.bam), "-p", str(out)] + list(extra_args)
        main(args)
        return out


class TestRefIdArg(MultiRefFixture):
    """``-r/--ref-id`` accepts zero, one, or multiple values."""

    def test_no_ref_id_uses_all_references(self):
        out = self._run("all", ["-t", str(self.vcf), str(self.bed)])
        html = Path(f"{out}.html").read_text()
        # master HTML (>=2 refs) => dropdown with both references
        self.assertIn("bamdash-ref-select", html)
        self.assertIn('value="refA"', html)
        self.assertIn('value="refB"', html)

    def test_single_ref_id_produces_master_html(self):
        out = self._run("one", ["-r", "refA", "-t", str(self.vcf), str(self.bed)])
        html = Path(f"{out}.html").read_text()
        # single reference now also uses the master HTML wrapper (dropdown
        # lists the one available reference) for consistency with multi-ref
        self.assertIn("bamdash-ref-select", html)
        self.assertIn('value="refA"', html)
        self.assertNotIn('value="refB"', html)

    def test_multiple_ref_ids_produce_master_html(self):
        out = self._run("two", ["-r", "refB", "refA", "-t", str(self.vcf), str(self.bed)])
        html = Path(f"{out}.html").read_text()
        self.assertIn("bamdash-ref-select", html)
        # order follows the user-supplied order
        self.assertLess(html.index('value="refB"'), html.index('value="refA"'))


class TestInvalidRef(MultiRefFixture):
    """An unknown reference id exits non-zero with a clear error."""

    def test_invalid_ref_id_exits_nonzero(self):
        out = self.tmp / "bad"
        with self.assertRaises(SystemExit) as ctx:
            main(["-b", str(self.bam), "-r", "NOPE", "-p", str(out)])
        self.assertEqual(ctx.exception.code, 1)


class TestStaticExport(MultiRefFixture):
    """Static image export writes one file per reference when multi-ref."""

    def test_multi_ref_static_one_image_per_reference(self):
        out = self._run("static", ["-t", str(self.vcf), str(self.bed), "--suffix", "png"])
        self.assertTrue(Path(f"{out}_refA.png").exists())
        self.assertTrue(Path(f"{out}_refB.png").exists())

    def test_single_ref_static_uses_plain_name(self):
        out = self._run("static1", ["-r", "refA", "-t", str(self.vcf), str(self.bed), "--suffix", "png"])
        # single ref => no ref token in the filename (backward compatible)
        self.assertTrue(Path(f"{out}.png").exists())
        self.assertFalse(Path(f"{out}_refA.png").exists())


class TestMultiRefRecordSplit(MultiRefFixture):
    """A multi-reference VCF/BED is split per reference in the dump output."""

    def test_vcf_records_split_per_reference(self):
        out = self._run("split", ["-t", str(self.vcf), str(self.bed), "--dump"])
        refA_vcf = Path(f"{out}_refA_vcf_data_0.tabular").read_text()
        refB_vcf = Path(f"{out}_refB_vcf_data_0.tabular").read_text()
        # refA has two variants (positions 5 and 10); refB has one (position 7)
        self.assertIn("5\tA\tT", refA_vcf)
        self.assertIn("10\tC\tG", refA_vcf)
        self.assertNotIn("7\tG\tA", refA_vcf)
        self.assertIn("7\tG\tA", refB_vcf)
        self.assertNotIn("5\tA\tT", refB_vcf)

    def test_bed_records_split_per_reference(self):
        out = self._run("splitbed", ["-t", str(self.vcf), str(self.bed), "--dump"])
        refA_bed = Path(f"{out}_refA_bed_data_0.tabular").read_text()
        refB_bed = Path(f"{out}_refB_bed_data_0.tabular").read_text()
        self.assertIn("regionA1", refA_bed)
        self.assertNotIn("regionB1", refA_bed)
        self.assertIn("regionB1", refB_bed)
        self.assertNotIn("regionA1", refB_bed)

    def test_single_ref_dump_uses_plain_names(self):
        out = self._run("dump1", ["-r", "refA", "-t", str(self.vcf), str(self.bed), "--dump"])
        # single ref => no ref token in sidecar names (backward compatible)
        self.assertTrue(Path(f"{out}_bam_stats.tabular").exists())
        self.assertTrue(Path(f"{out}_vcf_data_0.tabular").exists())
        self.assertFalse(Path(f"{out}_refA_bam_stats.tabular").exists())


class TestMasterHtmlStructure(MultiRefFixture):
    """The master HTML contains one panel per reference and inlined plotly.js."""

    def test_master_html_has_panels_and_inlined_plotly(self):
        out = self._run("struct", ["-t", str(self.vcf), str(self.bed)])
        html = Path(f"{out}.html").read_text()
        # two ref-panel divs, first one active
        self.assertIn('id="plot-refA" class="ref-panel active"', html)
        self.assertIn('id="plot-refB" class="ref-panel"', html)
        # plotly.js is inlined (offline-capable)
        self.assertIn("Plotly.newPlot", html)
        # no CDN script src (offline)
        self.assertNotIn('src="https://cdn', html)

    def test_master_html_has_no_dark_mode(self):
        """Dark mode and its button were removed from the master HTML."""
        out = self._run("nodark", ["-t", str(self.vcf), str(self.bed)])
        html = Path(f"{out}.html").read_text()
        # no dark template registration
        self.assertNotIn("plotly_dark_custom", html)
        # no per-figure updatemenu buttons in the figure layout JSON. The
        # plotly.js bundle itself contains the string "updatemenus" (it is
        # a plotly component name), so check the per-figure Plotly.newPlot
        # layout payload instead, which is where our buttons would live.
        import re
        layouts = re.findall(r'"updatemenus":\s*\[', html)
        self.assertEqual(layouts, [], "no per-figure updatemenu buttons should be present")
        self.assertNotIn('>dark<', html)

    def test_master_html_has_global_log_linear_toggle(self):
        """The master HTML has a global log/linear toggle next to the dropdown."""
        out = self._run("toggle", ["-t", str(self.vcf), str(self.bed)])
        html = Path(f"{out}.html").read_text()
        # the toggle button exists
        self.assertIn('id="bamdash-scale-toggle"', html)
        # per-reference axis ranges are injected for the toggle JS
        self.assertIn("bamdashAxisRanges", html)
        self.assertIn('"linear"', html)
        self.assertIn('"log"', html)
        # the applyScale helper relayouts every panel
        self.assertIn("applyScale", html)

    def test_single_ref_html_uses_master_with_toggle(self):
        """Single-reference output is the master HTML with the global toggle."""
        out = self._run("onebtn", ["-r", "refA", "-t", str(self.vcf), str(self.bed)])
        html = Path(f"{out}.html").read_text()
        # single ref also uses the master wrapper with the global toggle
        self.assertIn("bamdash-ref-select", html)
        self.assertIn('id="bamdash-scale-toggle"', html)
        # no per-figure updatemenu buttons in the figure layout JSON
        import re
        layouts = re.findall(r'"updatemenus":\s*\[', html)
        self.assertEqual(layouts, [], "no per-figure updatemenu buttons should be present")
        # no dark button
        self.assertNotIn("plotly_dark_custom", html)
        self.assertNotIn('>dark<', html)


if __name__ == "__main__":
    unittest.main()
