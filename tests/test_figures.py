"""
Figure-creation correctness tests for ``bamdash.scripts.plotting`` and the
master-HTML assembly in ``bamdash.scripts.html``.

Figures are produced end-to-end via the CLI (``bamdash.command.main``) and
verified at the **HTML level**: the rendered ``{prefix}.html`` is read and the
embedded ``Plotly.newPlot("<div_id>", <data>, <layout>)`` JSON payload is
parsed (see ``conftest.parse_plotly_payload``) so we can assert on trace
types, x/y data, marker colors, legendgroups, hovertemplates, layout shapes
(stems), axis ranges and subplot structure — without introspecting live
figure objects.

Pure-function helpers (``adjust_array_min_distance``, ``split_vcf_df``) are
unit-tested directly.
"""
from __future__ import annotations

import math
import statistics
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import pytest
from plotly.subplots import make_subplots

from bamdash.scripts import config, plotting
from tests.conftest import (
    make_bam,
    make_bed,
    make_gb,
    make_gb_record,
    make_vcf,
    parse_plotly_payload,
    run_cli,
)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _fig_and_payload(bam, prefix, extra_args, ref="refA", div_id=None):
    """Run the CLI and return ``(data, layout)`` parsed from the HTML for *ref*."""
    from bamdash.scripts import html as html_mod
    if div_id is None:
        div_id = f"plotly-{html_mod._sanitize_ref_id(ref)}"
    out = run_cli(bam, prefix, extra_args)
    html_text = Path(f"{out}.html").read_text(encoding="utf-8")
    return parse_plotly_payload(html_text, div_id)


# --------------------------------------------------------------------------- #
# adjust_array_min_distance (pure function)
# --------------------------------------------------------------------------- #
class TestAdjustArrayMinDistance:
    def test_output_is_non_decreasing(self):
        out = plotting.adjust_array_min_distance([1, 1.1, 1.2, 5], min_distance=1.0,
                                                 max_values=[0, 10])
        assert out == sorted(out)

    def test_adjacent_gaps_meet_min_distance_when_bounds_allow(self):
        out = plotting.adjust_array_min_distance([1, 1.1, 1.2], min_distance=1.0,
                                                 max_values=[0, 100])
        for a, b in zip(out, out[1:]):
            assert b - a >= 1.0 - 1e-6

    def test_values_stay_within_bounds(self):
        out = plotting.adjust_array_min_distance([1, 1.1, 1.2], min_distance=2.0,
                                                 max_values=[0, 10])
        for v in out:
            assert -max_values_default_lo(out) <= v <= 10 + 1e-6

    def test_idempotent(self):
        vals = [1, 1.1, 1.2, 5]
        once = plotting.adjust_array_min_distance(list(vals), min_distance=1.0,
                                                  max_values=[0, 10])
        twice = plotting.adjust_array_min_distance(list(once), min_distance=1.0,
                                                   max_values=[0, 10])
        assert once == twice

    def test_empty_and_single_unchanged(self):
        assert plotting.adjust_array_min_distance([], 1.0, [0, 10]) == []
        assert plotting.adjust_array_min_distance([5], 1.0, [0, 10]) == [5]


def max_values_default_lo(mv):
    # the lower bound used by the function is -max_values[0]/50; helper for clarity
    return mv[0] / 50


# --------------------------------------------------------------------------- #
# split_vcf_df (pure function)
# --------------------------------------------------------------------------- #
class TestSplitVcfDf:
    def test_no_duplicates_returns_single_element_list(self):
        df = pd.DataFrame({"position": [1, 2, 3], "type": ["SNP"] * 3})
        out = plotting.split_vcf_df(df)
        assert len(out) == 1
        assert len(out[0]) == 3

    def test_duplicates_round_robin_into_max_count_subdfs(self):
        # positions 1 (x3), 2 (x1) => max_n=3
        df = pd.DataFrame({"position": [1, 1, 1, 2], "type": ["SNP"] * 4})
        out = plotting.split_vcf_df(df)
        assert len(out) == 3
        # union of all rows preserved
        total = sum(len(sub) for sub in out)
        assert total == 4
        # each sub-df has at most one row per position
        for sub in out:
            assert sub["position"].is_unique


# --------------------------------------------------------------------------- #
# Coverage trace correctness (via HTML payload)
# --------------------------------------------------------------------------- #
class TestCoverageTraces:
    def test_two_coverage_traces_area_and_average(self, two_ref_bam, tmp_out):
        data, _layout = _fig_and_payload(two_ref_bam, tmp_out, [], ref="refA")
        scatters = [t for t in data if t.get("type", "scatter") == "scatter"]
        # exactly 2 scatters in a coverage-only figure: area + average
        assert len(scatters) == 2
        area = next(t for t in scatters if t.get("fill") == "tonexty")
        avg = next(t for t in scatters if t.get("mode") == "lines+text")
        # area trace styling from config
        assert area["fillcolor"] == config.coverage_fill_color
        assert area["line"]["color"] == config.coverage_line_color
        assert area["legendgroup"] == "coverage"
        # average trace styling
        assert avg["line"]["color"] == config.average_line_color
        assert avg["line"]["dash"] == "dash"
        assert avg["legendgroup"] == "average"

    def test_area_x_y_match_coverage_df(self, two_ref_bam, tmp_out):
        # refA: 50bp, 8 reads of length 20 starting at i*2 -> rebuild coverage independently
        import pysam
        bam = pysam.AlignmentFile(str(two_ref_bam), "rb")
        cov = bam.count_coverage("refA", quality_threshold=0)
        total = [sum(c[i] for c in cov) for i in range(50)]
        bam.close()
        data, _ = _fig_and_payload(two_ref_bam, tmp_out, [], ref="refA")
        area = next(t for t in data if t.get("fill") == "tonexty")
        assert list(area["x"]) == list(range(1, 51))
        assert list(area["y"]) == total

    def test_average_line_uses_unbinned_mean(self, two_ref_bam, tmp_out):
        import pysam
        bam = pysam.AlignmentFile(str(two_ref_bam), "rb")
        cov = bam.count_coverage("refA", quality_threshold=0)
        total = [sum(c[i] for c in cov) for i in range(50)]
        bam.close()
        expected_mean = statistics.mean(total)
        data, _ = _fig_and_payload(two_ref_bam, tmp_out, [], ref="refA")
        avg = next(t for t in data if t.get("mode") == "lines+text")
        # x spans min..max position; y is [mean, mean]
        assert avg["x"] == [1, 50]
        assert avg["y"] == [expected_mean, expected_mean]
        # text label includes the rounded mean with 'x'
        assert any(f"{round(expected_mean)}x" in t for t in avg["text"])

    def test_hovertemplate_references_position_coverage_bases(self, two_ref_bam, tmp_out):
        data, _ = _fig_and_payload(two_ref_bam, tmp_out, [], ref="refA")
        area = next(t for t in data if t.get("fill") == "tonexty")
        ht = area["hovertemplate"]
        for field in ["position", "coverage", "percentage A", "percentage C",
                      "percentage G", "percentage T"]:
            assert field in ht


# --------------------------------------------------------------------------- #
# Binning correctness
# --------------------------------------------------------------------------- #
class TestBinning:
    def test_binsize_one_is_unbinned(self, two_ref_bam, tmp_out):
        data, _ = _fig_and_payload(two_ref_bam, tmp_out, ["--binsize", "1"], ref="refA")
        area = next(t for t in data if t.get("fill") == "tonexty")
        # unbinned => 50 points
        assert len(area["x"]) == 50

    def test_binsize_three_averages_per_bin(self, two_ref_bam, tmp_out):
        import pysam
        bam = pysam.AlignmentFile(str(two_ref_bam), "rb")
        cov = bam.count_coverage("refA", quality_threshold=0)
        total = [sum(c[i] for c in cov) for i in range(50)]
        bam.close()
        data, _ = _fig_and_payload(two_ref_bam, tmp_out, ["--binsize", "3"], ref="refA")
        area = next(t for t in data if t.get("fill") == "tonexty")
        # Binning reduces the number of reported points relative to unbinned (50)
        # and every reported coverage value is an integer (bin means are rounded).
        assert len(area["x"]) < 50
        assert all(float(v).is_integer() for v in area["y"])
        # Each reported y must equal the independently-computed mean of the raw
        # coverage over the corresponding 3bp bin (rounded to 0 dp). The
        # implementation reports x as the 1-based *end* position of each bin, so
        # the bin covers 0-based indices [x-3, x-2, x-1] == total[x-3:x].
        for x_val, y_val in zip(area["x"], area["y"]):
            x = int(x_val)
            bin_vals = total[x - 3:x]
            expected = round(sum(bin_vals) / len(bin_vals))
            assert y_val == expected


# --------------------------------------------------------------------------- #
# upper / log_upper via axis_ranges JSON in master HTML
# --------------------------------------------------------------------------- #
class TestAxisRanges:
    def test_linear_and_log_upper_match_independent_computation(self, two_ref_bam, tmp_out):
        import json
        import pysam
        bam = pysam.AlignmentFile(str(two_ref_bam), "rb")
        cov = bam.count_coverage("refA", quality_threshold=0)
        total = [sum(c[i] for c in cov) for i in range(50)]
        bam.close()
        upper = max(total) + max(total) / 10
        log_upper = max(1, math.ceil(math.log10(upper)))
        out = run_cli(two_ref_bam, tmp_out, [])
        html_text = Path(f"{out}.html").read_text()
        # extract bamdashAxisRanges JSON
        import re
        m = re.search(r"var bamdashAxisRanges = (\{.*?\});", html_text, re.DOTALL)
        assert m, "bamdashAxisRanges not found"
        ranges = json.loads(m.group(1))
        assert ranges["refA"]["linear"] == [0, upper]
        assert ranges["refA"]["log"] == [0, log_upper]


# --------------------------------------------------------------------------- #
# VCF marker / stem correctness
# --------------------------------------------------------------------------- #
class TestVcfPlot:
    def test_one_marker_scatter_per_present_mutation_type(self, two_ref_bam, two_ref_vcf, tmp_out):
        # the two_ref_vcf has only SNPs on refA => only the SNP marker trace
        data, _ = _fig_and_payload(two_ref_bam, tmp_out, ["-t", str(two_ref_vcf)], ref="refA")
        markers = [t for t in data if t.get("mode") == "markers" and t.get("legendgroup") in ("SNP", "INS", "DEL")]
        types = {m["legendgroup"] for m in markers}
        assert types == {"SNP"}
        # SNP marker uses config color
        assert markers[0]["marker"]["color"] == config.snp_color

    def test_marker_y_uses_AF_when_present(self, two_ref_bam, two_ref_vcf, tmp_out):
        data, _ = _fig_and_payload(two_ref_bam, tmp_out, ["-t", str(two_ref_vcf)], ref="refA")
        snp = next(t for t in data if t.get("legendgroup") == "SNP")
        # refA variants have AF 0.5 and 0.3 (rounded to 2dp)
        assert sorted(snp["y"]) == [0.3, 0.5]

    def test_marker_y_ones_when_AF_absent(self, tmp_path, tmp_out):
        # VCF without AF INFO
        bam = tmp_path / "b.bam"
        make_bam(bam, refs=[("ref", 20)], reads=[{"ref": "ref", "start": 0, "length": 10, "mapq": 60}])
        vcf = tmp_path / "v.vcf"
        make_vcf(vcf, records=[{"chrom": "ref", "pos": 5, "ref": "A", "alt": "T"}])
        data, _ = _fig_and_payload(bam, tmp_out, ["-t", str(vcf)], ref="ref")
        snp = next(t for t in data if t.get("legendgroup") == "SNP")
        assert snp["y"] == [1]

    def test_stems_present_as_layout_shapes(self, two_ref_bam, two_ref_vcf, tmp_out):
        _data, layout = _fig_and_payload(two_ref_bam, tmp_out, ["-t", str(two_ref_vcf)], ref="refA")
        shapes = layout.get("shapes", [])
        # 2 variants * 3 segments each = 6 stem shapes
        assert len(shapes) == 6
        for s in shapes:
            assert s["line"]["color"] == config.stem_color
            assert s["line"]["width"] == config.stem_width
            assert s["layer"] == "below"

    def test_position_jittered_used_for_marker_x(self, two_ref_bam, two_ref_vcf, tmp_out):
        data, _ = _fig_and_payload(two_ref_bam, tmp_out, ["-t", str(two_ref_vcf)], ref="refA")
        snp = next(t for t in data if t.get("legendgroup") == "SNP")
        # raw positions are 5 and 10; jittered values should be close but may differ
        assert len(snp["x"]) == 2


# --------------------------------------------------------------------------- #
# Track box correctness (gb / bed)
# --------------------------------------------------------------------------- #
class TestTrackPlot:
    def test_bed_box_traces_fill_toself_with_cycled_alpha(self, two_ref_bam, two_ref_bed, tmp_out):
        data, _ = _fig_and_payload(two_ref_bam, tmp_out, ["-t", str(two_ref_bed)], ref="refA")
        # box traces are lines with fill=toself
        boxes = [t for t in data if t.get("fill") == "toself"]
        assert len(boxes) >= 1
        for b in boxes:
            assert b["mode"] == "lines"
            # fillcolor is an rgba string derived from config.track_color_single
            assert b["fillcolor"].startswith("rgba(")

    def test_strand_marker_symbol_matches_strand(self, two_ref_bam, two_ref_bed, tmp_out):
        data, _ = _fig_and_payload(two_ref_bam, tmp_out, ["-t", str(two_ref_bed)], ref="refA")
        # refA bed region is '+' => triangle-right
        markers = [t for t in data if t.get("mode") == "markers" and "hovertext" in t]
        plus = [m for m in markers if m["marker"]["symbol"] == config.strand_types[0]]
        assert len(plus) == 1
        assert plus[0]["showlegend"] is False

    def test_strand_marker_at_feature_midpoint(self, two_ref_bam, two_ref_bed, tmp_out):
        data, _ = _fig_and_payload(two_ref_bam, tmp_out, ["-t", str(two_ref_bed)], ref="refA")
        markers = [t for t in data if t.get("mode") == "markers" and "hovertext" in t]
        # regionA1: start 1, stop 15 (0-based 1 -> 1-based 2? bed_to_dict adds +1)
        # midpoint = global_start + (global_stop-global_start)/2
        m = markers[0]
        start, stop = 2, 15  # 0-based start 1 -> 1-based 2; stop 15 unchanged
        expected_x = start + (stop - start) / 2
        assert m["x"] == [expected_x]

    def test_single_feature_type_uses_track_color_single(self, two_ref_bam, two_ref_bed, tmp_out):
        # bed has a single "bed annotations" feature type => single color
        data, _ = _fig_and_payload(two_ref_bam, tmp_out, ["-t", str(two_ref_bed)], ref="refA")
        boxes = [t for t in data if t.get("fill") == "toself"]
        # the base color is config.track_color_single (rgb(145,145,145)) with alpha applied
        assert "145, 145, 145" in boxes[0]["fillcolor"]

    def test_gb_box_traces_present(self, two_ref_bam, make_gb_path, tmp_out):
        gbA = make_gb_path([make_gb_record("refA", 50, features=[{"type": "CDS", "location": [(1, 10, 1)]}])], "refA.gb")
        data, _ = _fig_and_payload(two_ref_bam, tmp_out, ["-t", str(gbA)], ref="refA")
        boxes = [t for t in data if t.get("fill") == "toself"]
        assert len(boxes) >= 1
        # legendgrouptitle text is the feature type
        assert any(b.get("legendgrouptitle", {}).get("text") == "CDS" for b in boxes)

    def test_track_yaxis_hidden_and_reversed(self, two_ref_bam, two_ref_bed, tmp_out):
        _data, layout = _fig_and_payload(two_ref_bam, tmp_out, ["-t", str(two_ref_bed)], ref="refA")
        # the bed track is row 2 => yaxis2; it should be hidden and autorange reversed
        yaxis2 = layout.get("yaxis2", {})
        assert yaxis2.get("visible") is False
        assert yaxis2.get("autorange") == "reversed"


# --------------------------------------------------------------------------- #
# build_figure layout
# --------------------------------------------------------------------------- #
class TestBuildFigureLayout:
    def test_template_and_hovermode(self, two_ref_bam, tmp_out):
        _data, layout = _fig_and_payload(two_ref_bam, tmp_out, [], ref="refA")
        assert layout["template"]["layout"].get("data") is not None or layout.get("template") is not None
        assert layout["hovermode"] == "closest"

    def test_stats_annotation_present(self, two_ref_bam, tmp_out):
        _data, layout = _fig_and_payload(two_ref_bam, tmp_out, [], ref="refA")
        anns = layout.get("annotations", [])
        assert len(anns) >= 1
        # the stats annotation is at y=1.14, yref=paper, no arrow
        stats_ann = next(a for a in anns if a.get("y") == 1.14 and a.get("yref") == "paper")
        assert stats_ann["showarrow"] is False

    def test_no_slider_sets_genome_position_title(self, two_ref_bam, tmp_out):
        _data, layout = _fig_and_payload(two_ref_bam, tmp_out, [], ref="refA")
        # coverage-only => 1 row => xaxis title is "genome position"
        assert layout["xaxis"]["title"]["text"] == "genome position"

    def test_slider_adds_rangeslider(self, two_ref_bam, tmp_out):
        _data, layout = _fig_and_payload(two_ref_bam, tmp_out, ["--slider"], ref="refA")
        rs = layout["xaxis"].get("rangeslider", {})
        assert rs.get("visible") is True
        assert rs.get("thickness") == 0.05

    def test_x_range_padded_by_max_position_over_50(self, two_ref_bam, tmp_out):
        _data, layout = _fig_and_payload(two_ref_bam, tmp_out, [], ref="refA")
        # refA length 50 => pad = 50/50 = 1 => range [-1, 51]
        assert layout["xaxis"]["range"] == [-1.0, 51.0]


# --------------------------------------------------------------------------- #
# prepare_static effect (svg is text/XML, introspectable)
# --------------------------------------------------------------------------- #
class TestPrepareStatic:
    def test_static_svg_reflects_log_axis_while_html_stays_linear(self, two_ref_bam, tmp_out):
        # run with both html and svg; the html figure keeps the interactive
        # (linear) axis, the static svg gets prepare_static (log axis, dtick=1)
        out = run_cli(two_ref_bam, tmp_out, ["--suffix", "html", "svg"])
        svg_path = Path(f"{out}_refA.svg")
        assert svg_path.exists()
        svg = svg_path.read_text()
        # a log axis with dtick=1 produces tick labels at 1, 10, 100, ...
        # (plotly renders these as text in the svg). At least "10" should appear
        # as a tick label for a coverage max around 8.8 (log10 upper = 1 => ticks 1,10).
        assert "<text" in svg  # svg has text elements
        # the interactive html figure should NOT force log on yaxis (linear default)
        data, layout = parse_plotly_payload(Path(f"{out}.html").read_text(), "plotly-refA")
        # yaxis type should not be "log" in the interactive figure
        assert layout.get("yaxis", {}).get("type") != "log"

    def test_static_svg_has_no_annotation_text(self, two_ref_bam, tmp_out):
        # prepare_static hides annotations; the stats title text should not render
        out = run_cli(two_ref_bam, tmp_out, ["--suffix", "svg"])
        svg = Path(f"{out}_refA.svg").read_text()
        # the stats annotation contains "<b>reference:</b>" — absent in static
        assert "<b>reference:</b>" not in svg
