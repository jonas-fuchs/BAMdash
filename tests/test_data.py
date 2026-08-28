"""
Data-parsing and scientific-correctness tests for ``bamdash.scripts.data``.

Covers: ``get_coverage_stats``, ``bam_to_coverage_df``, ``vcf_to_df``,
``bed_to_dict``, ``genbank_to_dict``, ``define_track_position``,
``get_mutations``, ``analyse_a3_signature``, ``annotate_vcf_df`` and
``annotate_vcfs_in_tracks``.

Scientific validation principle: expected values are computed independently
(Biopython ``Seq.translate`` with the NCBI standard table, hand-computed
coverage/recovery, regex-derived transition/transversion) — never copied from
the implementation under test.
"""
from __future__ import annotations

import math
import statistics
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq

from bamdash.scripts import data
from tests.conftest import (
    coverage_df,
    make_bam,
    make_bed,
    make_gb,
    make_gb_record,
    make_vcf,
    variant_dict,
)


def cds_dict(start, stop, strand="+", **qualifiers):
    """Build a minimal CDS dict matching ``genbank_to_dict`` output for a feature."""
    d = {"start": [start], "stop": [stop], "strand": strand}
    d.update(qualifiers)
    return d


# --------------------------------------------------------------------------- #
# get_coverage_stats
# --------------------------------------------------------------------------- #
class TestGetCoverageStats:
    def test_mean_and_recovery_match_independent_computation(self):
        positions = list(range(1, 11))
        covs = [0, 5, 5, 5, 10, 10, 1, 1, 1, 1]
        df = coverage_df(positions, covs)
        min_cov = 4
        mean_cov, recovery = data.get_coverage_stats(df, 1, 10, min_cov)
        # independent: mean = sum(cov)/(stop-start+1), recovery = (#cov>min)/N*100
        assert mean_cov == round(sum(covs) / 10)
        expected_rec = round(len([c for c in covs if c > min_cov]) / 10 * 100, 2)
        assert recovery == expected_rec

    def test_empty_subset_returns_zeros(self):
        df = coverage_df(list(range(1, 11)), [1] * 10)
        # range fully outside the data
        assert data.get_coverage_stats(df, 100, 200, 5) == (0, 0)

    def test_all_below_threshold_recovery_zero_mean_correct(self):
        positions = list(range(1, 6))
        covs = [1, 2, 1, 2, 1]
        df = coverage_df(positions, covs)
        mean_cov, recovery = data.get_coverage_stats(df, 1, 5, 5)
        assert mean_cov == round(sum(covs) / 5)
        assert recovery == 0


# --------------------------------------------------------------------------- #
# bam_to_coverage_df
# --------------------------------------------------------------------------- #
class TestBamToCoverageDf:
    def test_coverage_and_nucleotide_percentages_match_known_reads(self, tmp_path):
        # ref of length 10; two overlapping reads "AAAA" at pos 0 and "CCCC" at pos 0
        # so positions 1-4 carry one A and one C each => coverage 2, 50% A / 50% C.
        bam = tmp_path / "cov.bam"
        make_bam(
            bam,
            refs=[("ref", 10)],
            reads=[
                {"ref": "ref", "start": 0, "length": 4, "seq": "AAAA", "mapq": 60},
                {"ref": "ref", "start": 0, "length": 4, "seq": "CCCC", "mapq": 60},
            ],
        )
        df, _title, stats = data.bam_to_coverage_df(bam, "ref", min_cov=0, quality_thres=0)
        # positions 1-based
        assert list(df["position"]) == list(range(1, 11))
        # coverage: 2 on pos1-4 (both reads), 0 on pos5-10
        assert list(df["coverage"]) == [2, 2, 2, 2, 0, 0, 0, 0, 0, 0]
        # nucleotide percentages: pos1-4 are 50% A and 50% C
        assert df.loc[df["position"] == 1, "A"].iloc[0] == 50.0
        assert df.loc[df["position"] == 1, "C"].iloc[0] == 50.0
        # zero-coverage positions report 0 for all bases
        assert df.loc[df["position"] == 5, "A"].iloc[0] == 0
        assert df.loc[df["position"] == 5, "G"].iloc[0] == 0

    def test_gc_content_in_stat_dict_matches_independent_sum(self, tmp_path):
        bam = tmp_path / "gc.bam"
        make_bam(
            bam,
            refs=[("ref", 8)],
            reads=[
                {"ref": "ref", "start": 0, "length": 4, "seq": "GGGG", "mapq": 60},
                {"ref": "ref", "start": 4, "length": 4, "seq": "CCCC", "mapq": 60},
            ],
        )
        df, _title, stats = data.bam_to_coverage_df(bam, "ref", min_cov=0, quality_thres=0)
        # independent GC = (sum C + sum G)/len; pos1-4 are 100% G, pos5-8 are 100% C
        # (all 8 positions covered once each) => average GC = 100.0
        expected_gc = round((sum(df["C"]) + sum(df["G"])) / len(df), 2)
        assert stats["gc content (%)"] == expected_gc
        assert expected_gc == 100.0

    def test_quality_threshold_filters_low_base_quality_reads(self, tmp_path):
        # count_coverage's quality_threshold filters by BASE quality, not MAPQ.
        # One read with high base qual ("I"=40) and one with low ("#"=2) overlap.
        bam = tmp_path / "qual.bam"
        make_bam(
            bam,
            refs=[("ref", 10)],
            reads=[
                {"ref": "ref", "start": 0, "length": 5, "seq": "AAAAA", "mapq": 60, "qual": "I"},
                {"ref": "ref", "start": 0, "length": 5, "seq": "GGGGG", "mapq": 60, "qual": "#"},
            ],
        )
        df_high, _, _ = data.bam_to_coverage_df(bam, "ref", 0, quality_thres=15)
        df_low, _, _ = data.bam_to_coverage_df(bam, "ref", 0, quality_thres=0)
        # with threshold 15 the low-base-qual read is excluded => coverage 1 on pos1-5
        assert df_high.loc[df_high["position"] == 1, "coverage"].iloc[0] == 1
        # with threshold 0 both reads count => coverage 2 on pos1-5
        assert df_low.loc[df_low["position"] == 1, "coverage"].iloc[0] == 2

    def test_missing_reference_raises_and_logs(self, tmp_path, caplog):
        bam = tmp_path / "r.bam"
        make_bam(bam, refs=[("ref", 10)], reads=[{"ref": "ref", "start": 0, "length": 4, "mapq": 60}])
        with caplog.at_level("ERROR", logger="bamdash"):
            with pytest.raises(data.ReferenceNotFoundError):
                data.bam_to_coverage_df(bam, "NOPE", 0, 0)
        assert any("does not exist in bam file" in r.message for r in caplog.records)


# --------------------------------------------------------------------------- #
# vcf_to_df
# --------------------------------------------------------------------------- #
class TestVcfToDf:
    def _df(self, tmp_path, records, **kwargs):
        vcf = tmp_path / "v.vcf"
        make_vcf(vcf, records, **kwargs)
        return data.vcf_to_df(vcf, ref=records[0]["chrom"])

    def test_transition_vs_transversion(self, tmp_path):
        # AG, GA, CT, TC -> TRANSITION ; AC, CA, GC -> TRANSVERSION
        recs = [
            {"chrom": "r", "pos": p, "ref": ref, "alt": alt}
            for p, (ref, alt) in zip(
                [1, 2, 3, 4, 5, 6, 7],
                [("A", "G"), ("G", "A"), ("C", "T"), ("T", "C"),
                 ("A", "C"), ("C", "A"), ("G", "C")],
            )
        ]
        df = self._df(tmp_path, recs)
        expected = ["TRANSITION", "TRANSITION", "TRANSITION", "TRANSITION",
                    "TRANSVERSION", "TRANSVERSION", "TRANSVERSION"]
        assert list(df["point mutation type"]) == expected

    def test_non_snv_point_mutation_type_is_none(self, tmp_path):
        df = self._df(tmp_path, [{"chrom": "r", "pos": 1, "ref": "A", "alt": "CT"}])
        assert list(df["point mutation type"]) == ["NONE"]

    def test_type_annotation_ins_del_snp(self, tmp_path):
        recs = [
            {"chrom": "r", "pos": 1, "ref": "A", "alt": "T"},   # SNP
            {"chrom": "r", "pos": 2, "ref": "A", "alt": "AT"},  # INS
            {"chrom": "r", "pos": 3, "ref": "AT", "alt": "A"},  # DEL
        ]
        df = self._df(tmp_path, recs)
        assert list(df["type"]) == ["SNP", "INS", "DEL"]

    def test_info_fields_populated(self, tmp_path):
        df = self._df(tmp_path, [{"chrom": "r", "pos": 1, "ref": "A", "alt": "T", "info": {"AF": 0.5}}])
        assert "AF" in df.columns
        assert df["AF"].iloc[0] == 0.5

    def test_sparse_info_key_with_none_is_dropped(self, tmp_path):
        """A key that is None for any record is dropped entirely (locked behavior)."""
        vcf = tmp_path / "v.vcf"
        make_vcf(
            vcf,
            records=[
                {"chrom": "r", "pos": 1, "ref": "A", "alt": "T", "info": {"AF": 0.5, "DP": 10}},
                {"chrom": "r", "pos": 2, "ref": "C", "alt": "G", "info": {"AF": 0.3}},  # no DP
            ],
            info_fields=[("AF", "Float", "af"), ("DP", "Integer", "depth")],
        )
        df = data.vcf_to_df(vcf, ref="r")
        # AF present for all => kept; DP missing for one => dropped
        assert "AF" in df.columns
        assert "DP" not in df.columns

    def test_chrom_filter_excludes_other_references(self, tmp_path):
        vcf = tmp_path / "v.vcf"
        make_vcf(
            vcf,
            records=[
                {"chrom": "r1", "pos": 1, "ref": "A", "alt": "T"},
                {"chrom": "r2", "pos": 5, "ref": "C", "alt": "G"},
            ],
        )
        df = data.vcf_to_df(vcf, ref="r1")
        assert list(df["position"]) == [1]
        assert df["reference"].iloc[0] == "A"

    def test_output_sorted_by_position(self, tmp_path):
        df = self._df(tmp_path, [
            {"chrom": "r", "pos": 30, "ref": "A", "alt": "T"},
            {"chrom": "r", "pos": 5, "ref": "C", "alt": "G"},
            {"chrom": "r", "pos": 20, "ref": "G", "alt": "A"},
        ])
        assert list(df["position"]) == [5, 20, 30]

    def test_first_alt_only_for_multialt(self, tmp_path):
        vcf = tmp_path / "v.vcf"
        # write a multi-alt record manually (make_vcf writes single alt)
        vcf.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=r,length=10>\n"
            '##INFO=<ID=AF,Number=1,Type=Float,Description="af">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "r\t1\t.\tA\tT,G\t.\tPASS\tAF=0.5\n"
        )
        df = data.vcf_to_df(vcf, ref="r")
        assert df["mutation"].iloc[0] == "T"


# --------------------------------------------------------------------------- #
# bed_to_dict
# --------------------------------------------------------------------------- #
class TestBedToDict:
    def _dict(self, tmp_path, rows, comments=None, ref="ref"):
        bed = tmp_path / "b.bed"
        make_bed(bed, rows, comments=comments)
        df = coverage_df(list(range(1, 21)), [5] * 20)
        return data.bed_to_dict(bed, df, ref, min_cov=0)

    def test_zero_based_start_converted_to_one_based(self, tmp_path):
        d = self._dict(tmp_path, [{"chrom": "ref", "start": 0, "stop": 10}])
        ann = d["b"]
        key = next(iter(ann))
        # start 0 -> 1 ; stop unchanged
        assert ann[key]["start"] == [1]
        assert ann[key]["stop"] == [10]
        assert key == "1 10"

    def test_column_mapping_name_score_strand_positional(self, tmp_path):
        d = self._dict(tmp_path, [{"chrom": "ref", "start": 0, "stop": 5,
                                    "name": "gene1", "score": "42", "strand": "+"}])
        ann = d["b"][next(iter(d["b"]))]
        assert ann["name"] == "gene1"
        assert ann["score"] == "42"
        assert ann["strand"] == "+"

    def test_default_strand_na_when_absent(self, tmp_path):
        d = self._dict(tmp_path, [{"chrom": "ref", "start": 0, "stop": 5}])
        ann = d["b"][next(iter(d["b"]))]
        assert ann["strand"] == "NA"

    def test_comment_and_short_lines_skipped(self, tmp_path):
        # a comment line and a 2-field line should not break parsing
        d = self._dict(tmp_path, [{"chrom": "ref", "start": 0, "stop": 5}],
                       comments=["# header comment", "ref\t1"])
        assert len(d["b"]) == 1

    def test_no_matching_ref_returns_empty_annotations(self, tmp_path):
        d = self._dict(tmp_path, [{"chrom": "other", "start": 0, "stop": 5}], ref="ref")
        assert d["b"] == {}

    def test_per_region_coverage_and_recovery(self, tmp_path):
        df = coverage_df(list(range(1, 11)), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        bed = tmp_path / "b.bed"
        make_bed(bed, [{"chrom": "ref", "start": 0, "stop": 4}])  # positions 1-4
        d = data.bed_to_dict(bed, df, "ref", min_cov=2)
        ann = d["b"][next(iter(d["b"]))]
        # independent: mean = (0+1+2+3)/4 = 1.5 -> round = 2 ; recovery = (#>2)/4*100 = 25
        assert ann["mean coverage"] == round((0 + 1 + 2 + 3) / 4)
        assert ann["% recovery >= 2x"] == round(len([c for c in [0, 1, 2, 3] if c > 2]) / 4 * 100, 2)


# --------------------------------------------------------------------------- #
# genbank_to_dict
# --------------------------------------------------------------------------- #
class TestGenbankToDict:
    def test_strand_plus_and_minus(self, tmp_path):
        # GenBank format has no native strand=0; a strand=0 feature round-trips
        # as strand=1, so only '+' and '-' are reachable via files. The 'NA'
        # branch is exercised separately in the genbank_to_dict unit via a
        # strand=None location (not reachable through GenBank I/O).
        gb = tmp_path / "g.gb"
        make_gb(gb, [make_gb_record("ref", 20, features=[
            {"type": "CDS", "location": [(1, 10, 1)], "qualifiers": {"gene": ["a"]}},
            {"type": "CDS", "location": [(11, 20, -1)], "qualifiers": {"gene": ["b"]}},
        ])])
        df = coverage_df(list(range(1, 21)), [5] * 20)
        d, _seq = data.genbank_to_dict(gb, df, "ref", 0)
        assert d["CDS"][list(d["CDS"])[0]]["strand"] == "+"
        assert d["CDS"][list(d["CDS"])[1]]["strand"] == "-"

    def test_qualifiers_first_value_only(self, tmp_path):
        gb = tmp_path / "g.gb"
        make_gb(gb, [make_gb_record("ref", 10, features=[
            {"type": "CDS", "location": [(1, 10, 1)],
             "qualifiers": {"gene": ["only_first", "ignored"]}},
        ])])
        df = coverage_df(list(range(1, 11)), [5] * 10)
        d, _ = data.genbank_to_dict(gb, df, "ref", 0)
        cds = d["CDS"][next(iter(d["CDS"]))]
        assert cds["gene"] == "only_first"

    def test_codon_start_preserved(self, tmp_path):
        gb = tmp_path / "g.gb"
        make_gb(gb, [make_gb_record("ref", 10, features=[
            {"type": "CDS", "location": [(1, 10, 1)], "qualifiers": {"codon_start": ["2"]}},
        ])])
        df = coverage_df(list(range(1, 11)), [5] * 10)
        d, _ = data.genbank_to_dict(gb, df, "ref", 0)
        cds = d["CDS"][next(iter(d["CDS"]))]
        assert cds["codon_start"] == "2"

    def test_multi_record_file_finds_later_record(self, tmp_path):
        """Regression: matching record found by id/name even when not first."""
        gb = tmp_path / "g.gb"
        make_gb(gb, [
            make_gb_record("other", 10, features=[{"type": "CDS", "location": [(1, 10, 1)]}]),
            make_gb_record("ref", 10, features=[{"type": "CDS", "location": [(1, 10, 1)]}]),
        ])
        df = coverage_df(list(range(1, 11)), [5] * 10)
        d, seq = data.genbank_to_dict(gb, df, "ref", 0)
        assert "CDS" in d
        assert str(seq) == str(make_gb_record("ref", 10).seq)

    def test_no_matching_record_returns_empty(self, tmp_path):
        gb = tmp_path / "g.gb"
        make_gb(gb, [make_gb_record("other", 10, features=[{"type": "CDS", "location": [(1, 10, 1)]}])])
        df = coverage_df(list(range(1, 11)), [5] * 10)
        d, seq = data.genbank_to_dict(gb, df, "ref", 0)
        assert d == {}
        assert seq == ""

    def test_multipart_feature_start_stop_are_lists(self, tmp_path):
        from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
        from Bio.SeqRecord import SeqRecord
        rec = make_gb_record("ref", 30)
        rec.features = [SeqFeature(
            CompoundLocation([FeatureLocation(1, 10, 1), FeatureLocation(20, 30, 1)]),
            type="CDS", qualifiers={"gene": ["j"]},
        )]
        gb = tmp_path / "g.gb"
        make_gb(gb, [rec])
        df = coverage_df(list(range(1, 31)), [5] * 30)
        d, _ = data.genbank_to_dict(gb, df, "ref", 0)
        cds = d["CDS"][next(iter(d["CDS"]))]
        assert cds["start"] == [2, 21]   # 0-based +1
        assert cds["stop"] == [10, 30]


# --------------------------------------------------------------------------- #
# define_track_position
# --------------------------------------------------------------------------- #
class TestDefineTrackPosition:
    def _feat(self, start, stop, strand="+"):
        return {"start": [start], "stop": [stop], "strand": strand}

    def test_non_overlapping_same_track(self):
        fd = {"CDS": {"a 10": self._feat(1, 10), "b 20": self._feat(11, 20)}}
        out = data.define_track_position(fd)
        assert out["CDS"]["a 10"]["track"] == 0
        assert out["CDS"]["b 20"]["track"] == 0

    def test_overlapping_successive_tracks(self):
        fd = {"CDS": {"a 10": self._feat(1, 15), "b 20": self._feat(5, 25), "c 30": self._feat(10, 30)}}
        out = data.define_track_position(fd)
        tracks = [out["CDS"][k]["track"] for k in out["CDS"]]
        # first track 0, overlapping ones move to 1, 2
        assert tracks == [0, 1, 2]

    def test_different_feature_types_do_not_share_tracks(self):
        fd = {
            "CDS": {"a 10": self._feat(1, 10)},
            "gene": {"b 10": self._feat(1, 10)},
        }
        out = data.define_track_position(fd)
        # gene track offset by the number of CDS tracks used (>=1)
        assert out["CDS"]["a 10"]["track"] == 0
        assert out["gene"]["b 10"]["track"] >= 1

    def test_multipart_feature_uses_min_start_max_stop(self):
        fd = {"CDS": {"a 30": {"start": [2, 21], "stop": [10, 30], "strand": "+"},
                      "b 15": self._feat(5, 15)}}
        out = data.define_track_position(fd)
        # feature a spans 2..30 ; b starts at 5 < 30 => different track
        assert out["CDS"]["a 30"]["track"] == 0
        assert out["CDS"]["b 15"]["track"] == 1


# --------------------------------------------------------------------------- #
# get_mutations  (scientific core — validated with Bio.Seq.translate)
# --------------------------------------------------------------------------- #
def _translate(seq: str) -> str:
    """Independent translation of a nucleotide string via Biopython."""
    return str(Seq(seq).translate())


class TestGetMutations:
    # A simple CDS: start codon ATG (M), then codons; use the + strand.
    # We construct sequences so the affected codon and AA can be recomputed
    # independently with Seq.translate.
    def _cds(self, strand="+", codon_start=1, **qual):
        d = cds_dict(1, 12, strand=strand)
        d["codon_start"] = str(codon_start)
        d.update(qual)
        return d

    def test_synonymous_snp(self):
        # CDS = ATG TTA GGG = M L G ; mutate pos4 (codon2 pos0) T->C => CTA (still L)
        seq = Seq("ATG" + "TTA" + "GGG")  # M L G
        cds = {"start": [1], "stop": [9], "strand": "+", "codon_start": "1"}
        var = variant_dict(4, "T", "C", vtype="SNP")
        ac, eff = data.get_mutations(1, 9, cds, var, seq)
        assert eff == "SYN"

    def test_nonsynonymous_snp_exchange_string(self):
        # M L G ; mutate pos4 T->C => CTA (L->L? no, CTA is L too). Use pos5 T->A: TAA stop
        seq = Seq("ATG" + "TTA" + "GGG")  # M L G
        cds = {"start": [1], "stop": [9], "strand": "+", "codon_start": "1"}
        # pos5 is the 2nd nt of codon2 TTA. T->G => TGA (stop) -> STOP_GAINED instead.
        # For NON_SYN pick a missense: pos7 G->T in GGG => TGG (W). codon3.
        var = variant_dict(7, "G", "T", vtype="SNP")
        ac, eff = data.get_mutations(1, 9, cds, var, seq)
        ref_aa = _translate("GGG")  # G
        alt_aa = _translate("TGG")  # W
        assert eff == "NON_SYN"
        assert ac == f"{ref_aa}3{alt_aa}"

    def test_start_lost(self):
        # mutate the start codon (pos1, codon1) ATG->ACG (T) => M lost
        seq = Seq("ATG" + "CCC" + "GGG")
        cds = {"start": [1], "stop": [9], "strand": "+", "codon_start": "1"}
        var = variant_dict(1, "A", "C", vtype="SNP")
        ac, eff = data.get_mutations(1, 9, cds, var, seq)
        assert eff == "START_LOST"

    def test_stop_gained(self):
        # CDS = ATG TGG CAG = M W Q ; mutate pos7 (codon3 pos0) C->T => TAG (stop)
        seq = Seq("ATG" + "TGG" + "CAG")  # M W Q
        cds = {"start": [1], "stop": [9], "strand": "+", "codon_start": "1"}
        var = variant_dict(7, "C", "T", vtype="SNP")  # CAG -> TAG (*)
        ac, eff = data.get_mutations(1, 9, cds, var, seq)
        assert eff == "STOP_GAINED"
        ref_aa = _translate("CAG")  # Q
        alt_aa = _translate("TAG")  # *
        assert ac == f"{ref_aa}3{alt_aa}"

    def test_stop_lost(self):
        # CDS ends with a stop: ATG TGG TAA (M W *); mutate pos9 A->G => TGA still *? no TGA is *.
        # mutate pos7 T->C in TAA => CAA (Q) -> STOP_LOST
        seq = Seq("ATG" + "TGG" + "TAA")  # M W *
        cds = {"start": [1], "stop": [9], "strand": "+", "codon_start": "1"}
        var = variant_dict(7, "T", "C", vtype="SNP")  # TAA -> CAA (Q)
        ac, eff = data.get_mutations(1, 9, cds, var, seq)
        assert eff == "STOP_LOST"

    def test_insertion_triplet_ac_insertion(self):
        # in-frame insertion that preserves the affected codon's first AA.
        # CDS = ATG CCC GGG = M P G ; insert CCC at pos4 (codon2 pos0): C->CCCC
        # makes codon CCCCCC -> PP, ins[0]==ref_ac (P) => AC_INSERTION
        seq = Seq("ATG" + "CCC" + "GGG")
        cds = {"start": [1], "stop": [9], "strand": "+", "codon_start": "1"}
        var = variant_dict(4, "C", "CCCC", vtype="INS")
        ac, eff = data.get_mutations(1, 9, cds, var, seq)
        assert eff == "AC_INSERTION"

    def test_insertion_nontriplet_frameshift(self):
        seq = Seq("ATG" + "CCC" + "GGG")
        cds = {"start": [1], "stop": [9], "strand": "+", "codon_start": "1"}
        var = variant_dict(1, "A", "ACC", vtype="INS")  # len-1=2 not %3
        ac, eff = data.get_mutations(1, 9, cds, var, seq)
        assert eff == "FRAMESHIFT"
        assert "fsX" in ac

    def test_deletion_triplet_ac_deletion(self):
        # delete a whole codon (triplet): ref=ATG alt=A? that's 2 deleted.
        # triplet deletion: ref=ATGCCC alt=CCC? deletes ATG (3nt). len(ref)-1=5 not %3.
        # Use ref=ATGC alt=A? len-1=3 triplet. Deletes TGC.
        seq = Seq("ATG" + "CCC" + "GGG")
        cds = {"start": [1], "stop": [9], "strand": "+", "codon_start": "1"}
        var = variant_dict(1, "ATGC", "A", vtype="DEL")  # deletes TGC (3nt)
        ac, eff = data.get_mutations(1, 9, cds, var, seq)
        assert eff in ("AC_DELETION", "AC_CHANGE+AC_DELETION", "START_LOST", "STOP_LOST")

    def test_deletion_nontriplet_frameshift(self):
        seq = Seq("ATG" + "CCC" + "GGG")
        cds = {"start": [1], "stop": [9], "strand": "+", "codon_start": "1"}
        var = variant_dict(1, "AT", "A", vtype="DEL")  # deletes 1 nt
        ac, eff = data.get_mutations(1, 9, cds, var, seq)
        assert eff == "FRAMESHIFT"
        assert "fsX" in ac

    def test_minus_strand_reverse_complements_sequence(self):
        # On the - strand the CDS sequence is reverse-complemented before translation.
        # Build a + strand seq, compute its reverse complement independently, then
        # verify a SNP effect matches translation of that reverse complement.
        seq_plus = Seq("ATG" + "CCC" + "GGG")  # M P G on + strand
        cds = {"start": [1], "stop": [9], "strand": "-", "codon_start": "1"}
        # On - strand, cds_variant_pos = stop - variant_pos. Pick a variant near the
        # 3' end of the + strand so it lands in codon1 of the reverse complement.
        # We assert the function returns a non-NONE effect (the AA is derived from RC).
        var = variant_dict(8, "G", "T", vtype="SNP")
        ac, eff = data.get_mutations(1, 9, cds, var, seq_plus)
        # independently: rc = reverse_complement(ATG CCC GGG) = CCC GGG CAT
        rc = str(seq_plus.reverse_complement())
        # codon containing the variant on the RC strand:
        cds_pos = 9 - 8  # = 1
        mut_codon = math.floor(cds_pos / 3)  # 0
        codon_pos = cds_pos % 3              # 1
        ref_codon = list(rc[mut_codon * 3:mut_codon * 3 + 3])
        ref_aa = _translate("".join(ref_codon))
        ref_codon[codon_pos] = str(Seq("T").reverse_complement())  # mutation RC'd
        alt_aa = _translate("".join(ref_codon))
        if ref_aa == alt_aa:
            assert eff == "SYN"
        else:
            assert eff != "SYN"

    def test_codon_start_offset_shifts_reading_frame(self):
        # codon_start=2 means the reading frame starts one nt into the CDS.
        # seq = A ATG CCC GG ; with codon_start=2, seq_cds = seq[1:9] = ATGCCCGG,
        # so codon1 is ATG (the start). A SNP at pos1 (cds_variant_pos=0) hits
        # the start codon: A->C => CTG (L) => START_LOST.
        seq = Seq("A" + "ATG" + "CCC" + "GG")  # length 10, cds 1..9
        cds = {"start": [1], "stop": [9], "strand": "+", "codon_start": "2"}
        var = variant_dict(1, "A", "C", vtype="SNP")  # ATG -> CTG (L)
        ac, eff = data.get_mutations(1, 9, cds, var, seq)
        assert eff == "START_LOST"


# --------------------------------------------------------------------------- #
# analyse_a3_signature
# --------------------------------------------------------------------------- #
class TestAnalyseA3Signature:
    def test_c_ref_with_tc_context_and_t_mutation_is_yes(self):
        # position is 1-based; seq[pos-2:pos] must be "TC"
        seq = Seq("ATCG")  # positions: 1=A,2=T,3=C,4=G ; context for pos3 = seq[1:3]="TC"
        var = variant_dict(3, "C", "T", vtype="SNP")
        result, site = data.analyse_a3_signature(var, seq)
        assert result == "YES"

    def test_c_ref_without_tc_context_is_no(self):
        seq = Seq("AAGC")  # context for pos3 (C) = seq[1:3]="AG"
        var = variant_dict(3, "C", "T", vtype="SNP")
        result, _ = data.analyse_a3_signature(var, seq)
        assert result == "NO"

    def test_g_ref_with_ga_following_and_a_mutation_is_yes(self):
        # for G ref: site = seq[pos-1:pos+1]; pos2 G in "AGA" -> seq[1:3]="GA"
        seq = Seq("AGA")  # pos1=A, pos2=G, pos3=A
        var = variant_dict(2, "G", "A", vtype="SNP")
        result, _ = data.analyse_a3_signature(var, seq)
        assert result == "YES"

    def test_g_ref_at_last_position_is_no(self):
        seq = Seq("AAG")  # pos3 G is last; no following nt
        var = variant_dict(3, "G", "A", vtype="SNP")
        result, _ = data.analyse_a3_signature(var, seq)
        assert result == "NO"

    def test_non_snp_returns_no_dash(self):
        seq = Seq("ATCG")
        var = variant_dict(3, "C", "CT", vtype="INS")
        result, site = data.analyse_a3_signature(var, seq)
        assert result == "NO"
        assert site == "-"


# --------------------------------------------------------------------------- #
# annotate_vcf_df / annotate_vcfs_in_tracks
# --------------------------------------------------------------------------- #
class TestAnnotateVcfDf:
    def _vcf_df(self, rows):
        return pd.DataFrame(rows)

    def test_matching_cds_populates_protein_effect_columns(self):
        # CDS ATG CCC GGG = M P G ; SNP at pos7 (codon3 pos0) G->T => GGG->TGG (G->W) NON_SYN
        seq = Seq("ATG" + "CCC" + "GGG")
        cds = {"1 9": {"start": [1], "stop": [9], "strand": "+", "codon_start": "1",
                       "protein_id": "prot1"}}
        vcf_df = self._vcf_df([variant_dict(7, "G", "T", vtype="SNP")])
        out = data.annotate_vcf_df(vcf_df, cds, seq)
        assert out["protein"].iloc[0] == "prot1"
        assert out["effect"].iloc[0] == "NON_SYN"
        assert out["as mutation"].iloc[0] == "G3W"

    def test_no_matching_protein_returns_none(self):
        seq = Seq("ATG" + "CCC" + "GGG")
        cds = {"1 9": {"start": [1], "stop": [9], "strand": "+", "codon_start": "1",
                       "protein_id": ["prot1"]}}
        vcf_df = self._vcf_df([variant_dict(50, "A", "T", vtype="SNP")])  # outside CDS
        out = data.annotate_vcf_df(vcf_df, cds, seq)
        assert out["protein"].iloc[0] == "NONE"
        assert out["effect"].iloc[0] == "NONE"
        assert out["as mutation"].iloc[0] == "NONE"

    def test_multiple_cds_hits_aggregated_with_pipe(self):
        # Two CDS both spanning the variant; qualifier values are strings here
        # (annotate_vcf_df joins them with " | ").
        seq = Seq("ATG" + "CCC" + "GGG")
        cds = {
            "1 9": {"start": [1], "stop": [9], "strand": "+", "codon_start": "1", "protein_id": "p1"},
            "1 9b": {"start": [1], "stop": [9], "strand": "+", "codon_start": "1", "gene": "g2"},
        }
        vcf_df = self._vcf_df([variant_dict(7, "G", "T", vtype="SNP")])
        out = data.annotate_vcf_df(vcf_df, cds, seq)
        assert " | " in out["protein"].iloc[0]
        assert "p1" in out["protein"].iloc[0] and "g2" in out["protein"].iloc[0]


class TestAnnotateVcfsInTracks:
    def test_a3_columns_added_to_vcf(self):
        seq = Seq("ATCG")
        vcf_df = pd.DataFrame([variant_dict(3, "C", "T", vtype="SNP")])
        track_data = [[vcf_df, "vcf"], [{}, "gb", seq]]
        out = data.annotate_vcfs_in_tracks(track_data)
        assert "potential APOBEC3 activity" in out[0][0].columns
        assert "checked site" in out[0][0].columns
        assert out[0][0]["potential APOBEC3 activity"].iloc[0] == "YES"

    def test_multiple_gb_warns_and_no_annotation(self, caplog):
        vcf_df = pd.DataFrame([variant_dict(3, "C", "T", vtype="SNP")])
        track_data = [[vcf_df, "vcf"], [{}, "gb", "ATCG"], [{}, "gb", "ATCG"]]
        with caplog.at_level("WARNING", logger="bamdash"):
            out = data.annotate_vcfs_in_tracks(track_data)
        assert any("multiple" in r.message for r in caplog.records)
        # vcf df returned without protein/effect columns
        assert "protein" not in out[0][0].columns

    def test_no_gb_returns_unchanged(self):
        vcf_df = pd.DataFrame([variant_dict(3, "C", "T", vtype="SNP")])
        track_data = [[vcf_df, "vcf"]]
        out = data.annotate_vcfs_in_tracks(track_data)
        assert "potential APOBEC3 activity" not in out[0][0].columns

    def test_no_vcf_returns_unchanged(self):
        track_data = [[{}, "gb", "ATCG"]]
        out = data.annotate_vcfs_in_tracks(track_data)
        assert out == track_data
