"""
Shared pytest fixtures and synthetic-data builders for the BAMdash test suite.

The suite is self-contained: all BAM/VCF/BED/GenBank fixtures are generated
programmatically here, so the tests do not depend on the shipped ``test/``
sample data. Builders are parameterisable so individual tests can construct
minimal inputs for the behaviour under test.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
import pysam
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqIO import write as seqio_write
from Bio.SeqRecord import SeqRecord


# --------------------------------------------------------------------------- #
# Synthetic file builders
# --------------------------------------------------------------------------- #
def make_bam(path: Path, refs: list[tuple[str, int]], reads: list[dict]) -> None:
    """Create and index a synthetic BAM.

    :param path: output ``.bam`` path (a ``.bai`` is created next to it)
    :param refs: list of ``(reference_name, length_bp)`` header entries
    :param reads: list of read specs, each a dict with keys:
        ``ref`` (name), ``start`` (0-based), ``length`` (int), ``mapq`` (int),
        ``seq`` (optional query sequence; defaults to ``"ACGT"`` repeats),
        ``qual`` (optional quality char, default ``"I"``). All reads are
        forward, perfectly matched (single CMATCH cigar).
    """
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": name, "LN": ln} for name, ln in refs]}
    ref_index = {name: i for i, (name, _) in enumerate(refs)}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for i, spec in enumerate(reads):
            length = spec["length"]
            seq = spec.get("seq", "ACGT" * (length // 4 + 1))[:length]
            r = pysam.AlignedSegment()
            r.query_name = spec.get("name", f"read{i}")
            r.query_sequence = seq
            r.flag = 0
            r.reference_id = ref_index[spec["ref"]]
            r.reference_start = spec["start"]
            r.mapping_quality = spec.get("mapq", 60)
            r.cigar = [(pysam.CMATCH, length)]
            r.query_qualities = pysam.qualitystring_to_array(spec.get("qual", "I") * length)
            bam.write(r)
    pysam.index(str(path))


def make_vcf(path: Path, records: list[dict], contigs: list[tuple[str, int]] | None = None,
             info_fields: list[tuple[str, str, str]] | None = None) -> None:
    """Create a synthetic VCFv4.2 file.

    :param records: list of dicts with keys ``chrom``, ``pos`` (1-based),
        ``ref``, ``alt``, ``info`` (dict of INFO key->value, optional).
    :param contigs: list of ``(id, length)``; inferred from records if None.
    :param info_fields: list of ``(ID, Type, Description)`` header INFO lines;
        inferred from records if None.
    """
    if contigs is None:
        lengths = {}
        for rec in records:
            lengths[rec["chrom"]] = max(lengths.get(rec["chrom"], 0), rec["pos"] + 10)
        contigs = [(c, lengths[c]) for c in lengths]
    if info_fields is None:
        seen = []
        for rec in records:
            for k in rec.get("info", {}):
                if k not in seen:
                    seen.append(k)
        info_fields = [(k, "Float", f"{k} description") for k in seen]
    lines = ["##fileformat=VCFv4.2"]
    for cid, ln in contigs:
        lines.append(f"##contig=<ID={cid},length={ln}>")
    for iid, itype, idesc in info_fields:
        lines.append(f'##INFO=<ID={iid},Number=1,Type={itype},Description="{idesc}">')
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    for rec in records:
        info_str = ";".join(f"{k}={v}" for k, v in rec.get("info", {}).items()) or "."
        lines.append(
            f"{rec['chrom']}\t{rec['pos']}\t.\t{rec['ref']}\t{rec['alt']}\t.\tPASS\t{info_str}"
        )
    path.write_text("\n".join(lines) + "\n")


def make_bed(path: Path, rows: list[dict], comments: list[str] | None = None) -> None:
    """Create a synthetic BED file (tab-delimited).

    :param rows: list of dicts with keys ``chrom``, ``start`` (0-based),
        ``stop``, and optional ``name``, ``score``, ``strand``, plus any extra
        columns under ``extra`` (list of strings appended after strand).
    :param comments: lines starting with ``#`` to insert at the top.
    """
    lines = list(comments or [])
    for row in rows:
        fields = [row["chrom"], str(row["start"]), str(row["stop"])]
        if "name" in row:
            fields.append(str(row["name"]))
        if "score" in row:
            fields.append(str(row["score"]))
        if "strand" in row:
            fields.append(str(row["strand"]))
        for extra in row.get("extra", []):
            fields.append(str(extra))
        lines.append("\t".join(fields))
    path.write_text("\n".join(lines) + "\n")


def make_gb_record(ref_id: str, length: int, features: list[dict] | None = None,
                   seq_str: str | None = None) -> SeqRecord:
    """Build a minimal GenBank SeqRecord.

    :param features: list of feature specs, each a dict with keys ``type``
        (e.g. ``"CDS"``), ``location`` (a list of ``(start, stop, strand)``
        tuples; a single tuple makes a simple FeatureLocation, multiple make a
        CompoundLocation/join), and ``qualifiers`` (dict of key->list[str]).
    """
    seq = Seq(seq_str or ("ATG" * (length // 3) + "A" * (length % 3)))[:length]
    rec = SeqRecord(seq, id=ref_id, name=ref_id, annotations={"molecule_type": "DNA"})
    rec.features = []
    for spec in features or []:
        parts = spec["location"]
        if len(parts) == 1:
            s, e, strand = parts[0]
            loc = FeatureLocation(s, e, strand=strand)
        else:
            loc = CompoundLocation([FeatureLocation(s, e, strand=strand) for s, e, strand in parts])
        feat = SeqFeature(loc, type=spec["type"])
        feat.qualifiers = {k: list(v) for k, v in spec.get("qualifiers", {}).items()}
        rec.features.append(feat)
    return rec


def make_gb(path: Path, records: list[SeqRecord]) -> None:
    """Write one or more SeqRecords to a GenBank file."""
    seqio_write(records, str(path), "genbank")


# --------------------------------------------------------------------------- #
# In-memory data builders (no file I/O)
# --------------------------------------------------------------------------- #
def coverage_df(positions: list[int], coverages: list[int],
                bases: list[tuple[int, int, int, int]] | None = None) -> pd.DataFrame:
    """Build an in-memory coverage DataFrame matching ``bam_to_coverage_df`` output.

    :param positions: 1-based positions
    :param coverages: per-position total coverage
    :param bases: optional per-position ``(A, C, G, T)`` counts; when None an
        even split is assumed (remainder to A). Percentages are rounded to 1dp.
    """
    rows = []
    for i, (pos, cov) in enumerate(zip(positions, coverages)):
        if bases is not None:
            a, c, g, t = bases[i]
        elif cov == 0:
            a = c = g = t = 0
        else:
            a = cov // 4
            c = cov // 4
            g = cov // 4
            t = cov - a - c - g
        pct = lambda n: round(n / cov * 100, 1) if cov else 0
        rows.append((pos, cov, pct(a), pct(c), pct(g), pct(t)))
    return pd.DataFrame(rows, columns=["position", "coverage", "A", "C", "G", "T"])


def variant_dict(position: int, ref: str, alt: str, vtype: str | None = None,
                 **extra) -> dict:
    """Build a minimal variant dict matching a row of a ``vcf_to_df`` DataFrame."""
    if vtype is None:
        if len(alt) > len(ref):
            vtype = "INS"
        elif len(alt) < len(ref):
            vtype = "DEL"
        else:
            vtype = "SNP"
    d = {"position": position, "reference": ref, "mutation": alt, "type": vtype}
    d.update(extra)
    return d


def cds_dict(start: int, stop: int, strand: str = "+", **qualifiers) -> dict:
    """Build a minimal CDS dict matching ``genbank_to_dict`` output for a feature."""
    d = {"start": [start], "stop": [stop], "strand": strand}
    d.update(qualifiers)
    return d


# --------------------------------------------------------------------------- #
# CLI + HTML helpers
# --------------------------------------------------------------------------- #
def run_cli(bam: Path, prefix: Path, extra_args: list[str] | None = None) -> Path:
    """Invoke ``bamdash.command.main`` in-process and return the output prefix."""
    from bamdash.command import main
    args = ["-b", str(bam), "-p", str(prefix)] + list(extra_args or [])
    main(args)
    return prefix


def read_html(prefix: Path) -> str:
    return Path(f"{prefix}.html").read_text(encoding="utf-8")


def parse_plotly_payload(html_text: str, div_id: str) -> tuple[list, dict]:
    """Extract the ``(data, layout)`` JSON from a ``Plotly.newPlot`` call.

    The plotly ``to_html`` snippet spans multiple lines, e.g.::

        Plotly.newPlot(
            "plotly-refA",
            [ {...}, ... ],
            { ...layout... },
            {config}
        );

    We locate the ``"div_id"`` string that is the first argument, then
    ``raw_decode`` the JSON data array and layout object that follow it.
    ``json.JSONDecoder.raw_decode`` handles nested arrays/objects without
    regex bracket matching.

    Numeric arrays are serialized by plotly as base64-encoded binary dicts
    (``{"dtype": "f8", "bdata": "..."}``); ``_decode_plotly_arrays`` walks the
    parsed payload and converts those back to plain Python lists so tests can
    assert on values directly.

    :returns: ``(data_list, layout_dict)`` with all plotly binary arrays decoded
    """
    import json

    # find the div_id string literal that begins the newPlot args
    id_marker = f'"{div_id}"'
    start = html_text.find(id_marker)
    assert start != -1, f"div id '{div_id}' not found in HTML"
    # the data array starts at the first '[' after the id marker
    decoder = json.JSONDecoder()
    idx = start + len(id_marker)
    # advance to the first '[' (start of the data array), skipping commas/whitespace
    while html_text[idx] != "[":
        idx += 1
        assert idx < len(html_text), "data array not found after div id"
    data, data_end = decoder.raw_decode(html_text, idx)
    # advance to the layout object (first '{'), skipping whitespace/commas
    j = data_end
    while html_text[j] != "{":
        j += 1
        assert j < len(html_text), "layout object not found after data array"
    layout, _ = decoder.raw_decode(html_text, j)
    return _decode_plotly_arrays(data), _decode_plotly_arrays(layout)


def _decode_plotly_arrays(obj):
    """Recursively convert plotly base64 array dicts to plain lists.

    Plotly serializes numpy arrays as ``{"dtype": "<code>", "bdata": "<b64>"}``.
    This walks dicts/lists and replaces those entries with the decoded numeric
    list. ``dtype`` codes: f8=float64, i1=int8, i2=int16, i4=int32, i8=int64,
    u1=uint8, u2=uint16, u4=uint32, u8=uint64, f4=float32.
    """
    import base64
    import struct

    _dtypes = {
        "f8": ("d", 8), "f4": ("f", 4),
        "i1": ("b", 1), "i2": ("h", 2), "i4": ("i", 4), "i8": ("q", 8),
        "u1": ("B", 1), "u2": ("H", 2), "u4": ("I", 4), "u8": ("Q", 8),
    }
    if isinstance(obj, dict):
        if "dtype" in obj and "bdata" in obj and obj["dtype"] in _dtypes:
            fmt, size = _dtypes[obj["dtype"]]
            raw = base64.b64decode(obj["bdata"])
            n = len(raw) // size
            return list(struct.unpack(f"{n}{fmt}", raw))
        return {k: _decode_plotly_arrays(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_decode_plotly_arrays(v) for v in obj]
    return obj


# --------------------------------------------------------------------------- #
# pytest fixtures
# --------------------------------------------------------------------------- #
@pytest.fixture()
def tmp_out(tmp_path: Path) -> Path:
    """A fresh output prefix path inside a per-test tmp directory."""
    return tmp_path / "out"


@pytest.fixture()
def two_ref_bam(tmp_path: Path) -> Path:
    """The canonical 2-reference BAM (refA=50bp, refB=30bp) with reads on both."""
    bam = tmp_path / "multi.bam"
    make_bam(
        bam,
        refs=[("refA", 50), ("refB", 30)],
        reads=(
            [
                {"ref": "refA", "start": i * 2, "length": 20, "mapq": 60}
                for i in range(8)
            ]
            + [
                {"ref": "refB", "start": i * 3, "length": 30, "mapq": 60}
                for i in range(4)
            ]
        ),
    )
    return bam


@pytest.fixture()
def two_ref_vcf(tmp_path: Path) -> Path:
    """A VCF with records on both refA and refB."""
    vcf = tmp_path / "multi.vcf"
    make_vcf(
        vcf,
        records=[
            {"chrom": "refA", "pos": 5, "ref": "A", "alt": "T", "info": {"AF": 0.5}},
            {"chrom": "refA", "pos": 10, "ref": "C", "alt": "G", "info": {"AF": 0.3}},
            {"chrom": "refB", "pos": 7, "ref": "G", "alt": "A", "info": {"AF": 0.8}},
        ],
    )
    return vcf


@pytest.fixture()
def two_ref_bed(tmp_path: Path) -> Path:
    """A BED with records on both refA and refB."""
    bed = tmp_path / "multi.bed"
    make_bed(
        bed,
        rows=[
            {"chrom": "refA", "start": 1, "stop": 15, "name": "regionA1", "score": "0", "strand": "+"},
            {"chrom": "refB", "start": 5, "stop": 25, "name": "regionB1", "score": "0", "strand": "-"},
        ],
    )
    return bed


def _simple_gb(ref_id: str, length: int, feature_type: str = "CDS") -> SeqRecord:
    return make_gb_record(
        ref_id, length, features=[{"type": feature_type, "location": [(1, 10, 1)]}]
    )


@pytest.fixture()
def make_gb_path(tmp_path: Path):
    """Factory fixture: write GenBank records to a tmp file and return its path."""
    def _factory(records: list[SeqRecord], name: str = "rec.gb") -> Path:
        path = tmp_path / name
        make_gb(path, records)
        return path
    return _factory
