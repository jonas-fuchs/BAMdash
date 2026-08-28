# python
"""
contains main workflow
"""

# BUILT-INS
import argparse
import copy
import json
import logging
import os
import sys

# LIBS
import pandas as pd
import pysam

from bamdash import __version__

# BAMDASH
from bamdash.scripts import config, data, plotting
from bamdash.scripts import html as html_mod

logger = logging.getLogger("bamdash")


def get_args(sysargs):
    """
    arg parsing for bamdash
    """
    parser = argparse.ArgumentParser(
        usage='''\tbamdash -b "bam file path" [additional arguments]''')
    parser.add_argument(
        "-b",
        "--bam",
        required=True,
        type=str,
        metavar="BAM",
        help="bam file location"
    )
    parser.add_argument(
        "-r",
        "--ref-id",
        required=False,
        default=None,
        type=str,
        nargs="*",
        metavar="REF_ID",
        help="seq reference id(s) for which to generate plots; default: all references in bam if argument is omitted"
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        default="./plot",
        metavar="./plot",
        help="path and partial filename for output files"
    )
    parser.add_argument(
        "-s",
        "--suffix",
        type=str,
        nargs="*",
        default=["html"],
        metavar="html",
        help="output file extensions appended to prefix "
             "(allowed: html, png, jpg, jpeg, webp, svg, pdf, eps), multiple values allowed, default: html"
    )
    parser.add_argument(
        "-q",
        "--quality-threshold",
        type=int,
        default=15,
        metavar="15",
        help="quality threshold for reads (default: 15)"
    )
    parser.add_argument(
        "-bs",
        "--binsize",
        default=1,
        type=int,
        metavar="N",
        help="bins the coverage data into N bp bins (default: 1, no binning)"
    )
    parser.add_argument(
        "-t",
        "--tracks",
        default=None,
        type=str,
        metavar="track_1",
        nargs="*",
        help="file location of tracks (accepted: *.vcf, *.bed, *.gb, multiple paths to files allowed)"
    )
    parser.add_argument(
        "-c",
        "--coverage",
        default=5,
        type=int,
        metavar="5",
        help="minimum coverage"
    )
    parser.add_argument(
        "--slider",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="show slider"
    )
    parser.add_argument(
        "--offline",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="inline the plotly.js bundle into the output html so the file is "
             "fully usable without an internet connection (default: True). "
             "Use --no-offline to load plotly.js from a CDN instead, which "
             "produces a much smaller html file but requires internet access "
             "to view"
    )
    parser.add_argument(
        "--custom-config",
        default=None,
        type=str,
        metavar="TOML",
        help="path to a user-supplied TOML config file whose settings override "
             "the shipped defaults (see config.toml). Only the keys present in "
             "the file are overwritten; all others keep their default values"
    )
    parser.add_argument(
        "-d",
        "--dimensions",
        type=int,
        metavar="px",
        default=None,
        nargs=2,
        help="width and height of static (non-html) output in px "
             "(default: 1920 1080; ignored when only html is requested)"
    )
    parser.add_argument(
        "--dump",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="dump annotated track data; filenames derive from --prefix (default: no dump)"
    )
    parser.add_argument(
        "--verbose",
        action="count",
        default=0,
        help="increase logging verbosity: --verbose for INFO, --verbose --verbose for DEBUG"
    )
    parser.add_argument(
        "-v",
        "--version",
        action='version',
        version=f"bamdash {__version__}"
    )
    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)
    # configure logging early so warnings emitted below are formatted
    # default WARNING, --verbose -> INFO, --verbose --verbose -> DEBUG
    level = logging.WARNING - 10 * min(args.verbose, 2)
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")
    # a bare "-s" with no values is almost certainly a mistake
    if not args.suffix:
        parser.error("the following arguments are required: at least one --suffix / -s value")
    # drop duplicate suffixes while preserving order
    args.suffix = list(dict.fromkeys(args.suffix))
    # validate suffixes against the supported set
    allowed_suffixes = {"html", "png", "jpg", "jpeg", "webp", "svg", "pdf", "eps"}
    bad = [s for s in args.suffix if s not in allowed_suffixes]
    if bad:
        parser.error(f"unsupported --suffix value(s): {', '.join(bad)} "
                     f"(allowed: {', '.join(sorted(allowed_suffixes))})")
    # --dimensions only affects static (non-html) output
    static_suffixes = [s for s in args.suffix if s != "html"]
    if args.dimensions is not None and not static_suffixes:
        logger.warning("--dimensions has no effect when no static (non-html) suffix is requested")
    if args.dimensions is None:
        args.dimensions = [1920, 1080]
    return args


def _ensure_bam_index(bam_path):
    """
    Create the bam index (``.bai``) if it is missing.

    pysam requires an index for the ``count_coverage`` / ``get_index_statistics``
    calls bamdash relies on. ``pysam.index`` builds it via samtools when no
    matching ``.bai`` exists next to the bam file.

    :param bam_path: path to the ``.bam`` file
    :raises FileNotFoundError: if the bam file itself does not exist
    :return: path to the index file that is now guaranteed to exist
    """
    if not os.path.isfile(bam_path):
        # let the downstream AlignmentFile open raise the canonical error; but
        # give a clear message first since pysam's message can be cryptic.
        raise FileNotFoundError(f"bam file not found: {bam_path}")

    # the canonical index path is the bam path with a trailing ".bai"
    index_path = bam_path + ".bai"
    if os.path.isfile(index_path):
        return index_path

    logger.info("no bam index found, creating %s", index_path)
    # pysam.index returns the index path on success; building the index in
    # place (no explicit output path) writes bam_path + ".bai".
    ret = pysam.index(bam_path)
    # older pysam versions return the path; newer ones return None on success
    return ret or index_path


def _load_reference(args, ref, refs, warned_no_match=None):
    """
    Load coverage and track data for a single reference and build its figure.

    :param args: parsed CLI args
    :param ref: reference id being processed
    :param refs: full list of requested reference ids; used to decide whether a
        non-matching track file warrants a warning (a file that matches *some*
        requested ref is silently skipped for the refs it does not contain)
    :param warned_no_match: set of track paths already warned about in
        ``main`` (files matching no requested ref); suppresses duplicate
        per-reference warnings
    :return: ``(ref_fig, upper, log_upper, track_data, stat_dict, number_of_tracks)``
        or ``None`` if the reference has no usable data (e.g. zero-length)
    """
    # coverage is required for every reference; a missing ref raises ReferenceNotFoundError
    coverage_df, title, stat_dict = data.bam_to_coverage_df(
        args.bam, ref, args.coverage, args.quality_threshold)

    track_heights = [1]
    track_data = []
    if args.tracks is not None:
        number_of_tracks = len(args.tracks) + 1
        for track in args.tracks:
            if track.endswith("vcf"):
                vcf_data = [data.vcf_to_df(track, ref), "vcf"]
                if vcf_data[0].empty:
                    logger.warning("vcf data does not contain the seq reference id '%s'", ref)
                    number_of_tracks -= 1
                else:
                    track_heights = track_heights + [config.vcf_track_proportion]
                    track_data.append(vcf_data)
            elif track.endswith("gb"):
                gb_dict, seq = data.genbank_to_dict(track, coverage_df, ref, args.coverage)
                if gb_dict:
                    track_heights = track_heights + [config.gb_track_proportion]
                    track_data.append([gb_dict, "gb", seq])
                else:
                    # Only warn if this gb file matches none of the requested
                    # references at all; otherwise it simply belongs to a
                    # different reference and is expected to be empty here.
                    # Files already warned about in main() are skipped to
                    # avoid one warning per reference.
                    if warned_no_match is None or track not in warned_no_match:
                        gb_ids = data.gb_record_ids(track)
                        if not (set(refs) & gb_ids):
                            logger.warning(
                                "gb data in %s does not contain any of the "
                                "requested reference id(s) %s", track, refs)
                    number_of_tracks -= 1
            elif track.endswith("bed"):
                bed_data = [data.bed_to_dict(track, coverage_df, ref, args.coverage), "bed"]
                # bed_to_dict keys the feature dict by the filename stem
                if bed_data[0] and next(iter(bed_data[0].values())):
                    track_heights = track_heights + [config.bed_track_proportion]
                    track_data.append(bed_data)
                else:
                    logger.warning("bed data does not contain the seq reference id '%s'", ref)
                    number_of_tracks -= 1
            else:
                logger.error("one of the track types is not supported (supported are *.vcf, *.bed and *.gb")
                sys.exit(1)
    else:
        number_of_tracks = 1

    # annotate if one gb and vcfs are in tracks
    track_data = data.annotate_vcfs_in_tracks(track_data)

    fig, upper, log_upper = plotting.build_figure(
        ref, coverage_df, track_data, number_of_tracks, track_heights, title, args)

    return fig, upper, log_upper, track_data, stat_dict, number_of_tracks


def _dump_reference(args, ref, track_data, stat_dict, multi):
    """
    Write the per-reference dump sidecars for ``--dump``.

    :param args: parsed CLI args
    :param ref: reference id
    :param track_data: track data list for this reference
    :param stat_dict: bam stats dict for this reference
    :param multi: whether more than one reference is being plotted; when False
        the original ``{prefix}_*`` names are used (no ref token) for backward
        compatibility
    """
    if multi:
        safe_ref = html_mod._sanitize_ref_id(ref)
        prefix = f"{args.prefix}_{safe_ref}"
    else:
        prefix = args.prefix
    vcf_track_count, bed_track_count, gb_track_count = 0, 0, 0
    pd.DataFrame.from_dict(stat_dict, orient="index").to_csv(
        f"{prefix}_bam_stats.tabular", sep="\t", header=False, index=True)
    if track_data:
        for track in track_data:
            if track[1] == "vcf":
                track[0] = track[0].drop(['position_jittered'], axis=1)  # do not report jittered position
                track[0].to_csv(f"{prefix}_vcf_data_{vcf_track_count}.tabular", sep="\t", header=True, index=False)
                vcf_track_count += 1
            elif track[1] == "bed":
                bed_key = next(iter(track[0]))
                bed_df = pd.DataFrame.from_dict(track[0][bed_key], orient="index")
                bed_df.drop("track", axis=1, inplace=True)
                bed_df.to_csv(f"{prefix}_bed_data_{bed_track_count}.tabular", sep="\t", header=True, index=False)
                bed_track_count += 1
            elif track[1] == "gb":
                with open(f"{prefix}_gb_data_{gb_track_count}.json", "w") as fp:
                    json.dump(track[0], fp)
                gb_track_count += 1


def main(sysargs=sys.argv[1:]):
    """
    main function for data extraction and plotting
    """
    # parse args
    args = get_args(sysargs)

    # apply the (possibly user-overridden) plotting config before any
    # plotting code runs. Unknown keys or a missing custom file are errors.
    if args.custom_config is not None:
        try:
            config.load_config(args.custom_config)
        except (FileNotFoundError, KeyError) as exc:
            logger.error("%s", exc)
            sys.exit(1)

    # ensure the bam index exists; pysam needs it for coverage/index stats
    try:
        _ensure_bam_index(args.bam)
    except FileNotFoundError as exc:
        logger.error("%s", exc)
        sys.exit(1)

    # resolve reference ids: default to all references in the bam header
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        available_refs = list(bam.references)
    if not args.ref_id:
        refs = available_refs
        logger.info("no ref-id given, using all %d reference(s) in bam", len(refs))
    else:
        refs = args.ref_id
        missing = [r for r in refs if r not in available_refs]
        if missing:
            logger.error(
                "ref id(s) %s do not exist in bam file. Available references are %s",
                missing, available_refs)
            sys.exit(1)

    # Pre-screen track files that cannot match any requested reference so we
    # can warn exactly once per offending file instead of once per reference.
    warned_no_match = set()
    if args.tracks:
        for track in args.tracks:
            if track.endswith("gb") and not (set(refs) & data.gb_record_ids(track)):
                logger.warning(
                    "gb data in %s does not contain any of the requested "
                    "reference id(s) %s", track, refs)
                warned_no_match.add(track)

    # load data and build a figure per reference
    ref_figures = {}
    for ref in refs:
        try:
            fig, upper, log_upper, track_data, stat_dict, number_of_tracks = _load_reference(
                args, ref, refs, warned_no_match)
        except data.ReferenceNotFoundError as exc:
            # already validated above, but guard against races / zero-length refs
            logger.error("%s", exc)
            sys.exit(1)
        ref_figures[ref] = {
            "fig": fig,
            "upper": upper,
            "log_upper": log_upper,
            "track_data": track_data,
            "stat_dict": stat_dict,
            "number_of_tracks": number_of_tracks,
        }
        logger.info("built figure for reference '%s'", ref)

    multi = len(ref_figures) >= 2

    # export one file per requested suffix
    for suffix in args.suffix:
        if suffix == "html":
            # Always build the master HTML, even for a single reference. This
            # keeps the output consistent: the dropdown just lists the one
            # available reference and the global log/linear toggle is always
            # present next to it.
            html_mod.build_master_html(ref_figures, args.prefix, args.offline)
        else:
            # static image: one file per reference as {prefix}_{ref}.{suffix}
            for ref, payload in ref_figures.items():
                fig = payload["fig"]
                log_upper = payload["log_upper"]
                # apply static-image layout cleanup on a copy so the HTML figure
                # (already written or about to be written) keeps its buttons.
                # deepcopy preserves the subplot grid (JSON round-trip does not).
                static_fig = copy.deepcopy(fig)
                plotting.prepare_static(static_fig, log_upper)
                if multi:
                    safe_ref = html_mod._sanitize_ref_id(ref)
                    out = f"{args.prefix}_{safe_ref}.{suffix}"
                else:
                    out = f"{args.prefix}.{suffix}"
                static_fig.write_image(out, width=args.dimensions[0], height=args.dimensions[1])

    # dump track data (per reference)
    if args.dump:
        for ref, payload in ref_figures.items():
            _dump_reference(args, ref, payload["track_data"], payload["stat_dict"], multi)