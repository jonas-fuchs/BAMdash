# python
"""
contains main workflow
"""

# BUILT-INS
import sys
import argparse
import logging
import math
import json

# LIBS
import plotly.io as pio
import pandas as pd
import pysam
from plotly.subplots import make_subplots
import kaleido

# BAMDASH
from bamdash.scripts import data
from bamdash.scripts import plotting
from bamdash.scripts import config
from bamdash import __version__

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
        metavar="REF_ID",
        help="seq reference id (default: first reference in bam)"
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
             "(allowed: html, png, jpg, jpeg, webp, svg, pdf, eps)"
    )
    parser.add_argument(
        "-q",
        "--quality-threshold",
        type=int,
        default=15,
        metavar="15",
        help="quality threshold for reads"
    )
    parser.add_argument(
        "-bs",
        "--binsize",
        default=1,
        type=int,
        metavar="N",
        help="bins for the coverage plot"
    )
    parser.add_argument(
        "-t",
        "--tracks",
        default=None,
        type=str,
        metavar="track_1",
        nargs="*",
        help="file location of tracks (accepted: *.vcf, *.bed, *.gb)"
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
        help="dump annotated track data; filenames derive from --prefix"
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


def main(sysargs=sys.argv[1:]):
    """
    main function for data extraction and plotting
    """
    # parse args
    args = get_args(sysargs)

    # resolve ref id: fall back to the first reference in the bam header
    if args.ref_id is None:
        with pysam.AlignmentFile(args.bam, "rb") as bam:
            args.ref_id = bam.references[0]
        logger.info("no ref-id given, using '%s'", args.ref_id)

    # define subplot number, track heights and parse data
    coverage_df, title, stat_dict = data.bam_to_coverage_df(args.bam, args.ref_id, args.coverage, args.quality_threshold)

    track_heights = [1]
    track_data = []
    # extract data and check if ref was found
    if args.tracks is not None:
        number_of_tracks = len(args.tracks)+1
        for track in args.tracks:
            if track.endswith("vcf"):
                vcf_data = [data.vcf_to_df(track, args.ref_id), "vcf"]
                if vcf_data[0].empty:
                    logger.warning("vcf data does not contain the seq reference id")
                    number_of_tracks -= 1
                else:
                    track_heights = track_heights + [config.vcf_track_proportion]
                    track_data.append(vcf_data)
            elif track.endswith("gb"):
                gb_dict, seq = data.genbank_to_dict(track, coverage_df, args.ref_id, args.coverage)
                if gb_dict:
                    track_heights = track_heights + [config.gb_track_proportion]
                    track_data.append([gb_dict, "gb", seq])
                else:
                    logger.warning("gb data does not contain the seq reference id")
                    number_of_tracks -= 1
            elif track.endswith("bed"):
                bed_data = [data.bed_to_dict(track, coverage_df, args.ref_id, args.coverage), "bed"]
                if bed_data[0]["bed annotations"]:
                    track_heights = track_heights + [config.bed_track_proportion]
                    track_data.append(bed_data)
                else:
                    logger.warning("bed data does not contain the seq reference id")
                    number_of_tracks -= 1
            else:
                logger.error("one of the track types is not supported (supported are *.vcf, *.bed and *.gb")
                sys.exit(1)
    else:
        number_of_tracks = 1

    # annotate if one gb and vcfs are in tracks
    track_data = data.annotate_vcfs_in_tracks(track_data)

    # define layout
    fig = make_subplots(
        rows=number_of_tracks,
        cols=1,
        shared_xaxes=True,
        row_heights=track_heights,
        vertical_spacing=config.plot_spacing,
    )
    # create coverage plot
    plotting.create_coverage_plot(fig, 1, coverage_df, args.binsize)
    # create track plots
    if track_data:
        for index, track in enumerate(track_data):
            row = index+2
            if track[1] == "vcf":
                plotting.create_vcf_plot(fig, row, track[0])
            elif track[1] == "gb":
                plotting.create_track_plot(fig, row, track[0], config.box_gb_size, config.box_gb_alpha)
            elif track[1] == "bed":
                plotting.create_track_plot(fig, row, track[0], config.box_bed_size, config.box_bed_alpha)

    # define own templates
    pio.templates["plotly_dark_custom"], pio.templates["plotly_white_custom"] = pio.templates["plotly_dark"], pio.templates["plotly_white"]
    # change params
    pio.templates["plotly_dark_custom"].update(
        layout=dict(yaxis=dict(linecolor="white", tickcolor="white", zerolinecolor="rgb(17,17,17)"),
                    xaxis=dict(linecolor="white", tickcolor="white", zerolinecolor="rgb(17,17,17)"),
                    updatemenudefaults=dict(bgcolor="rgb(115, 115, 115)")
                    )
    )
    pio.templates["plotly_white_custom"].update(
        layout=dict(yaxis=dict(linecolor="black", tickcolor="black", zerolinecolor="white"),
                    xaxis=dict(linecolor="black", tickcolor="black", zerolinecolor="white"),
                    updatemenudefaults=dict(bgcolor="rgb(204, 204, 204)")
                    )
    )

    # compute upper thresholds for y axis
    upper = max(coverage_df["coverage"]+max(coverage_df["coverage"]/10))
    log_upper = max(1, math.ceil(math.log10(upper)))

    # global formatting
    fig.update_layout(
        template="plotly_white_custom",
        hovermode="closest",
        font=dict(
            family=config.font,
            size=config.font_size,
        ),
        # Add buttons
        updatemenus=[
            dict(
                type="buttons",
                direction="left",
                buttons=[
                    dict(
                        args=[{"yaxis.type": "linear", "yaxis.range": [0, upper], "yaxis.autorange": False}],
                        label="linear",
                        method="relayout"
                    ),
                    dict(
                        args=[{"yaxis.type": "log", "yaxis.range": [0, log_upper], "yaxis.autorange": False}],
                        label="log",
                        method="relayout"
                    ),
                    dict(args=[{"template": pio.templates["plotly_dark_custom"], "visible": True}],
                         label="dark",
                         method="relayout"),
                    dict(args=[{"template": pio.templates["plotly_white_custom"], "visible": True}],
                         label="light",
                         method="relayout"),
                ],
                pad={"r": 10, "t": 1},
                showactive=False,
                xanchor="left",
                y=1.15,
                yanchor="top"
            )
        ],
        # add global stats as annotation
        annotations=[
            dict(text=title, y=1.14, yref="paper",
                 align="center", showarrow=False)
        ]
    )
    # global x axes
    fig.update_xaxes(
        mirror=False,
        showline=True,
        linewidth=1,
        ticks="outside",
        minor_ticks="outside",
        range=[0-max(coverage_df["position"])/50, max(coverage_df["position"])+max(coverage_df["position"])/50],
        showgrid=False,
    )
    # global y axis
    fig.update_yaxes(
        mirror=False,
        zeroline=False,
        showline=True,
        linewidth=1,
        ticks="outside",
        minor_ticks="outside",
        showgrid=False
    )
    # if a range slider is shown, do not display the xaxis title
    # (will be shown underneath)
    if args.slider:
        # last y axis
        fig.update_xaxes(
            rangeslider=dict(
                visible=True,
                thickness=0.05
            ),
            row=number_of_tracks
        )
    else:
        # last x axis
        fig.update_xaxes(title_text="genome position", row=number_of_tracks, col=1)
    # export one file per requested suffix
    static_prepared = False
    for suffix in args.suffix:
        if suffix == "html":
            fig.write_html(f"{args.prefix}.html")
        else:
            # apply static-image layout cleanup once before the first static export
            if not static_prepared:
                if config.show_log:  # correct log layout
                    fig.update_yaxes(type="log", dtick=1, row=1, range=[0, log_upper], autorange=False)
                fig.update_layout(updatemenus=[dict(visible=False)])  # no buttons
                fig.update_layout(annotations=[dict(visible=False)])  # no annotations
                static_prepared = True
            fig.write_image(f"{args.prefix}.{suffix}", width=args.dimensions[0], height=args.dimensions[1])

    # dump track data
    vcf_track_count, bed_track_count, gb_track_count = 0, 0, 0
    if args.dump:
        pd.DataFrame.from_dict(stat_dict, orient="index").to_csv(f"{args.prefix}_bam_stats.tabular", sep="\t", header=False, index=True)
        if track_data:
            for track in track_data:
                if track[1] == "vcf":
                    track[0] = track[0].drop(['position_jittered'], axis=1)  # do not report jittered position
                    track[0].to_csv(f"{args.prefix}_vcf_data_{vcf_track_count}.tabular", sep="\t", header=True, index=False)
                    vcf_track_count += 1
                elif track[1] == "bed":
                    bed_df = pd.DataFrame.from_dict(track[0]["bed annotations"], orient="index")
                    bed_df.drop("track", axis=1, inplace=True)
                    bed_df.to_csv(f"{args.prefix}_bed_data_{bed_track_count}.tabular", sep="\t", header=True, index=False)
                    bed_track_count += 1
                elif track[1] == "gb":
                    with open(f"{args.prefix}_gb_data_{gb_track_count}.json", "w") as fp:
                        json.dump(track[0], fp)
                    gb_track_count += 1