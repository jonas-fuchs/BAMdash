"""
contains main workflow
"""

# BUILT-INS
import sys
import argparse
import math

# LIBS
from plotly.subplots import make_subplots

# BAMDASH
from bamdash.scripts import data
from bamdash.scripts import plotting
from bamdash.scripts import config
from bamdash import __version__
from . import _program


def get_args(sysargs):
    """
    arg parsing for bamdash
    """
    parser = argparse.ArgumentParser(
        prog=_program,
        usage='''\tbamdash -b "bam file path" -r "reference_id" [additional arguments]''')
    parser.add_argument(
        "-b",
        "--bam",
        required=True,
        type=str,
        metavar=" ",
        help="bam file location"
    )
    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        type=str,
        metavar=" ",
        help="reference id"
    )
    parser.add_argument(
        "-t",
        "--tracks",
        default=None,
        action="store",
        type=str,
        metavar="track_1",
        nargs="*",
        help="file location of tracks"
    )
    parser.add_argument(
        "--slider",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="show slider"
    )
    parser.add_argument(
        "-e",
        "--export_static",
        type=str,
        metavar="None",
        default=None,
        help="export as png, jpg, pdf, svg"
    )
    parser.add_argument(
        "-d",
        "--dimensions",
        type=int,
        metavar="px",
        action="store",
        default=[1920, 1080],
        nargs=2,
        help="width and height of the static image in px"
    )
    parser.add_argument(
        "-v",
        "--version",
        action='version',
        version=f"virheat {__version__}"
    )
    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        return parser.parse_args(sysargs)


def main(sysargs=sys.argv[1:]):
    """
    main function for data extraction and plotting
    """
    # parse args
    args = get_args(sysargs)
    # define subplot number and track heights
    track_heights = [1]
    if args.tracks is not None:
        number_of_tracks = len(args.tracks)+1
        for track in args.tracks:
            if track.endswith("vcf"):
                track_heights = track_heights + [config.vcf_track_proportion]
            elif track.endswith("gb"):
                track_heights = track_heights + [config.gb_track_proportion]
            elif track.endswith("bed"):
                track_heights = track_heights + [config.bed_track_proportion]
            else:
                sys.exit("one of the track types is not supported")
    else:
        number_of_tracks = 1
    # define layout
    fig = make_subplots(
        rows=number_of_tracks,
        cols=1,
        shared_xaxes=True,
        row_heights=track_heights,
        vertical_spacing=config.plot_spacing,
    )
    # create coverage plot
    coverage_df, title = data.bam_to_coverage_df(args.bam, args.reference)
    plotting.create_coverage_plot(fig, 1, coverage_df)
    # create track plots
    if args.tracks is not None:
        for index, track in enumerate(args.tracks):
            row = index+2
            if track.endswith("vcf"):
                vcf_df = data.vcf_to_df(track, args.reference)
                plotting.create_vcf_plot(fig, row, vcf_df)
            elif track.endswith("gb"):
                gb_dict = data.genbank_to_dict(track, coverage_df, args.reference)
                plotting.create_track_plot(fig, row, gb_dict, config.box_gb_size, config.box_gb_alpha)
            elif track.endswith("bed"):
                bed_dict = data.bed_to_dict(track, coverage_df, args.reference)
                plotting.create_track_plot(fig, row, bed_dict, config.box_bed_size, config.box_bed_alpha)

    # global formatting
    fig.update_layout(
        plot_bgcolor="white",
        hovermode="x unified",
        title=dict(
            text=title,
            x=1,
            font=dict(
                family="Arial",
                size=16,
                color='#000000'
            )
        ),
        font=dict(
            family="Arial",
            size=16,
        )
    )
    # global x axes
    fig.update_xaxes(
        mirror=False,
        ticks="outside",
        showline=True,
        linecolor="black",
        range=[0, max(coverage_df["position"])]
    )
    # global y axis
    fig.update_yaxes(
        mirror=False,
        ticks="outside",
        showline=True,
        linecolor="black"
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

    fig.write_html(f"{args.reference}_plot.html")
    if args.export_static is not None:
        # static image specific options
        if config.show_log:
            fig["layout"]["yaxis"]["type"] = "log"
            fig["layout"]["yaxis"]["range"] = (0, math.log(fig["layout"]["yaxis"]["range"][1], 10))
        fig.update_layout(updatemenus=[dict(visible=False)])
        # write static image
        fig.write_image(f"{args.reference}_plot.{args.export_static}", width=args.dimensions[0], height=args.dimensions[1])
