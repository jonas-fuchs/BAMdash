"""
contains main workflow
"""

# BUILT-INS
import sys
import argparse
import math
import json

# LIBS
import plotly.io as pio
import pandas as pd
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
        help="seq reference id"
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
        "--dump",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="dump annotated track data"
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

    # define subplot number, track heights and parse data
    coverage_df, title = data.bam_to_coverage_df(args.bam, args.reference, args.coverage)
    track_heights = [1]
    track_data = []
    # extract data and check if ref was found
    if args.tracks is not None:
        number_of_tracks = len(args.tracks)+1
        for track in args.tracks:
            if track.endswith("vcf"):
                vcf_data = [data.vcf_to_df(track, args.reference), "vcf"]
                if vcf_data[0].empty:
                    print("WARNING: vcf data does not contain the seq reference id")
                    number_of_tracks -= 1
                else:
                    track_heights = track_heights + [config.vcf_track_proportion]
                    track_data.append(vcf_data)
            elif track.endswith("gb"):
                gb_dict, seq = data.genbank_to_dict(track, coverage_df, args.reference, args.coverage)
                if gb_dict:
                    track_heights = track_heights + [config.gb_track_proportion]
                    track_data.append([gb_dict, "gb", seq])
                else:
                    print("WARNING: gb data does not contain the seq reference id")
                    number_of_tracks -= 1
            elif track.endswith("bed"):
                bed_data = [data.bed_to_dict(track, coverage_df, args.reference, args.coverage), "bed"]
                if bed_data[0]["bed annotations"]:
                    track_heights = track_heights + [config.bed_track_proportion]
                    track_data.append(bed_data)
                else:
                    print("WARNING: bed data does not contain the seq reference id")
                    number_of_tracks -= 1
            else:
                sys.exit("one of the track types is not supported (supported are *.vcf, *.bed and *.gb")
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
    plotting.create_coverage_plot(fig, 1, coverage_df)
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

    # global formatting
    fig.update_layout(
        template="plotly_white_custom",
        hovermode="x unified",
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
                        args=["yaxis.type", "linear"],
                        label="linear",
                        method="relayout"
                    ),
                    dict(
                        args=["yaxis.type", "log"],
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
        range=[0, max(coverage_df["position"])],
        showgrid=False,
    )
    # global y axis
    fig.update_yaxes(
        mirror=False,
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
    # html export
    fig.write_html(f"{args.reference}_plot.html")
    # static image export
    if args.export_static is not None:
        # static image specific options
        if config.show_log:  # correct log layout
            fig["layout"]["yaxis"]["type"] = "log"
            fig["layout"]["yaxis"]["range"] = (0, math.log(fig["layout"]["yaxis"]["range"][1], 10))
            fig.update_yaxes(dtick=1, row=1)
        fig.update_layout(updatemenus=[dict(visible=False)])  # no buttons
        fig.update_layout(annotations=[dict(visible=False)])  # no annotations
        # write static image
        pio.kaleido.scope.mathjax = None  # fix so no weird box is shown
        fig.write_image(f"{args.reference}_plot.{args.export_static}", width=args.dimensions[0], height=args.dimensions[1])

    # dump track data
    vcf_track_count, bed_track_count, gb_track_count = 0, 0, 0
    if args.dump and track_data:
        for track in track_data:
            if track[1] == "vcf":
                track[0].to_csv(f"{args.reference}_vcf_data_{vcf_track_count}.tabular", sep="\t", header=True, index=False)
                vcf_track_count += 1
            elif track[1] == "bed":
                bed_df = pd.DataFrame.from_dict(track[0]["bed annotations"], orient="index")
                bed_df.drop("track", axis=1, inplace=True)
                bed_df.to_csv(f"{args.reference}_bed_data_{bed_track_count}.tabular", sep="\t", header=True, index=False)
                bed_track_count += 1
            elif track[1] == "gb":
                with open(f"{args.reference}_gb_data_{gb_track_count}.json", "w") as fp:
                    json.dump(track[0], fp)
                gb_track_count += 1
