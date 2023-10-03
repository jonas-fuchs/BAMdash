# BUILT-INS
import sys
import argparse

# bamcov
from bamcov import __version__
from . import _program


def get_args(sysargs):
    """
    arg parsing for virheat
    """
    parser = argparse.ArgumentParser(
        prog=_program,
        usage='''\tbamcov -b "bam file path" -r "reference_id" [additional arguments]''')
    parser.add_argument(
        "-b",
        "--bam",
        type=str,
        metavar=" ",
        help="bam file location"
    )
    parser.add_argument(
        "-r",
        "--reference",
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
        help="export as png, pdf, svg"
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
    print("todo")

