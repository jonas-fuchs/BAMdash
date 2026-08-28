"""
BAMdash plotting configuration.

The defaults are stored in the shipped ``config.toml`` (next to this module)
and loaded into module-level attributes at import time, so the rest of the
codebase can keep using ``config.coverage_fill_color`` etc.

A user can override individual settings by supplying their own TOML file via
the ``--custom-config`` CLI argument. ``load_config`` overlays only the keys
present in the user file on top of the shipped defaults; everything else is
left untouched.

Use ``load_config(custom_path)`` from ``command.main`` *before* any plotting
code runs so the overlay takes effect. Calling it with ``custom_path=None``
(or no argument) re-applies the shipped defaults, which is useful for tests
that need a clean slate.
"""
import logging
import sys
import tomllib
from importlib.resources import files
from typing import Any

logger = logging.getLogger("bamdash")


# --------------------------------------------------------------------------- #
# defaults
# --------------------------------------------------------------------------- #
# These mirror config.toml and act as the ultimate fallback if the shipped
# TOML is ever missing or incomplete. They are also the source of the canonical
# key set: only keys listed here are recognised, anything else in a user TOML
# is rejected with an error so typos do not silently get ignored.
_DEFAULTS: dict[str, Any] = {
    # pdf settings
    "show_log": True,
    # overall layout
    "vcf_track_proportion": 0.3,
    "gb_track_proportion": 0.5,
    "bed_track_proportion": 0.2,
    "plot_spacing": 0.05,
    "font": "Arial",
    "font_size": 16,
    # coverage customize
    "coverage_fill_color": "rgba(255, 212, 135, 0.4)",
    "coverage_line_color": "rgba(224, 168, 68, 1)",
    "average_line_color": "grey",
    "average_line_width": 1,
    # track customize
    "track_color_scheme": "agsunset",
    "track_color_single": "rgb(145, 145, 145)",
    "strand_types": ["triangle-right", "triangle-left", "diamond-wide"],
    "strand_marker_size": 8,
    "strand_marker_line_width": 1,
    "strand_marker_line_color": "rgba(0, 0, 0, 0.2)",
    "box_bed_alpha": [0.7, 0.7],
    "box_bed_size": [0.4, 0.4],
    "box_gb_alpha": [0.7, 0.8],
    "box_gb_size": [0.4, 0.3],
    # variant customize
    "variant_marker_size": 10,
    "variant_marker_line_width": 1,
    "variant_line_color": "black",
    "stem_color": "grey",
    "stem_width": 1,
    "snp_color": "grey",
    "ins_color": "blue",
    "del_color": "red",
}


def _read_toml(path) -> dict[str, Any]:
    """Read a TOML file (``str`` / ``Path`` / resource path) into a flat dict."""
    with open(path, "rb") as fp:
        return tomllib.load(fp)


def _shipped_toml_path():
    """Return the path to the shipped ``config.toml`` resource."""
    return files("bamdash.scripts").joinpath("config.toml")


def _apply(values: dict[str, Any]) -> None:
    """Set every key in ``values`` as a module-level attribute.

    ``values`` must already be validated (only known keys, correct types).
    """
    module = sys.modules[__name__]
    for key, value in values.items():
        setattr(module, key, value)


def _validate(overrides: dict[str, Any], source: str) -> dict[str, Any]:
    """Reject unknown keys so typos in a user TOML do not get silently ignored.

    :param overrides: parsed TOML dict to check
    :param source: human-readable description of the file (for the error msg)
    :return: ``overrides`` unchanged if valid
    :raises KeyError: if an unknown key is present
    """
    unknown = [k for k in overrides if k not in _DEFAULTS]
    if unknown:
        raise KeyError(
            f"unknown config key(s) in {source}: {', '.join(sorted(unknown))}. "
            f"Allowed keys: {', '.join(sorted(_DEFAULTS))}"
        )
    return overrides


def load_config(custom_path: str | None = None) -> None:
    """(Re)load the configuration into module-level attributes.

    The shipped ``config.toml`` is always applied first, establishing the
    defaults. When ``custom_path`` is given, only the keys present in that
    file are overlaid on top — all other settings keep their default values.

    Call this from ``command.main`` *before* any plotting code runs. Calling
    with ``custom_path=None`` resets to the shipped defaults (useful in tests).

    :param custom_path: optional path to a user-supplied TOML file whose keys
        override the shipped defaults
    :raises FileNotFoundError: if ``custom_path`` does not exist
    :raises KeyError: if the user TOML contains an unknown config key
    """
    # start from the hard-coded fallback defaults
    _apply(_DEFAULTS)

    # overlay the shipped config.toml (overrides the fallbacks; should match).
    # config.toml is packaged as required package-data, so a missing file is a
    # real packaging bug — let the FileNotFoundError surface rather than
    # silently falling back to _DEFAULTS.
    shipped = _read_toml(_shipped_toml_path())
    _apply(_validate(shipped, "shipped config.toml"))

    # overlay the user-supplied config, if any
    if custom_path is not None:
        overrides = _read_toml(custom_path)
        _apply(_validate(overrides, f"custom config {custom_path}"))
        logger.info("applied custom config from %s", custom_path)


# load the shipped defaults at import time so `from bamdash.scripts import config`
# yields a fully-populated module without an explicit load_config() call.
load_config()
