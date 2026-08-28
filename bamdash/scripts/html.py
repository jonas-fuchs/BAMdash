"""
contains defs for assembling the multi-reference master HTML
"""
# BUILT INS
import html

# LIBS
import json
import logging
import re
from importlib.resources import files
from string import Template

import plotly
import plotly.offline as pyo

logger = logging.getLogger("bamdash")


def _sanitize_ref_id(ref):
    """
    Produce a filesystem- and DOM-safe identifier from a reference id.

    :param ref: reference id string (may contain spaces, slashes, etc.)
    :return: sanitized identifier suitable for filenames and HTML element ids
    """
    # collapse any run of non-alphanumeric characters into a single underscore
    return re.sub(r"[^A-Za-z0-9]+", "_", ref).strip("_")


def _panel_for_ref(fig, ref):
    """
    Render a single reference figure as a wrapped panel div.

    :param fig: plotly figure for this reference
    :param ref: reference id (used to build the panel id)
    :return: ``(panel_html, safe_id)`` where panel_html is a ``<div>`` containing
        the plotly ``Plotly.newPlot`` snippet (no plotly.js bundle), and safe_id
        is the sanitized id used in the ``id`` attribute
    """
    safe_id = _sanitize_ref_id(ref)
    # responsive: true lets plotly fill the panel container; combined with the
    # absolutely-positioned panels in master.html every reference resolves its
    # 100% size against the same viewport box.
    plot_html = fig.to_html(
        full_html=False,
        include_plotlyjs=False,
        div_id=f"plotly-{safe_id}",
        config={"responsive": True},
    )
    panel = (f'<div id="plot-{safe_id}" class="ref-panel">\n{plot_html}\n</div>')
    return panel, safe_id


def build_master_html(ref_figures, prefix, offline=True):
    """
    Assemble a master HTML document containing one panel per reference and a
    dropdown to switch between them.

    :param ref_figures: ordered mapping ``{ref: {"fig": fig, ...}}``; order is
        preserved for the dropdown
    :param prefix: output path/prefix; the master HTML is written to
        ``f"{prefix}.html"``
    :param offline: when True (default) the plotly.js bundle is inlined into
        the document so it is fully usable without an internet connection.
        When False, plotly.js is loaded from a CDN instead, which produces a
        much smaller html file but requires internet access to view.
    """
    if offline:
        # inline the plotly.js bundle once (offline, ~3MB)
        plotly_js = pyo.get_plotlyjs()
        plotly_tag = f"<script>{plotly_js}</script>"
    else:
        # load plotly.js from the CDN; keeps the html small but needs internet.
        # NOTE: use the bundled plotly.js version, NOT plotly.__version__ (the
        # python package version). The two diverged long ago: plotly 6.5.0
        # ships plotly.js v3.3.0. Using the python version produces a CDN URL
        # for a JS release that does not exist (403), so plotly never loads and
        # no figure renders. get_plotlyjs_version() returns the exact version
        # of the JS bundle that fig.to_html() targets, so the CDN URL always
        # matches the figure's expected plotly.js.
        plotly_version = pyo.get_plotlyjs_version()
        plotly_tag = (
            f'<script src="https://cdn.plot.ly/plotly-{plotly_version}.min.js">'
            f"</script>"
        )

    dropdown_options = []
    panels = []
    # per-reference linear/log y-axis ranges so the global log/linear toggle
    # can apply the correct range to every subfigure (each reference has its
    # own coverage maximum and thus its own linear/log upper bound).
    axis_ranges = {}
    first = True
    for ref in ref_figures:
        fig = ref_figures[ref]["fig"]
        upper = ref_figures[ref]["upper"]
        log_upper = ref_figures[ref]["log_upper"]
        panel, safe_id = _panel_for_ref(fig, ref)
        # mark the first panel as active so it is visible on load
        if first:
            panel = panel.replace('class="ref-panel"', 'class="ref-panel active"', 1)
            first = False
        panels.append(panel)
        axis_ranges[safe_id] = {"linear": [0, upper], "log": [0, log_upper]}
        # escape the displayed label but keep the value as the safe id
        dropdown_options.append(f'    <option value="{safe_id}">{html.escape(ref)}</option>')

    template_path = files("bamdash.scripts.templates").joinpath("master.html")
    template = Template(template_path.read_text(encoding="utf-8"))
    master = template.safe_substitute(
        title=html.escape(prefix),
        plotly_js=plotly_tag,
        dropdown_options="\n".join(dropdown_options),
        panels="\n".join(panels),
        axis_ranges=json.dumps(axis_ranges),
    )

    with open(f"{prefix}.html", "w", encoding="utf-8") as fp:
        fp.write(master)
    logger.info("wrote master html with %d reference(s) to %s.html", len(ref_figures), prefix)
