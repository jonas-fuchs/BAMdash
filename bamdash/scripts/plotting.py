"""
contains defs for plotting
"""

# BUILT-INS
import logging
import math
import statistics
import sys
from collections import Counter

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

from bamdash.scripts import config

logger = logging.getLogger("bamdash")


def _ensure_custom_templates():
    """Register the plotly_white_custom template (idempotent).

    """
    pio.templates["plotly_white_custom"] = pio.templates["plotly_white"]
    pio.templates["plotly_white_custom"].update(
        layout={"yaxis": {"linecolor": "black", "tickcolor": "black", "zerolinecolor": "white"},
                    "xaxis": {"linecolor": "black", "tickcolor": "black", "zerolinecolor": "white"},
                    "updatemenudefaults": {"bgcolor": "rgb(204, 204, 204)"}
                    }
    )


def build_figure(ref, coverage_df, track_data, number_of_tracks, track_heights,
                 title, args):
    """
    Build a complete plotly subplot figure for a single reference.

    Encapsulates the figure construction, track plotting, custom template,
    stats annotation and axis formatting that previously lived inline in
    ``command.main``. Called once per reference. The linear/log y-axis
    toggle is provided globally by the master HTML page (next to the
    reference dropdown) and applies to all subfigures at once, so no
    per-figure updatemenu buttons are added here.

    :param ref: reference id (used only for logging context)
    :param coverage_df: per-position coverage dataframe for this reference
    :param track_data: list of ``[data, type]`` / ``[data, "gb", seq]`` tracks for this ref
    :param number_of_tracks: total subplot rows (coverage + non-empty tracks) for this ref
    :param track_heights: relative row heights for ``make_subplots``
    :param title: stats annotation string (ref-specific) for this reference
    :param args: parsed CLI args (uses binsize, slider)
    :return: ``(fig, upper, log_upper)`` where upper/log_upper are the linear/log
        y-axis thresholds for the coverage axis of this reference
    """
    _ensure_custom_templates()

    fig = make_subplots(
        rows=number_of_tracks,
        cols=1,
        shared_xaxes=True,
        row_heights=track_heights,
        vertical_spacing=config.plot_spacing,
    )
    # coverage plot
    create_coverage_plot(fig, 1, coverage_df, args.binsize)
    # track plots
    if track_data:
        for index, track in enumerate(track_data):
            row = index + 2
            if track[1] == "vcf":
                create_vcf_plot(fig, row, track[0])
            elif track[1] == "gb":
                create_track_plot(fig, row, track[0], config.box_gb_size, config.box_gb_alpha)
            elif track[1] == "bed":
                create_track_plot(fig, row, track[0], config.box_bed_size, config.box_bed_alpha)

    # compute upper thresholds for y axis
    upper = max(coverage_df["coverage"] + max(coverage_df["coverage"] / 10))
    log_upper = max(1, math.ceil(math.log10(upper)))

    # global formatting. No per-figure linear/log updatemenu buttons are
    # added: the master HTML page provides a single global log/linear
    # toggle (next to the reference dropdown) that relayouts every
    # subfigure at once, so per-figure buttons would be redundant.
    fig.update_layout(
        template="plotly_white_custom",
        hovermode="closest",
        font={
            "family": config.font,
            "size": config.font_size,
        },
        # add global stats as annotation
        annotations=[
            {"text": title, "y": 1.14, "yref": "paper",
                 "align": "center", "showarrow": False}
        ]
    )
    # global x axes
    fig.update_xaxes(
        mirror=False,
        showline=True,
        linewidth=1,
        ticks="outside",
        minor_ticks="outside",
        range=[0 - max(coverage_df["position"]) / 50, max(coverage_df["position"]) + max(coverage_df["position"]) / 50],
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
    if args.slider:
        fig.update_xaxes(
            rangeslider={
                "visible": True,
                "thickness": 0.05
            },
            row=number_of_tracks
        )
    else:
        fig.update_xaxes(title_text="genome position", row=number_of_tracks, col=1)

    return fig, upper, log_upper


def prepare_static(fig, log_upper):
    """
    Apply the static-image layout cleanup to a figure in place: log y-axis if
    configured, hide updatemenus and annotations. Used before ``write_image``.

    :param fig: plotly figure to mutate
    :param log_upper: log y-axis upper bound for the coverage axis
    """
    if config.show_log:  # correct log layout
        fig.update_yaxes(type="log", dtick=1, row=1, range=[0, log_upper], autorange=False)
    fig.update_layout(updatemenus=[{"visible": False}])  # no buttons
    fig.update_layout(annotations=[{"visible": False}])  # no annotations


def create_coverage_plot(fig, row, coverage_df, bin_size):
    """
    :param fig: plotly fig
    :param row: where to plot
    :param coverage_df: coverage dataframe
    :param bin_size: bin size for the coverage calculation
    :return: updated figure
    """
    # average data if there is a bin size > 1
    if bin_size > 1:
        positions, coverage, a_count, c_count, g_count, t_count = [], [], [], [], [], []
        for pos in coverage_df["position"][::bin_size]:
            if pos == 1:
                last_index = pos-1
                continue
            # get the index of the prior row (pos is one more than index)
            index = pos - 2
            positions.append(pos-1)
            coverage.append(round(coverage_df.loc[last_index:index, "coverage"].mean(), 0))
            a_count.append(round(coverage_df.loc[last_index:index, "A"].mean(), 2))
            c_count.append(round(coverage_df.loc[last_index:index, "C"].mean(), 2))
            g_count.append(round(coverage_df.loc[last_index:index, "G"].mean(), 2))
            t_count.append(round(coverage_df.loc[last_index:index, "T"].mean(), 2))
            # remember the index for the next bin start
            last_index = pos - 1
        # create new df for cov plot
        coverage_df_plot = pd.DataFrame(
            list(zip(positions, coverage, a_count, c_count, g_count, t_count)),
            columns=["position", "coverage", "A", "C", "G", "T"]
        )
    elif bin_size == 1:
        coverage_df_plot = coverage_df
    else:
        logger.error("bin size below 1 is not valid")
        sys.exit(1)

    # define hover template
    h_template = ""
    for index, description in enumerate(["position", "coverage", "percentage A", "percentage C", "percentage G", "percentage T"]):
        h_template = h_template + f"<b>{description}: </b>%" + "{customdata" + f"[{index}" + "]}<br>"
    h_template = h_template + "<extra></extra>"  # remove trace name
    # add dots with info
    fig.add_trace(
        go.Scatter(
            x=coverage_df_plot["position"],
            y=coverage_df_plot["coverage"],
            customdata=coverage_df_plot,
            fill="tonexty",
            fillcolor=config.coverage_fill_color,
            line={"color": config.coverage_line_color},
            hovertemplate=h_template,
            legendgroup="coverage",
            name="coverage",
            showlegend=True
        ),
        row=row,
        col=1
    )
    # add average info
    average_cov = statistics.mean(coverage_df["coverage"])
    fig.add_trace(
        go.Scatter(
            x=[min(coverage_df["position"]), max(coverage_df["position"])],
            y=[average_cov]*2,
            text=["", f"{round(average_cov)}x"],
            textposition="bottom left",
            mode="lines+text",
            line={"color": config.average_line_color, "width": config.average_line_width, "dash": "dash"},
            showlegend=True,
            legendgroup="average",
            name="average",
        ),
        row=row,
        col=1
    )


def split_vcf_df(df):
    """
    splits the vcf dataframe into multiple dfs if multiple mutations are on the same pos
    :param df: vcf_df
    :return: list of individual dfs
    """
    list_index = 0
    unique, counter, max_n = list(Counter(df["position"]).keys()), [0]*len(set(df["position"])), max(Counter(df["position"]).values())

    if max_n > 1:
        sub_dfs = [[]]*max_n

        for row in df.values.tolist():
            for index, pos in enumerate(unique):
                if row[0] == pos:
                    list_index = counter[index]
                    counter[index] += 1
            sub_dfs[list_index] = sub_dfs[list_index] + [row]

        return [pd.DataFrame(x, columns=df.columns) for x in sub_dfs]
    else:
        return [df]


def adjust_array_min_distance(
    values: list[float],
    min_distance: float,
    max_iterations: int = 1000,
    tolerance: float = 1e-6,
) -> list[float]:
    """
    Adjust 1D values to keep a minimum distance with minimal deviation.

    Ported from ResistanceProfiler's implementation: an iterative algorithm
    that walks *sorted* values and separates neighbouring entries when they
    are closer than ``min_distance``. The array is re-sorted every iteration
    because value mutations from previous passes can change the ordering.
    Convergence is detected when the largest per-iteration adjustment falls
    below ``tolerance``.

    :param values: values to adjust
    :param min_distance: required minimum distance between neighbouring values
    :param max_iterations: upper bound on refinement iterations
    :param tolerance: convergence threshold for the max per-iteration adjustment
    :return: adjusted values with the same ordering as the input
    """
    arr = list(values)
    n = len(arr)
    if n == 0:
        return []
    if n == 1:
        return [float(arr[0])]

    for _ in range(max_iterations):
        # re-sort each iteration: value mutations from previous passes change
        # the ordering, so the neighbour pairs to check change too.
        sorted_indices = sorted(range(n), key=lambda i: arr[i])
        max_adjustment = 0.0
        for i in range(1, n):
            idx1 = sorted_indices[i - 1]
            idx2 = sorted_indices[i]
            gap = arr[idx2] - arr[idx1]
            if gap < min_distance - tolerance:
                push = (min_distance - gap) / 2.0
                arr[idx1] -= push
                arr[idx2] += push
                max_adjustment = max(max_adjustment, push)
        if max_adjustment < tolerance:
            break

    return [float(v) for v in arr]


def create_vcf_plot(fig, row, vcf_df):
    """
    :param fig: plotly fig
    :param row: where to plot
    :param vcf_df: df with vcf info
    :return: updated figure
    """
    # disable false positive slice warning
    pd.options.mode.chained_assignment = None

    # draw stems from x value to jittered x value
    # get max x-value to calculate min distance to jitter
    for trace in fig['data']:
        if 'x' in trace:
            x_values = trace['x']
            if x_values is not None:
                max_x = max(x_values)

    # Jitter the top-dot x positions so overlapping variants are dispersed.
    # Ported from ResistanceProfiler's _apply_top_jitter: sort by (position,
    # allele frequency) so ties break deterministically, then enforce a minimum
    # separation of max_x / 200 — the same 1/200-of-the-visible-range ratio the
    # reference uses (feature_length / 200), scaled to the coverage x-axis. The
    # stem stays anchored at the true genomic position; only the top dot shifts.
    sorted_indices = sorted(
        range(len(vcf_df)),
        key=lambda i: (vcf_df["position"].iloc[i],
                       vcf_df["AF"].iloc[i] if "AF" in vcf_df else 0),
    )
    sorted_positions = [vcf_df["position"].iloc[i] for i in sorted_indices]
    min_distance = max_x / 200
    jittered_sorted = adjust_array_min_distance(
        sorted_positions, min_distance=min_distance)

    # scatter back into original row order
    jittered = [0.0] * len(vcf_df)
    for orig_idx, jit_x in zip(sorted_indices, jittered_sorted):
        jittered[orig_idx] = jit_x
    vcf_df["position_jittered"] = jittered

    # add stem independent of upper layers
    if "AF" in vcf_df:
        y_data = vcf_df["AF"]
    else:
        y_data = [1] * len(vcf_df["position"])
    # add lines
    shapes = []
    for x_value, x_value_jittered, y_value in zip(vcf_df["position"], vcf_df["position_jittered"], y_data):
        for coordinates in [(x_value, -0.3, x_value, -0.15),
                            (x_value, -0.15, x_value_jittered, 0),
                            (x_value_jittered, 0, x_value_jittered, y_value)]:
            shapes.append(
                {
                    "type": "line",
                    "xref": f"x{row}",
                    "yref": f"y{row}",
                    "x0": coordinates[0],
                    "y0": coordinates[1],
                    "x1": coordinates[2],
                    "y1": coordinates[3],
                    "line": {
                        "color": config.stem_color,
                        "width": config.stem_width
                    },
                    "layer": 'below'
                }
            )

    # plot shape in each subplot
    if fig["layout"]["shapes"]:
        for shape in shapes:
            fig.add_shape(shape)
    else:
        fig.update_layout(shapes=shapes)

    # plot the respective jittered x values as scatter
    for mut, color in zip(["SNP", "INS", "DEL"], [config.snp_color, config.ins_color, config.del_color]):
        if mut not in list(vcf_df["type"]):
            continue
        # plot a single layer for each mutation type and again a layer
        # if there are multiple mutations of the same type at the same position
        vcf_subset_temp = vcf_df[vcf_df["type"] == mut]
        vcf_subsets = split_vcf_df(vcf_subset_temp)
        for vcf_subset, show_legend in zip(vcf_subsets, [True]+[False]*(len(vcf_subsets)-1)):
            # check if allelic frequency is present otherwise set to one
            if "AF" in vcf_subset:
                vcf_subset["AF"] = vcf_subset["AF"].round(2)
                y_data = vcf_subset["AF"]
            else:
                y_data = [1] * len(vcf_subset["position"])
            # create hover template
            h_template = ""
            for index, description in enumerate(list(vcf_subset.columns)):
                # do not show position
                if index == len(vcf_subset.columns)-1:  # last one excludes the jittered data
                    continue
                h_template = h_template + f"<b>{description}: </b>%" + "{customdata" + f"[{index}" + "]}<br>"
            h_template = h_template + "<extra></extra>"  # remove trace name

            # add trace
            fig.add_trace(
                go.Scatter(
                    x=vcf_subset["position_jittered"],
                    y=y_data,
                    name=mut,
                    legendgroup=mut,
                    mode="markers",
                    customdata=vcf_subset,
                    showlegend=show_legend,
                    hovertemplate=h_template,
                    marker={
                        "color": color,
                        "size": config.variant_marker_size,
                        "line": {
                            "color": config.variant_line_color,
                            "width": config.variant_marker_line_width
                        }
                    }
                ),
            row=row,
            col=1
            )

    # update the y axis
    fig.update_yaxes(
        range=[-0.3, 1.15],
        tickvals=[0, 0.5, 1],
        ticktext=['0', '0.5', '1'],
        col=1,
        row=row
    )

    # not need to show yaxis if af is not in vcf
    if "AF" not in vcf_df:
        fig.update_yaxes(visible=False, row=row, col=1)
    else:
        fig.update_yaxes(title_text="frequency", row=row, col=1)


def create_track_plot(fig, row, feature_dict, box_size, box_alpha):
    """
    :param fig: plotly fig
    :param row: where to plot
    :param feature_dict: all infos for tracks as dictionary
    :param box_size: list of box sizes
    :param box_alpha: list of box alpha values
    :param subplot: subplot index of the plot
    :return: updated figure
    """
    # define colors
    n_colors = len(feature_dict)
    if n_colors > 1:
        colors = px.colors.sample_colorscale(config.track_color_scheme, [n / (n_colors - 1) for n in range(n_colors)])
    else:
        colors = [config.track_color_single]

    for feature, color in zip(feature_dict, colors):
        # define colors with 2 different alpha values, box size and cycle counter
        color_thes, b_size, cycle = ["rgba(" + color[4:-1] + f", {box_alpha[0]})", "rgba(" + color[4:-1] + f", {box_alpha[1]})"], box_size, 0
        # iterate over the different seq features
        for annotation, legend_vis in zip(feature_dict[feature], [True] + [False] * (len(feature_dict[feature]) - 1)):
            # define current cycle
            if cycle == 2:
                cycle = 0
            # get various plot info
            track = feature_dict[feature][annotation]["track"]
            # define strand marker
            if feature_dict[feature][annotation]["strand"] == "+":
                marker_type = config.strand_types[0]
            elif feature_dict[feature][annotation]["strand"] == "-":
                marker_type = config.strand_types[1]
            else:
                marker_type = config.strand_types[2]
            # define a hover text
            h_text = f"<b>type: </b> {feature}<br>"
            for classifier in feature_dict[feature][annotation]:
                if classifier == "track" or classifier == "translation":
                    continue
                h_text = h_text + f"<b>{classifier}: </b>{feature_dict[feature][annotation][classifier]} <br>"

            # define place for hover info
            global_start, global_stop = min(feature_dict[feature][annotation]["start"]), max(feature_dict[feature][annotation]["stop"])
            x = global_start + (global_stop - global_start) / 2
            # add hover info
            fig.add_trace(
                go.Scatter(
                    x=[x],
                    y=[track],
                    legendgroup=feature,
                    mode="markers",
                    marker={
                        "size": config.strand_marker_size,
                        "symbol": marker_type,
                        "color": color,
                        "line": {
                            "width": config.strand_marker_line_width,
                            "color": config.strand_marker_line_color
                        }
                    },
                    name="",
                    showlegend=False,
                    hoverinfo="text",
                    hovertext=h_text
                ),
                row=row,
                col=1
            )
            # plot the features
            previous_loc, previous_legend_vis = None, None
            # plot parts separately
            for start, stop in zip(feature_dict[feature][annotation]["start"], feature_dict[feature][annotation]["stop"]):
                # plot a line to indicate that the feature parts belong together
                if previous_loc is not None:
                    if feature_dict[feature][annotation]["strand"] == '-':
                        line_x = [min(previous_loc), max(start, stop)]
                    else:
                        line_x = [max(previous_loc), min(start, stop)]
                    fig.add_trace(
                        go.Scatter(
                            x=[line_x[0], line_x[1]],
                            y=[track, track],  # y-coordinate of the hline
                            mode='lines',
                            line={"color": color_thes[cycle]},
                            legendgroup=feature,
                            hoverinfo='skip',
                            showlegend=False,
                        ),
                        row=row,
                        col=1
                    )
                # plot the feature
                fig.add_trace(
                    go.Scatter(
                        x=[start, stop, stop, start, start],
                        y=[track + b_size[cycle], track + b_size[cycle], track - b_size[cycle], track - b_size[cycle], track + b_size[cycle]],
                        mode="lines",
                        fill="toself",
                        fillcolor=color_thes[cycle],
                        line={"color": color_thes[cycle]},
                        showlegend=legend_vis if previous_legend_vis is None else False,  # edge case for tracks that start with part features and result in legend duplication
                        hoverinfo='skip',
                        name=feature,
                        legendgroup=feature,
                    ),
                    row=row,
                    col=1
                )
                previous_loc, previous_legend_vis = [start, stop], legend_vis
            # switch to next cycle
            cycle += 1
    # reverse yaxis and hide it
    fig.update_yaxes(visible=False, row=row, col=1, autorange="reversed")


