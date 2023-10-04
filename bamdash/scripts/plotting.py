"""
contains defs for plotting
"""

# BUILT-INS
import statistics
# LIBS
import plotly.graph_objects as go
import plotly.express as px
# BAMDASH
from bamdash.scripts import config


def create_coverage_plot(fig, row, coverage_df):
    """
    :param fig: plotly fig
    :param row: where to plot
    :param coverage_df: coverage dataframe
    :return: updated figure
    """

    # define hover template
    h_template = ""
    for index, description in enumerate(["coverage", "percentage A", "percentage C", "percentage G", "percentage T"]):
        h_template = h_template + f"<b>{description}: </b>%" + "{customdata" + f"[{index+1}" + "]}<br>"
    h_template = h_template + "<extra></extra>"  # remove trace name
    # add dots with info
    fig.add_trace(
        go.Scatter(
            x=coverage_df["position"],
            y=coverage_df["coverage"],
            customdata=coverage_df,
            fill="tonexty",
            fillcolor=config.coverage_fill_color,
            line=dict(color=config.coverage_line_color),
            hovertemplate=h_template,
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
            textposition="top left",
            mode="lines+text",
            line=dict(color=config.average_line_color, width=config.average_line_width, dash="dash"),
            showlegend=True,
            name="average"
        ),
        row=row,
        col=1
    )
    # y axis title
    fig.update_yaxes(title_text="genome coverage", range=[0, max(coverage_df["coverage"])], row=row, col=1)
    # Add dropdown
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="left",
                buttons=list([
                    dict(
                        args=[f"yaxis{row}.type" if row > 1 else f"yaxis.type", "linear"],
                        label="linear",
                        method="relayout"
                    ),
                    dict(
                        args=[f"yaxis{row}.type" if row > 1 else f"yaxis.type", "log"],
                        label="log",
                        method="relayout"
                    )
                ]),
                pad={"r": 10, "t": 10},
                showactive=True,
                xanchor="left",
                y=1.1,
                yanchor="top"
            ),
        ],
    )


def create_vcf_plot(fig, row, vcf_df):
    """
    :param fig: plotly fig
    :param row: where to plot
    :param vcf_df: df with vcf info
    :return: updated figure
    """
    # check if allelic frequency is present,
    # otherwise set to one
    if "AF" in vcf_df:
        y_data = vcf_df["AF"]
    else:
        y_data = [1] * len(vcf_df["position"])
    # create hover template
    h_template = ""
    for index, description in enumerate(list(vcf_df.columns)):
        # do not show position
        if index == 0:
            continue
        h_template = h_template + f"<b>{description}: </b>%" + "{customdata" + f"[{index}" + "]}<br>"
    h_template = h_template + "<extra></extra>"  # remove trace name

    for mut, color in zip(["SNP", "INS", "DEL"], [config.snp_color, config.ins_color, config.del_color]):
        if mut not in list(vcf_df["type"]):
            continue
        # add trace
        fig.add_trace(
            go.Scatter(
                x=vcf_df["position"],
                y=y_data,
                name=mut,
                mode="markers",
                customdata=vcf_df,
                showlegend=True,
                hovertemplate=h_template,
                marker=dict(
                    color=color,
                    size=config.variant_marker_size,
                    line=dict(
                        color=config.variant_line_color,
                        width=config.variant_marker_line_width
                    )
                )
            ),
        row=row,
        col=1
        )
    # add lines
    fig.update_layout(
        shapes=[dict(
            type="line",
            xref=f"x{row}",
            yref=f"y{row}",
            x0=x,
            y0=0,
            x1=x,
            y1=y-0.05,
            line=dict(
                color=config.stem_color,
                width=config.stem_width
            )
        ) for x, y in zip(vcf_df["position"], y_data)]
    )
    # not need to show yaxis if af is not in vcf
    if "AF" not in vcf_df:
        fig.update_yaxes(visible=False, row=row, col=1)
    else:
        fig.update_yaxes(title_text="frequency", range=[0, 1], row=row, col=1)


def create_track_plot(fig, row, feature_dict, box_size, box_alpha):
    """
    :param fig: plotly fig
    :param row: where to plot
    :param feature_dict: all infos for tracks as dictionary
    :param box_size: list of box sizes
    :param box_alpha: list of box alpha values
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
            positions = [int(x) for x in annotation.split(" ")]
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
            x = positions[0] + (positions[1] - positions[0]) / 2
            # add hover info
            fig.add_trace(
                go.Scatter(
                    x=[x],
                    y=[track],
                    legendgroup=feature,
                    mode="markers",
                    marker=dict(
                        size=config.strand_marker_size,
                        symbol=marker_type,
                        color=color,
                        line=dict(
                            width=config.strand_marker_line_width,
                            color=config.strand_marker_line_color
                        )
                    ),
                    name="",
                    showlegend=False,
                    hoverinfo="text",
                    hovertext=h_text
                ),
                row=row,
                col=1
            )
            # add the track rectangle
            fig.add_trace(
                go.Scatter(
                    x=[positions[0], positions[1], positions[1], positions[0], positions[0]],
                    y=[track + b_size[cycle], track + b_size[cycle], track - b_size[cycle], track - b_size[cycle], track + b_size[cycle]],
                    mode="lines",
                    fill="toself",
                    fillcolor=color_thes[cycle],
                    line=dict(color=color_thes[cycle]),
                    showlegend=legend_vis,
                    hoverinfo='skip',
                    name="",
                    legendgroup=feature,
                    legendgrouptitle_text=feature,
                ),
                row=row,
                col=1
            )
            # switch to next cycle
            cycle += 1
    # reverse yaxis and hide it
    fig.update_yaxes(visible=False, row=row, col=1, autorange="reversed")


