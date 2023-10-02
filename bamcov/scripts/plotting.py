"""
contains defs for plotting
"""

# BUILT-INS
import statistics
# LIBS
import plotly.graph_objects as go
import plotly.express as px

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
            fillcolor="rgba(255, 212, 135, 0.2)",
            line=dict(color="rgba(224, 168, 68, 1)"),
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
            text=[f"{round(average_cov)}x", ""],
            textposition="top right",
            mode="lines+text",
            line=dict(color="grey", width=1, dash="dash"),
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

    for mut, color in zip(["SNP", "INS", "DEL"], ["grey", "blue", "red"]):
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
                    size=14,
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
                color="grey",
                width=1
            )
        ) for x, y in zip(vcf_df["position"], y_data)]
    )
    # not need to show y axis if af is not in vcf
    if "AF" not in vcf_df:
        fig.update_yaxes(visible=False, row=row, col=1)
    else:
        fig.update_yaxes(title_text="frequency", range=[0, 1], row=row, col=1)


def create_gb_plot(fig, row, feature_dict):
    """
    :param fig: plotly fig
    :param row: where to plot
    :param feature_dict: all infos from gb file
    :return: updated figure
    """
    # define colors
    n_colors = len(feature_dict)
    colors = px.colors.sample_colorscale("agsunset", [n / (n_colors - 1) for n in range(n_colors)])

    for feature, color in zip(feature_dict, colors):
        # define colors with 2 different alpha values, box size and cycle counter
        color_thes, b_size, cycle = ["rgba(" + color[4:-1] + ", 0.6)", "rgba(" + color[4:-1] + ", 0.8)"], [0.4,
                                                                                                           0.3], 0
        # iterate over the different seq features
        for annotation, legend_vis in zip(feature_dict[feature], [True] + [False] * (len(feature_dict) - 1)):
            # define current cycle
            if cycle == 2:
                cycle = 0
            # get various plot info
            positions = [int(x) for x in annotation.split(" ")]
            track = feature_dict[feature][annotation]["track"]
            single_pos = list(range(positions[0], positions[1] + 1))
            # define a hover template
            h_template = f"<b>type: </b> {feature}<br>"
            for index, classifier in enumerate(feature_dict[feature][annotation]):
                if classifier == "track" or classifier == "translation":
                    continue
                h_template = h_template + f"<b>{classifier}: </b>%" + "{customdata" + f"[{index}" + "]}<br>"
            h_template = h_template + "<extra></extra>"
            # and create the custom data
            custom_data = [list(feature_dict[feature][annotation].values())[:-1]] * len(single_pos)
            # add upper line
            fig.add_trace(
                go.Scatter(
                    x=positions,
                    y=[track + b_size[cycle], track + b_size[cycle]],
                    mode="lines",
                    line=dict(color="grey"),
                    showlegend=False,
                    legendgroup=feature,
                    hoverinfo="skip"
                ),
                row=row,
                col=1
            )
            # fill the square
            fig.add_trace(
                go.Scatter(
                    x=positions,
                    y=[track - b_size[cycle], track - b_size[cycle]],
                    mode="none",
                    line=dict(color="grey"),
                    fill="tonexty",
                    fillcolor=color_thes[cycle],
                    showlegend=legend_vis,
                    name="",
                    legendgroup=feature,
                    legendgrouptitle_text=feature,
                    hoverinfo="skip"
                ),
                row=row,
                col=1
            )
            # plot remaining lines to form a square
            for x, y in zip(
                    [positions, [positions[0]] * 2, [positions[1]] * 2],
                    [[track - b_size[cycle], track - b_size[cycle]], [track + b_size[cycle], track - b_size[cycle]],
                     [track + b_size[cycle], track - b_size[cycle]]]
            ):
                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        mode="lines",
                        line=dict(color="grey"),
                        showlegend=False,
                        legendgroup=feature,
                        hoverinfo="skip"
                    ),
                    row=row,
                    col=1
                )
            # hover info
            fig.add_trace(
                go.Scatter(
                    x=single_pos,
                    y=[track] * len(single_pos),
                    mode=None,
                    line=dict(color=color),
                    opacity=0,
                    showlegend=False,
                    legendgroup=feature,
                    hovertemplate=h_template,
                    customdata=custom_data
                ),
                row=row,
                col=1
            )
            fig.update_yaxes(visible=False, row=row, col=1)
            # switch to next cycle
            cycle += 1

