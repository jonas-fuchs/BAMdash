"""
contains defs for plotting
"""

# BUILT-INS
import statistics
# LIBS
import plotly.graph_objects as go

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
            y=[average_cov, average_cov],
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

