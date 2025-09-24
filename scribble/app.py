"""
Interactive visualization for multi-code count tables using Dash and Plotly.

This small web‑app reads a tab separated text file containing count data
indexed by one or more ``code_#`` columns along with a final ``count``
column.  The number of ``code_#`` columns may vary from file to file.

Users can select which codes should be plotted along the x, y and
optionally z axes.  When there are more codes in the table than axes
selected, the counts are marginalised over the unused codes (i.e. the
counts are summed across the unused dimensions).  Two‑dimensional
selections are visualised with a scatter plot whose bubble size and
colour reflect the total count.  Three‑dimensional selections are
visualised as a 3D scatter where marker size and colour encode the
marginal counts.

To run this app install the required dependencies with:

    pip install dash==2.* pandas plotly

Then execute ``python app.py`` and open the URL printed in the console
in your web browser.  By default it reads ``counts.txt`` from the
same directory as this script; adjust ``DATA_FILE`` below if your
dataset is stored elsewhere.
"""

from __future__ import annotations

import os
import sys
from typing import List, Dict, Tuple

import pandas as pd
from dash import Dash, dcc, html, Output, Input, callback, ctx
import plotly.graph_objs as go
import plotly.express as px


# Path to the data file.  Adjust this if your counts file has a
# different name or location.
DATA_FILE = os.path.join("/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/selections/AG24_1/counts.txt")
DATA_FILE = os.path.join("/Users/adrianomartinelli/projects/delt/delt-core/scribble/random_counts.csv")


def load_counts(path: str) -> Tuple[pd.DataFrame, List[str]]:
    """Load the count table from a tab‑separated text file.

    Parameters
    ----------
    path : str
        Location of the counts file.  The file must have a header row
        containing one or more ``code_#`` columns followed by a final
        ``count`` column.  Columns are expected to be separated by
        tabs, but any whitespace separator is accepted.

    Returns
    -------
    tuple
        A tuple containing the DataFrame with the loaded counts and a
        list of the names of the ``code_#`` columns.  The DataFrame
        always contains a numeric ``count`` column.
    """
    # Attempt to parse the file assuming tab separation first.  If this
    # fails, fall back to whitespace delimitation.  Using pandas here
    # ensures that large files are loaded efficiently and numeric
    # columns remain numeric.
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        df = pd.read_csv(path, delim_whitespace=True)

    # Normalise column names (strip whitespace) and identify the
    # code columns.  A code column is any column whose name begins
    # with "code_" (case sensitive).  We enforce a lowercase count
    # column name for consistency.
    df.columns = [c.strip() for c in df.columns]
    code_columns = [col for col in df.columns if col.startswith("code_")]
    if not code_columns:
        raise ValueError(
            "No code columns were detected.  Column names must start with 'code_'."
        )
    if "count" not in df.columns:
        raise ValueError("No 'count' column detected.  The file must have a count column.")

    # Convert all code columns to numeric if possible.  This makes
    # sorting and plotting behave sensibly.  If conversion fails,
    # pandas will leave the data as object dtype (string) which is
    # still acceptable for categorisation.  We do not force numeric
    # conversion globally because codes might be alphanumeric.
    for col in code_columns:
        df[col] = pd.to_numeric(df[col], errors="ignore")

    # Ensure the count column is numeric for plotting and summing.
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0)

    return df, code_columns


def create_app(df: pd.DataFrame, code_columns: List[str]) -> Dash:
    """Create and configure the Dash application.

    Parameters
    ----------
    df : pandas.DataFrame
        The count table with code columns and a numeric 'count'.
    code_columns : list
        Names of the code columns in ``df``.

    Returns
    -------
    dash.Dash
        A Dash application ready to be run.
    """
    app = Dash(__name__)

    # Prepare options for the dropdowns.  Each option has a human
    # friendly label equal to the column name and a value equal to
    # that name.  For the z‑axis we append a 'None' option to allow
    # two‑dimensional plotting on datasets with three or more codes.
    dropdown_options = [
        {"label": col.replace("_", " ").title(), "value": col} for col in code_columns
    ]
    z_options = dropdown_options + [
        {"label": "None (2D)", "value": "__none__"}
    ]

    # When only one code is present we cannot build a meaningful plot.
    if len(code_columns) < 2:
        raise ValueError(
            "The dataset contains only one code column. At least two codes are required to plot."
        )

    # Define default selections: pick the first two codes for x and y.
    default_x = code_columns[0]
    default_y = code_columns[1] if len(code_columns) > 1 else code_columns[0]
    # If a third code exists, use it for z by default; otherwise
    # default to 2D plotting.
    default_z = code_columns[2] if len(code_columns) > 2 else "__none__"

    app.layout = html.Div(
        [
            html.H1("Code Count Visualiser", style={"textAlign": "center"}),
            html.Div(
                [
                    html.Label("X‑axis:"),
                    dcc.Dropdown(
                        id="x-axis-dropdown",
                        options=dropdown_options,
                        value=default_x,
                        clearable=False,
                        style={"width": "200px"},
                    ),
                ],
                style={"display": "inline-block", "margin": "0 20px"},
            ),
            html.Div(
                [
                    html.Label("Y‑axis:"),
                    dcc.Dropdown(
                        id="y-axis-dropdown",
                        options=dropdown_options,
                        value=default_y,
                        clearable=False,
                        style={"width": "200px"},
                    ),
                ],
                style={"display": "inline-block", "margin": "0 20px"},
            ),
            html.Div(
                [
                    html.Label("Z‑axis (optional):"),
                    dcc.Dropdown(
                        id="z-axis-dropdown",
                        options=z_options,
                        value=default_z,
                        clearable=False,
                        style={"width": "200px"},
                    ),
                ],
                style={"display": "inline-block", "margin": "0 20px"},
            ),
            dcc.Graph(id="count-graph", style={"height": "700px"}),
            html.Div(
                [
                    html.P(
                        "Select codes for each axis to explore the data. When more than "
                        "three codes are present in the dataset, choosing only two or three "
                        "axes will sum counts over the remaining codes.",
                        style={"marginTop": "20px"},
                    )
                ],
                style={"width": "80%", "margin": "auto", "textAlign": "center"},
            ),
        ]
    )

    @app.callback(
        Output("count-graph", "figure"),
        Input("x-axis-dropdown", "value"),
        Input("y-axis-dropdown", "value"),
        Input("z-axis-dropdown", "value"),
    )
    def update_graph(x_axis: str, y_axis: str, z_axis: str):
        """Update the plotted figure when the user selects different axes.

        Parameters
        ----------
        x_axis, y_axis, z_axis : str
            Column names selected for the x, y and z axes.  The special
            value ``"__none__"`` indicates that the z axis is not used.

        Returns
        -------
        plotly.graph_objs.Figure
            The updated figure reflecting the user's selections.
        """
        selected_axes = [axis for axis in [x_axis, y_axis] if axis]
        if z_axis != "__none__":
            selected_axes.append(z_axis)

        # Group by the selected axes and sum the counts over all other
        # codes.  This automatically marginalises over unused codes.
        grouped = (
            df.groupby(selected_axes)["count"]
            .sum()
            .reset_index()
        )

        # Determine the maximum count for marker scaling.  Avoid
        # division by zero by using max(max_count, 1).
        max_count = grouped["count"].max() or 1

        if z_axis == "__none__":
            # Two‑dimensional scatter: size and colour encode the counts.
            fig = px.scatter(
                grouped,
                x=x_axis,
                y=y_axis,
                size="count",
                color="count",
                hover_name="count",
                title=f"Counts grouped by {x_axis} and {y_axis}",
                labels={x_axis: x_axis.replace("_", " ").title(),
                        y_axis: y_axis.replace("_", " ").title(),
                        "count": "Count"},
            )
            # Provide a consistent scaling for marker sizes.  Plotly
            # automatically scales sizes but specifying a range helps
            # maintain readability across different datasets.
            fig.update_traces(
                marker=dict(sizemode="area", sizeref=(2. * max_count / (40 ** 2)), sizemin=4)
            )
            fig.update_layout(
                xaxis_title=x_axis.replace("_", " ").title(),
                yaxis_title=y_axis.replace("_", " ").title(),
                legend_title="Count",
                coloraxis_colorbar=dict(title="Count"),
            )
        else:
            # Three‑dimensional scatter.
            fig = go.Figure(
                data=[
                    go.Scatter3d(
                        x=grouped[x_axis],
                        y=grouped[y_axis],
                        z=grouped[z_axis],
                        mode="markers",
                        marker=dict(
                            size=grouped["count"] / max_count * 30 + 4,
                            color=grouped["count"],
                            colorscale="Viridis",
                            showscale=True,
                            colorbar=dict(title="Count"),
                        ),
                        text=[f"Count: {c}" for c in grouped["count"]],
                        hovertemplate=
                            f"{x_axis}: %{{x}}<br>{y_axis}: %{{y}}<br>{z_axis}: %{{z}}<br>Count: %{{text}}"
                    )
                ]
            )
            fig.update_layout(
                scene=dict(
                    xaxis_title=x_axis.replace("_", " ").title(),
                    yaxis_title=y_axis.replace("_", " ").title(),
                    zaxis_title=z_axis.replace("_", " ").title(),
                ),
                title=f"Counts grouped by {x_axis}, {y_axis} and {z_axis}",
            )

        return fig

    return app


def main():
    """Entry point to run the Dash server."""
    # Load data and create app.  If the file cannot be read or the
    # structure is invalid, raise a descriptive exception.
    df, code_columns = load_counts(DATA_FILE)
    app = create_app(df, code_columns)
    # Run the server.  Setting debug to True enables automatic
    # reloading on code changes during development.  Use host="0.0.0.0"
    # to listen on all interfaces; by default Dash listens on
    # localhost only.
    app.run_server(debug=False)


if __name__ == "__main__":
    main()

import numpy as np
import pandas as pd
random_counts = pd.DataFrame({
    'code_1': np.random.random(100),
    'code_2': np.random.random(100),
    'code_3': np.random.random(100),
    'count': np.random.random(100),
})
random_counts.to_csv("/Users/adrianomartinelli/projects/delt/delt-core/scribble/random_counts.csv",
                     sep='\t',
                     index=False)
