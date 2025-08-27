import marimo

__generated_with = "0.15.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import polars as pl
    import altair as alt
    from pathlib import Path

    # read prepared DE data
    results_directory = Path("./")
    df = (
        pl.read_csv(results_directory / "turboid_differential_expression.csv")
        .with_columns(pl.col("time_point").fill_null(" "))
        .with_columns(
            pl.concat_str(
                pl.col("Cre"),
                pl.col("sample_type"),
                pl.col("comparison"),
                pl.col("time_point"),
                separator=" ",
            )
            .str.replace_all(".", "", literal=True)
            .alias("comparison_concat")
        )
        .filter(pl.col("comparison").str.contains(pattern="cre").not_())
    )


    x = mo.md("##The organ-specific secretomes of negative and positive energy balance")
    x
    return alt, df, mo, pl


@app.cell
def _(df, mo):
    experiment_1 = mo.ui.dropdown(
        options=sorted(list(set(df["comparison_concat"]))),
        label="Choose your first comparison (x axis)",
        value=sorted(list(set(df["comparison_concat"])))[0],
        searchable=True,
    )

    experiment_2 = mo.ui.dropdown(
        options=sorted(list(set(df["comparison_concat"]))),
        label="Choose your second comparison (y axis)",
        value=sorted(list(set(df["comparison_concat"])))[1],
        searchable=True,
    )
    mo.vstack([experiment_1, experiment_2])
    return experiment_1, experiment_2


@app.cell
def _(alt, df, experiment_1, experiment_2, mo, pl):
    def scatter_2(df):
        x = f"log2_FC_{experiment_1.selected_key}"
        y = f"log2_FC_{experiment_2.selected_key}"

        # Filter out inf and NaN values
        valid_df = df.filter(pl.col(x).is_finite() & pl.col(y).is_finite())

        domain_min = min(valid_df[x].min(), valid_df[y].min())
        domain_max = max(valid_df[x].max(), valid_df[y].max()) * 1.1
        # Define the main plot domain and margin position
        margin_pos = domain_min * 1.15 

        # Create significance categories based on p-values
        x_pval = f"p_value_{experiment_1.selected_key}"
        y_pval = f"p_value_{experiment_2.selected_key}"

        # Calculate overlap count (points with both x and y values)
        overlap_count = df.filter(
            (~pl.col(x).is_null()) & (~pl.col(y).is_null())
        ).height

        # Calculate Pearson correlation coefficient for overlapping data
        overlap_data = df.filter(
            (~pl.col(x).is_null()) & (~pl.col(y).is_null())
        )

        if overlap_count > 1:
            # Calculate correlation using numpy
            import numpy as np
            x_vals = overlap_data[x].to_numpy()
            y_vals = overlap_data[y].to_numpy()
            correlation = np.corrcoef(x_vals, y_vals)[0, 1]
        else:
            correlation = None

        # Identify top 10 points based on absolute fold change (x or y)
        df_with_abs_fc = df.with_columns([
            pl.coalesce([pl.col(x).abs(), pl.lit(0)]).alias("abs_x_fc"),
            pl.coalesce([pl.col(y).abs(), pl.lit(0)]).alias("abs_y_fc"),
        ]).with_columns([
            pl.max_horizontal(["abs_x_fc", "abs_y_fc"]).alias("max_abs_fc")
        ])

        # Get top 10 indices based on maximum absolute fold change
        top_10_indices = (
            df_with_abs_fc
            .with_row_index()
            .sort("max_abs_fc", descending=True)
            .head(10)["index"]
            .to_list()
        )

        # Create data with margin positions for missing values and significance categories
        df_with_margins = df.with_columns(
            [
                pl.when(pl.col(x).is_null())
                .then(margin_pos)
                .otherwise(pl.col(x))
                .alias(f"{x}_display"),
                pl.when(pl.col(y).is_null())
                .then(margin_pos)
                .otherwise(pl.col(y))
                .alias(f"{y}_display"),
                # Create significance category
                pl.when((pl.col(x_pval) < 0.05) & (pl.col(y_pval) < 0.05))
                .then(pl.lit("Both significant"))
                .when(pl.col(x_pval) < 0.05)
                .then(pl.lit("X-axis significant"))
                .when(pl.col(y_pval) < 0.05)
                .then(pl.lit("Y-axis significant"))
                .otherwise(pl.lit("Neither significant"))
                .alias("significance_category"),
                # Add label for top 10 points
                pl.col("entry_name").str.split("_").list.first().alias("short_name"),
            ]
        ).with_row_index().with_columns([
            pl.when(pl.col("index").is_in(top_10_indices))
            .then(pl.col("short_name"))
            .otherwise(pl.lit(""))
            .alias("label")
        ])

        # Calculate the actual plot domain
        plot_domain_min = domain_min * 1.3
        plot_domain_max = domain_max

        # Main scatter plot with extended domain and no grid lines
        scatter_plot = (
            alt.Chart(df_with_margins)
            .mark_circle()
            .encode(
                x=alt.X(
                    f"{x}_display:Q",
                    scale=alt.Scale(domain=(plot_domain_min, plot_domain_max)),
                    axis=alt.Axis(grid=False),
                ),
                y=alt.Y(
                    f"{y}_display:Q",
                    scale=alt.Scale(domain=(plot_domain_min, plot_domain_max)),
                    axis=alt.Axis(grid=False),
                ),
                color=alt.Color(
                    "significance_category:N",
                    scale=alt.Scale(
                        domain=[
                            "Neither significant",
                            "X-axis significant",
                            "Y-axis significant",
                            "Both significant",
                        ],
                        range=[
                            "#d3d3d3",
                            "#d6e5b3",
                            "#75c489",
                            "#5d8f70",
                        ],
                    ),
                    legend=alt.Legend(title="Significance (p < 0.05)"),
                ),
                tooltip=[
                    "entry_name:N",
                    "uniprot_id:N",
                    "Mass:Q",
                    "uniprot_function:N",
                ],
            )
        )

        # Origin lines (x=0 and y=0)
        x_origin_line = (
            alt.Chart(pl.DataFrame({"x": [0]}))
            .mark_rule(color="gray", strokeWidth=1)
            .encode(
                x=alt.X(
                    "x:Q",
                    scale=alt.Scale(domain=(plot_domain_min, plot_domain_max)),
                )
            )
        )

        y_origin_line = (
            alt.Chart(pl.DataFrame({"y": [0]}))
            .mark_rule(color="gray", strokeWidth=1)
            .encode(
                y=alt.Y(
                    "y:Q",
                    scale=alt.Scale(domain=(plot_domain_min, plot_domain_max)),
                )
            )
        )

        # Vertical line separating x margin 
        x_margin_line = (
            alt.Chart(pl.DataFrame({"x": [domain_min]}))
            .mark_rule(color="black", strokeWidth=1)
            .encode(
                x=alt.X(
                    "x:Q",
                    scale=alt.Scale(domain=(plot_domain_min, plot_domain_max)),
                )
            )
        )

        # Horizontal line separating y margin 
        y_margin_line = (
            alt.Chart(pl.DataFrame({"y": [domain_min]}))
            .mark_rule(color="black", strokeWidth=1)
            .encode(
                y=alt.Y(
                    "y:Q",
                    scale=alt.Scale(domain=(plot_domain_min, plot_domain_max)),
                )
            )
        )

        # Create text annotations for correlation and overlap count
        annotation_data = pl.DataFrame({
            "x": [plot_domain_max * 0.95, plot_domain_max * 0.95],
            "y": [plot_domain_max * 0.95, plot_domain_max * 0.85],
            "text": [
                f"r = {correlation:.3f}" if correlation is not None else "r = N/A",
                f"n = {overlap_count}",
            ]
        })

        text_annotations = (
            alt.Chart(annotation_data)
            .mark_text(
                align="right",
                baseline="top",
                fontSize=12,
                color="black"
            )
            .encode(
                x=alt.X(
                    "x:Q",
                    scale=alt.Scale(domain=(plot_domain_min, plot_domain_max)),
                ),
                y=alt.Y(
                    "y:Q",
                    scale=alt.Scale(domain=(plot_domain_min, plot_domain_max)),
                ),
                text="text:N"
            )
        )

        # Add labels for top 10 points
        point_labels = (
            alt.Chart(df_with_margins.filter(pl.col("label") != ""))
            .mark_text(
                align="left",
                baseline="middle",
                dx=5,  # offset to the right of the point
                fontSize=10,
                color="black"
            )
            .encode(
                x=alt.X(
                    f"{x}_display:Q",
                    scale=alt.Scale(domain=(plot_domain_min, plot_domain_max)),
                ),
                y=alt.Y(
                    f"{y}_display:Q",
                    scale=alt.Scale(domain=(plot_domain_min, plot_domain_max)),
                ),
                text="label:N"
            )
        )

        return (
            x_origin_line
            + y_origin_line
            + scatter_plot
            + x_margin_line
            + y_margin_line
            + text_annotations
            + point_labels
        ).properties(width=500, height=500)


    plot_data_2 = (
        df.filter(
            pl.col("comparison_concat").str.contains_any(
                patterns=[experiment_1.selected_key, experiment_2.selected_key]
            )
        )
        .drop(
            [
                "Regulation",
                "-log10_pval",
                "number_of_proteins",
                "experiment",
                "comparison",
                "Cre",
                "time_point",
                "sample_type",
                "pep_num",
            ]
        )
        .pivot(on="comparison_concat", values=["log2_FC", "p_value"])
    )

    chart = mo.ui.altair_chart(scatter_2(plot_data_2))
    chart
    return (plot_data_2,)


@app.cell
def _(plot_data_2):
    plot_data_2
    return


if __name__ == "__main__":
    app.run()
