#!/usr/bin/env python3
"""Generate publication-ready revision figures as vector PDF and 600 dpi JPEG.

Figures 1 (workflow artwork) and 4 (dashboard screenshot) are intentionally
excluded. Figure 3 is exported as two standalone QC panels rather than a
manuscript composite.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw, rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D
import seaborn as sns
import yaml

from delt_hit.quality_control.plot_codon_hits import get_stats
from delt_hit.cli.library.api import close_figure, prepare_graph_bundle, visualize_reaction_graph
from generate_composite_figures import (
    build_enumeration_composite,
    build_properties_and_examples_composite,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "paper" / "revision" / "figures"
EXAMPLE_ROOT = PROJECT_ROOT / "supporting_material" / "experiments" / "example-single-display"
CAMPAIGN_ROOT = EXAMPLE_ROOT / "campaign"
VISUALIZATION_ROOT = CAMPAIGN_ROOT / "library" / "visualization"
PROPERTIES_ROOT = CAMPAIGN_ROOT / "library" / "properties" / "AG24_4_top_hits"
REPORT_ROOT = CAMPAIGN_ROOT / "demultiplex" / "cutadapt_output_files"
BENCHMARK_ROOT = PROJECT_ROOT / "benchmarks" / "demultiplex" / "tools"
FAVALLI_ROOT = PROJECT_ROOT / "supporting_material" / "experiments" / "favalli"

FULL_WIDTH_IN = 180 / 25.4
JPEG_DPI = 600
DATASET_PATTERN = re.compile(
    r"^synthetic_(?P<cycles>\d+)cycle_(?P<bbpc>\d+)bbpc_(?P<depth>\d+)m(?:_err=\d+)?$"
)
TOOL_STYLES = {
    "deli": {"label": "DELi", "color": "#1f77b4", "marker": "o"},
    "delt": {"label": "DELT-Hit", "color": "#d62728", "marker": "s"},
}
METHOD_ORDER = ["counts", "edgeR", "z_score"]
METHOD_PALETTE = {"counts": "#4C78A8", "edgeR": "#F58518", "z_score": "#54A24B"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def configure_style() -> None:
    matplotlib.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 6,
            "axes.titlesize": 7,
            "axes.labelsize": 6,
            "xtick.labelsize": 5.5,
            "ytick.labelsize": 5.5,
            "legend.fontsize": 5.5,
            "legend.title_fontsize": 6,
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = output_dir / f"{stem}.pdf"
    jpeg_path = output_dir / f"{stem}.jpeg"
    facecolor = fig.get_facecolor()
    fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.03, facecolor=facecolor)
    fig.savefig(
        jpeg_path,
        dpi=JPEG_DPI,
        bbox_inches="tight",
        pad_inches=0.03,
        facecolor=facecolor,
        pil_kwargs={"quality": 95},
    )
    plt.close(fig)
    print(f"Wrote {pdf_path}")
    print(f"Wrote {jpeg_path}")


def load_report_rows() -> list[dict[str, float | int | str]]:
    rows = []
    for report_path in sorted(REPORT_ROOT.glob("*.cutadapt.json")):
        report = json.loads(report_path.read_text())
        input_reads = int(report["read_counts"]["input"])
        output_reads = int(report["read_counts"]["output"])
        rows.append(
            {
                "region": report_path.name.removesuffix(".cutadapt.json"),
                "input": input_reads,
                "output": output_reads,
                "discarded": input_reads - output_reads,
                "with_fraction": output_reads / input_reads,
                "discarded_fraction": 1 - output_reads / input_reads,
            }
        )
    if not rows:
        raise FileNotFoundError(f"No Cutadapt JSON reports found under {REPORT_ROOT}")
    return rows


def build_figure2(output_dir: Path) -> None:
    """Render the poor-selection report as a reproducible terminal-style panel."""
    rows = load_report_rows()
    fig, ax = plt.subplots(figsize=(FULL_WIDTH_IN, 4.25))
    background = "#202124"
    grid = "#b9bec5"
    foreground = "#e4e7eb"
    cyan = "#1fc7c9"
    danger = "#ff5a5f"
    fig.patch.set_facecolor(background)
    ax.set_facecolor(background)
    ax.axis("off")
    ax.text(
        0.5,
        0.965,
        "Cutadapt Pipeline Report",
        ha="center",
        va="top",
        color=foreground,
        family="monospace",
        fontsize=7,
        fontstyle="italic",
        transform=ax.transAxes,
    )

    headers = ["Region", "Input", "With adapters", "Discarded", "% with", "% discarded"]
    table_rows = [
        [
            str(row["region"]),
            f"{int(row['input']):,}",
            f"{int(row['output']):,}",
            f"{int(row['discarded']):,}",
            f"{float(row['with_fraction']):.2%}",
            f"{float(row['discarded_fraction']):.2%}",
        ]
        for row in rows
    ]
    table = ax.table(
        cellText=table_rows,
        colLabels=headers,
        cellLoc="right",
        colLoc="center",
        colWidths=[0.12, 0.18, 0.21, 0.17, 0.12, 0.18],
        bbox=[0.01, 0.27, 0.98, 0.60],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(6.2)
    for (row_index, column_index), cell in table.get_celld().items():
        cell.set_facecolor(background)
        cell.set_edgecolor(grid)
        cell.set_linewidth(0.8)
        cell.get_text().set_color(foreground)
        cell.get_text().set_fontfamily("monospace")
        if row_index == 0:
            cell.get_text().set_fontweight("bold")
        elif column_index == 0:
            cell.get_text().set_color(cyan)
            cell.get_text().set_ha("left")
        elif column_index == 5 and float(rows[row_index - 1]["discarded_fraction"]) > 0.10:
            cell.get_text().set_color(danger)
            cell.get_text().set_fontweight("bold")

    overall_input = int(rows[0]["input"])
    overall_output = int(rows[-1]["output"])
    overall_discarded = overall_input - overall_output
    ax.text(0.01, 0.18, "Overall:", color=foreground, family="monospace", weight="bold", transform=ax.transAxes)
    ax.text(0.04, 0.115, "With adapters :", color=foreground, family="monospace", transform=ax.transAxes)
    ax.text(
        0.28,
        0.115,
        f"{overall_output:,} ({overall_output / overall_input:.2%})",
        color=cyan,
        family="monospace",
        weight="bold",
        transform=ax.transAxes,
    )
    ax.text(0.04, 0.055, "Discarded     :", color=foreground, family="monospace", transform=ax.transAxes)
    ax.text(
        0.28,
        0.055,
        f"{overall_discarded:,} ({overall_discarded / overall_input:.2%})",
        color=cyan,
        family="monospace",
        weight="bold",
        transform=ax.transAxes,
    )
    save_figure(fig, output_dir, "figure2-qc-poor-selection-report")


def qc_plot_data(region_name: str) -> tuple[pd.DataFrame, int, int]:
    report_path = REPORT_ROOT / f"{region_name}.cutadapt.json"
    stats = get_stats(report_path)
    frame = pd.DataFrame(stats)
    plot_data = pd.concat(frame["error_counts"].tolist())
    plot_data = (
        plot_data[["index", "number_of_errors", "error_counts"]]
        .groupby(["index", "number_of_errors"])
        .sum()
        .reset_index()
        .pivot(index="index", columns="number_of_errors", values="error_counts")
        .fillna(0)
    )
    return plot_data, int(frame.reads_in.iloc[0]), int(frame.reads_out.iloc[0])


def build_qc_component(output_dir: Path, region_name: str, stem: str) -> None:
    plot_data, reads_in, reads_out = qc_plot_data(region_name)
    fig, axes = plt.subplots(1, 2, figsize=(FULL_WIDTH_IN, 3.45), constrained_layout=True)
    colors = ["#287bb5", "#f28e2b", "#59a14f"]
    plot_data.plot(kind="bar", stacked=True, ax=axes[0], width=0.86, color=colors[: len(plot_data.columns)])
    plot_data.plot(kind="bar", stacked=True, ax=axes[1], width=0.86, color=colors[: len(plot_data.columns)])
    axes[1].set_yscale("log")
    axes[1].set_ylim(1, max(axes[1].get_ylim()) * 1.25)
    for index, axis in enumerate(axes):
        axis.set_xlabel("Barcode index")
        axis.set_ylabel("Read count" if index == 0 else "Read count (log scale)")
        axis.set_xticks([])
        axis.grid(axis="y", color="#d9d9d9", linewidth=0.4, alpha=0.8)
        axis.legend(title="Number of errors", loc="upper right")
    save_figure(fig, output_dir, stem)


def build_qc_composite(output_dir: Path) -> None:
    """Build the manuscript-ready S0/C0 linear-count QC figure."""
    s0_data, _s0_reads_in, _s0_reads_out = qc_plot_data("0-S0")
    c0_data, _c0_reads_in, _c0_reads_out = qc_plot_data("1-C0")
    fig, axes = plt.subplots(1, 2, figsize=(FULL_WIDTH_IN, 3.45), constrained_layout=True)
    colors = ["#287bb5", "#f28e2b", "#59a14f"]

    for axis, plot_data, x_label in [
        (axes[0], s0_data, "S0 barcode index"),
        (axes[1], c0_data, "C0 index"),
    ]:
        plot_data.plot(
            kind="bar",
            stacked=True,
            ax=axis,
            width=0.86,
            color=colors[: len(plot_data.columns)],
        )
        axis.set_xlabel(x_label)
        axis.set_ylabel("Read count")
        axis.set_xticks([])
        axis.grid(axis="y", color="#d9d9d9", linewidth=0.4, alpha=0.8)
        axis.legend(title="Number of errors", loc="upper right")

    save_figure(fig, output_dir, "figure3-qc-s0-c0-composite")


def render_molecular_weight_source(output_path: Path) -> None:
    """Render the original histogram without its redundant Matplotlib title."""
    properties_path = PROPERTIES_ROOT / "properties.parquet"
    if not properties_path.exists():
        raise FileNotFoundError(properties_path)
    properties = pd.read_parquet(properties_path)
    fig, axis = plt.subplots(figsize=(8, 6), constrained_layout=True)
    sns.histplot(properties["prop_mw"].dropna(), kde=False, ax=axis, color="#4C8CB5")
    axis.set_xlabel("Molecular weight (Da)", fontsize=18)
    axis.set_ylabel("Frequency", fontsize=18)
    axis.tick_params(axis="both", labelsize=16)
    axis.yaxis.get_offset_text().set_fontsize(16)
    axis.grid(True)
    fig.savefig(output_path, dpi=600, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def save_original_composite(source_path: Path, output_dir: Path, stem: str) -> None:
    """Export an original fixed-pixel composite at 180 mm publication width."""
    image = Image.open(source_path).convert("RGB")
    effective_dpi = image.width / FULL_WIDTH_IN
    jpeg_path = output_dir / f"{stem}.jpeg"
    pdf_path = output_dir / f"{stem}.pdf"
    image.save(jpeg_path, format="JPEG", quality=95, dpi=(effective_dpi, effective_dpi))
    image.save(pdf_path, format="PDF", resolution=effective_dpi)
    print(f"Wrote {pdf_path}")
    print(f"Wrote {jpeg_path}")


def build_box_b0_code0(output_dir: Path) -> None:
    """Export the B0 code 0 structure used in the manuscript box."""
    config_path = CAMPAIGN_ROOT / "config.yaml"
    config = yaml.safe_load(config_path.read_text())
    entry = next(item for item in config["whitelists"]["B0"] if int(item["index"]) == 0)
    molecule = Chem.MolFromSmiles(entry["smiles"])
    if molecule is None:
        raise ValueError(f"Could not parse B0 code 0 SMILES: {entry['smiles']}")

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "box-b0-code0-structure.pdf"
    jpeg_path = output_dir / "box-b0-code0-structure.jpeg"
    png_path = output_dir / "box-b0-code0-structure.png"
    with tempfile.TemporaryDirectory(prefix="delt-hit-b0-code0-") as temporary_dir:
        svg_path = Path(temporary_dir) / "b0-code0.svg"
        drawer = rdMolDraw2D.MolDraw2DSVG(1200, 800)
        drawer.drawOptions().padding = 0.08
        drawer.drawOptions().bondLineWidth = 3.0
        drawer.drawOptions().minFontSize = 28
        drawer.drawOptions().clearBackground = False
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, molecule)
        drawer.FinishDrawing()
        svg_path.write_text(drawer.GetDrawingText())
        subprocess.run(
            [
                "rsvg-convert",
                "--format=pdf",
                "--width=510",
                f"--output={output_path}",
                str(svg_path),
            ],
            check=True,
        )
        subprocess.run(
            [
                "rsvg-convert",
                "--format=png",
                "--width=3600",
                f"--output={png_path}",
                str(svg_path),
            ],
            check=True,
        )
        transparent_image = Image.open(png_path).convert("RGBA")
        transparent_image.save(png_path, format="PNG", dpi=(600, 600))
        white_background = Image.new("RGBA", transparent_image.size, "white")
        white_background.alpha_composite(transparent_image)
        white_background.convert("RGB").save(
            jpeg_path,
            format="JPEG",
            quality=95,
            dpi=(600, 600),
        )
    print(f"Wrote {output_path}")
    print(f"Wrote {jpeg_path}")
    print(f"Wrote {png_path}")


def build_box3_reaction_template(output_dir: Path) -> None:
    """Export the Box 3 amide-bond-formation SMIRKS as reaction artwork."""
    config = yaml.safe_load((CAMPAIGN_ROOT / "config.yaml").read_text())
    smirks = config["catalog"]["reactions"]["ABF1"]["smirks"]
    reaction = rdChemReactions.ReactionFromSmarts(smirks)
    if reaction is None:
        raise ValueError(f"Could not parse Box 3 reaction SMIRKS: {smirks}")

    output_dir.mkdir(parents=True, exist_ok=True)
    stem = "box3-amide-bond-formation-reaction"
    pdf_path = output_dir / f"{stem}.pdf"
    png_path = output_dir / f"{stem}.png"
    jpeg_path = output_dir / f"{stem}.jpeg"
    draw_options = rdMolDraw2D.MolDrawOptions()
    draw_options.clearBackground = False
    draw_options.bondLineWidth = 3.0
    draw_options.minFontSize = 28

    with tempfile.TemporaryDirectory(prefix="delt-hit-box3-reaction-") as temporary_dir:
        svg_path = Path(temporary_dir) / f"{stem}.svg"
        svg = Draw.ReactionToImage(
            reaction,
            subImgSize=(900, 450),
            useSVG=True,
            drawOptions=draw_options,
        )
        svg_path.write_text(svg)
        subprocess.run(
            [
                "rsvg-convert",
                "--format=pdf",
                "--width=510",
                f"--output={pdf_path}",
                str(svg_path),
            ],
            check=True,
        )
        subprocess.run(
            [
                "rsvg-convert",
                "--format=png",
                "--width=3600",
                f"--output={png_path}",
                str(svg_path),
            ],
            check=True,
        )

    transparent_image = Image.open(png_path).convert("RGBA")
    transparent_image.save(png_path, format="PNG", dpi=(600, 600))
    white_background = Image.new("RGBA", transparent_image.size, "white")
    white_background.alpha_composite(transparent_image)
    white_background.convert("RGB").save(
        jpeg_path,
        format="JPEG",
        quality=95,
        dpi=(600, 600),
    )
    print(f"Wrote {pdf_path}")
    print(f"Wrote {png_path}")
    print(f"Wrote {jpeg_path}")


def render_reaction_sources(config: dict, output_dir: Path) -> None:
    """Render publication-resolution reaction schemes from their SMIRKS."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in ("ABF1", "SR", "ABF2"):
        smirks = config["catalog"]["reactions"][name]["smirks"]
        reaction = rdChemReactions.ReactionFromSmarts(smirks)
        if reaction is None:
            raise ValueError(f"Could not parse reaction SMIRKS for {name}: {smirks}")
        svg = Draw.ReactionToImage(
            reaction,
            subImgSize=(900, 450),
            useSVG=True,
        )
        svg_path = output_dir / f"{name}.svg"
        png_path = output_dir / f"{name}.png"
        svg_path.write_text(svg)
        subprocess.run(
            [
                "rsvg-convert",
                "--format=png",
                "--width=6000",
                f"--output={png_path}",
                str(svg_path),
            ],
            check=True,
        )


def build_figure5(output_dir: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="delt-hit-figure5-") as temporary_dir:
        temporary_root = Path(temporary_dir)
        graph_path = temporary_root / "reaction_graph.png"
        config = yaml.safe_load((CAMPAIGN_ROOT / "config.yaml").read_text())
        graph = prepare_graph_bundle(cfg=config)["G"]
        graph_axis = visualize_reaction_graph(
            graph,
            node_size=1400,
            reaction_node_size=1600,
            font_size=16,
        )
        graph_axis.figure.savefig(graph_path, dpi=300, bbox_inches="tight")
        close_figure(graph_axis.figure)

        reaction_root = temporary_root / "reactions"
        render_reaction_sources(config, reaction_root)

        source_path = temporary_root / "enumeration-summary.png"
        build_enumeration_composite(
            VISUALIZATION_ROOT,
            source_path,
            reaction_graph_path=graph_path,
            reaction_root=reaction_root,
        )
        save_original_composite(source_path, output_dir, "figure5-enumeration-summary")


def build_figure6(output_dir: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="delt-hit-figure6-") as temporary_dir:
        temporary_root = Path(temporary_dir)
        property_root = temporary_root / "properties"
        property_root.mkdir()
        render_molecular_weight_source(property_root / "prop_mw.png")
        source_path = temporary_root / "top-hit-properties-and-structures.png"
        build_properties_and_examples_composite(
            VISUALIZATION_ROOT,
            property_root,
            source_path,
        )
        save_original_composite(source_path, output_dir, "figure6-top-hit-properties-and-structures")


def parse_dataset_name(dataset_name: str) -> tuple[int, int, int] | None:
    match = DATASET_PATTERN.match(dataset_name)
    if match is None:
        return None
    return int(match["cycles"]), int(match["bbpc"]), int(match["depth"]) * 1_000_000


def format_reads(reads: int) -> str:
    return f"{reads // 1_000_000}M" if reads % 1_000_000 == 0 else str(reads)


def load_timings() -> dict[tuple[int, int], dict[str, list[tuple[int, float]]]]:
    grouped = defaultdict(lambda: defaultdict(list))
    for path in sorted(BENCHMARK_ROOT.glob("*/**/timing.json")):
        report = json.loads(path.read_text())
        parsed = parse_dataset_name(report["dataset_name"])
        if parsed:
            cycles, bbpc, reads = parsed
            grouped[(cycles, bbpc)][report["tool"]].append((reads, float(report["timings"]["total_s"])))
    return grouped


def load_memory() -> dict[tuple[int, int], dict[str, list[tuple[int, float]]]]:
    grouped = defaultdict(lambda: defaultdict(list))
    for path in sorted(BENCHMARK_ROOT.glob("*/**/job-stats.json")):
        report = json.loads(path.read_text())
        peak_rss = report.get("peak_rss_bytes")
        parts = path.relative_to(BENCHMARK_ROOT).parts
        parsed = parse_dataset_name(parts[1])
        if peak_rss is not None and parsed:
            cycles, bbpc, reads = parsed
            grouped[(cycles, bbpc)][parts[0]].append((reads, float(peak_rss) / 1024**3))
    return grouped


def plot_tool_lines(axis: plt.Axes, tool_points: dict[str, list[tuple[int, float]]]) -> None:
    for tool in ("deli", "delt"):
        points = sorted(tool_points.get(tool, []))
        if not points:
            continue
        style = TOOL_STYLES[tool]
        axis.plot(
            [point[0] for point in points],
            [point[1] for point in points],
            label=style["label"],
            color=style["color"],
            marker=style["marker"],
            markersize=3,
            linewidth=1,
        )
    all_reads = sorted({point[0] for points in tool_points.values() for point in points})
    axis.set_xscale("log")
    axis.set_xticks(all_reads, [format_reads(value) for value in all_reads])
    axis.set_xlabel("Number of reads")
    axis.grid(True, which="both", linestyle="--", linewidth=0.35, alpha=0.5)
    axis.legend(frameon=False)


def build_supplementary_figure1(output_dir: Path) -> None:
    timings = load_timings()
    fig, axes = plt.subplots(1, 3, figsize=(FULL_WIDTH_IN, 2.35), constrained_layout=True)
    for axis, cycles, panel in zip(axes, (2, 3, 4), ("a", "b", "c")):
        plot_tool_lines(axis, timings[(cycles, 10)])
        axis.set_ylabel("Runtime (s)")
        axis.set_title(f"{panel}  {cycles}-cycle library", loc="left", weight="bold")
    save_figure(fig, output_dir, "supplementary-figure1-runtime-cycles")


def build_supplementary_figure2(output_dir: Path) -> None:
    timings = load_timings()
    fig, axis = plt.subplots(figsize=(FULL_WIDTH_IN, 4.2), constrained_layout=True)
    line_styles = {10: "-", 100: "--", 1000: ":"}
    all_reads = set()
    for tool in ("deli", "delt"):
        for bbpc in (10, 100, 1000):
            points = sorted(timings[(2, bbpc)].get(tool, []))
            if not points:
                continue
            style = TOOL_STYLES[tool]
            axis.plot(
                [point[0] for point in points],
                [point[1] for point in points],
                label=f"{style['label']} {bbpc} BB/cycle",
                color=style["color"],
                marker=style["marker"],
                linestyle=line_styles[bbpc],
                markersize=3,
                linewidth=1,
            )
            all_reads.update(point[0] for point in points)
    axis.set_xscale("log")
    reads = sorted(all_reads)
    axis.set_xticks(reads, [format_reads(value) for value in reads])
    axis.set_xlabel("Number of reads")
    axis.set_ylabel("Runtime (s)")
    axis.grid(True, which="both", linestyle="--", linewidth=0.35, alpha=0.5)
    axis.legend(frameon=False, ncol=2)
    save_figure(fig, output_dir, "supplementary-figure2-runtime-building-blocks")


def build_supplementary_figure3(output_dir: Path) -> None:
    memory = load_memory()
    fig, axes = plt.subplots(1, 3, figsize=(FULL_WIDTH_IN, 2.35), constrained_layout=True)
    for axis, cycles, panel in zip(axes, (2, 3, 4), ("a", "b", "c")):
        plot_tool_lines(axis, memory[(cycles, 10)])
        axis.set_ylabel("Peak RSS (GiB)")
        axis.set_title(f"{panel}  {cycles}-cycle library", loc="left", weight="bold")
    save_figure(fig, output_dir, "supplementary-figure3-peak-memory")


def build_supplementary_figure4(output_dir: Path) -> None:
    hits = pd.read_csv(FAVALLI_ROOT / "enrichment" / "ca9_ds" / "hits_50.csv")
    frames = []
    for method in ("counts", "edgeR"):
        subset = hits[(hits["method"] == method) & (hits["replicate"] == "aggregate")]
        counts = subset["code_1"].value_counts().rename_axis("code_value").reset_index(name="count")
        counts["method"] = method
        counts["replicate"] = "aggregate"
        frames.append(counts)
    for replicate, subset in hits[(hits["method"] == "z_score") & (hits["replicate"] != "aggregate")].groupby("replicate"):
        counts = subset["code_1"].value_counts().rename_axis("code_value").reset_index(name="count")
        counts["method"] = "z_score"
        counts["replicate"] = replicate
        frames.append(counts)
    count_data = pd.concat(frames, ignore_index=True)
    top_codes = (
        count_data.groupby(["method", "code_value"], as_index=False)["count"]
        .mean()
        .groupby("code_value", as_index=False)["count"]
        .sum()
        .sort_values(["count", "code_value"], ascending=[False, True])
        .head(10)["code_value"]
    )
    count_data = count_data[count_data["code_value"].isin(top_codes)]
    fig, axis = plt.subplots(figsize=(FULL_WIDTH_IN, 4.0), constrained_layout=True)
    sns.barplot(
        data=count_data,
        x="code_value",
        y="count",
        hue="method",
        hue_order=METHOD_ORDER,
        palette=METHOD_PALETTE,
        estimator="mean",
        errorbar="sd",
        ax=axis,
    )
    axis.set_xlabel("Code 1 building block")
    axis.set_ylabel("Occurrences among top 50 compounds")
    axis.legend(title="Method", frameon=False)
    save_figure(fig, output_dir, "supplementary-figure4-favalli-top50-code1")


def main() -> None:
    args = parse_args()
    configure_style()
    build_figure2(args.output_dir)
    build_qc_component(
        args.output_dir,
        "2-B0",
        "figure3a-qc-building-block-0",
    )
    build_qc_component(
        args.output_dir,
        "1-C0",
        "figure3b-qc-constant-region-0",
    )
    build_qc_component(
        args.output_dir,
        "0-S0",
        "figure3c-qc-selection-region-0",
    )
    build_qc_composite(args.output_dir)
    build_box_b0_code0(args.output_dir)
    build_box3_reaction_template(args.output_dir)
    build_figure5(args.output_dir)
    build_figure6(args.output_dir)
    build_supplementary_figure1(args.output_dir)
    build_supplementary_figure2(args.output_dir)
    build_supplementary_figure3(args.output_dir)
    build_supplementary_figure4(args.output_dir)


if __name__ == "__main__":
    main()
