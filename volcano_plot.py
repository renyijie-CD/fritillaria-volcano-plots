#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Manuscript-style volcano plot script for processed volcano-plot data tables.

Input table should contain:
    cluster_id / Name / mz / rt / log2FC / FC / p_value

This script reproduces the manuscript volcano-plot style:
    - non-significant features: grey
    - significant features: blue-to-lightgrey-to-orange-to-deepred color scale by log2FC
    - dashed threshold lines at log2FC = +/-1 and p = 0.05
    - optional Down/Up arrows
    - optional axis labels and colorbar label
    - each comparison can use its own x/y limits
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import FancyArrowPatch


LOG2FC_ALIASES = [
    "log2FC", "log2fc", "log2_fc", "Log2FC", "LOG2FC",
    "log2(FoldChange)", "log2FoldChange", "logFC", "LogFC"
]

PVALUE_ALIASES = [
    "pvalue", "p_value", "Pvalue", "PValue", "P.Value", "p.value",
    "P", "p", "pval", "Pval", "PVAL"
]

FEATURE_ALIASES = [
    "cluster_id", "Cluster_ID", "Metabolite", "metabolite", "Name", "name",
    "Compound", "compound", "Feature", "feature", "ID", "id"
]


def read_table(input_file, sheet_name=0):
    input_file = Path(input_file)
    suffix = input_file.suffix.lower()
    if suffix in [".xlsx", ".xls"]:
        return pd.read_excel(input_file, sheet_name=sheet_name)
    if suffix == ".csv":
        return pd.read_csv(input_file)
    if suffix in [".tsv", ".txt"]:
        return pd.read_csv(input_file, sep="\t")
    raise ValueError(f"Unsupported file type: {suffix}. Use .xlsx, .xls, .csv, .tsv, or .txt")


def find_first_existing(columns, aliases):
    for alias in aliases:
        if alias in columns:
            return alias

    cleaned = {
        str(c).replace(" ", "").replace("_", "").replace("-", "").replace(".", "").lower(): c
        for c in columns
    }
    for alias in aliases:
        key = alias.replace(" ", "").replace("_", "").replace("-", "").replace(".", "").lower()
        if key in cleaned:
            return cleaned[key]
    return None


def calculate_from_stats_table(df, log2fc_col=None, pvalue_col=None, feature_col=None):
    columns = list(df.columns)

    log2fc_col = log2fc_col or find_first_existing(columns, LOG2FC_ALIASES)
    pvalue_col = pvalue_col or find_first_existing(columns, PVALUE_ALIASES)
    feature_col = feature_col or find_first_existing(columns, FEATURE_ALIASES)

    if log2fc_col is None or pvalue_col is None:
        raise ValueError("Cannot find log2FC/pvalue columns. Specify --log2fc-col and --pvalue-col.")

    out = pd.DataFrame()
    if feature_col and feature_col in df.columns:
        out["Feature"] = df[feature_col].astype(str)
    else:
        out["Feature"] = [f"Feature_{i+1}" for i in range(len(df))]

    out["log2FC"] = pd.to_numeric(df[log2fc_col], errors="coerce")
    out["pvalue"] = pd.to_numeric(df[pvalue_col], errors="coerce")
    return out


def add_regulation(df, fc_cutoff=1.0, p_cutoff=0.05):
    out = df.copy()
    out = out.dropna(subset=["log2FC", "pvalue"]).copy()
    out = out[out["pvalue"] > 0].copy()

    out["neg_log10_p"] = -np.log10(out["pvalue"])

    out["Regulation"] = "NS"
    out.loc[(out["log2FC"] >= fc_cutoff) & (out["pvalue"] < p_cutoff), "Regulation"] = "Up"
    out.loc[(out["log2FC"] <= -fc_cutoff) & (out["pvalue"] < p_cutoff), "Regulation"] = "Down"
    return out


def parse_range(text):
    parts = [float(x.strip()) for x in str(text).split(",")]
    if len(parts) != 2:
        raise ValueError("Range should be like -4.0,6.6")
    return parts[0], parts[1]


def plot_volcano(
    volcano_df,
    output_prefix,
    fc_cutoff=1.0,
    p_cutoff=0.05,
    xlim=(-4.0, 6.6),
    ylim=(0, 14.0),
    color_vmin=-3.5,
    color_vmax=4.0,
    point_size=42,
    sig_alpha=0.96,
    ns_alpha=0.90,
    dpi=600,
    show_labels=True,
    show_arrows=True,
    show_colorbar_label=True,
):
    plot_df = volcano_df.dropna(subset=["log2FC", "pvalue", "neg_log10_p"]).copy()

    ns_df = plot_df[plot_df["Regulation"] == "NS"]
    sig_df = plot_df[plot_df["Regulation"] != "NS"]

    up_n = int((plot_df["Regulation"] == "Up").sum())
    down_n = int((plot_df["Regulation"] == "Down").sum())

    # Calibrated against the user's original volcano plot:
    # blue -> light grey -> orange-red -> deep red
    cmap = LinearSegmentedColormap.from_list(
        "blue_grey_orangered_deepred",
        [
            (0.00, "#3F5FD0"),
            (0.48, "#D8D8D8"),
            (0.70, "#F28A63"),
            (1.00, "#C0002B"),
        ]
    )
    norm = TwoSlopeNorm(vmin=color_vmin, vcenter=0, vmax=color_vmax)

    fig, ax = plt.subplots(figsize=(4.8, 4.8), dpi=dpi)

    # Non-significant features
    ax.scatter(
        ns_df["log2FC"],
        ns_df["neg_log10_p"],
        c="#D0D0D0",
        s=point_size,
        alpha=ns_alpha,
        edgecolors="none",
        zorder=1,
    )

    # Significant features
    sc = ax.scatter(
        sig_df["log2FC"],
        sig_df["neg_log10_p"],
        c=sig_df["log2FC"],
        cmap=cmap,
        norm=norm,
        s=point_size,
        alpha=sig_alpha,
        edgecolors="none",
        zorder=2,
    )

    # Threshold lines
    ax.axvline(-fc_cutoff, linestyle="--", color="#808080", linewidth=0.8, zorder=0)
    ax.axvline(fc_cutoff, linestyle="--", color="#808080", linewidth=0.8, zorder=0)
    ax.axhline(-np.log10(p_cutoff), linestyle="--", color="#808080", linewidth=0.8, zorder=0)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks(np.arange(np.ceil(xlim[0]), np.floor(xlim[1]) + 1, 2))
    ax.set_yticks(np.arange(0, np.floor(ylim[1]) + 1, 2))

    if show_labels:
        ax.set_xlabel(r"log$_2$FC", fontsize=16, fontname="Times New Roman")
        ax.set_ylabel(r"-log$_{10}$Pvalue", fontsize=16, fontname="Times New Roman")
    else:
        ax.set_xlabel("")
        ax.set_ylabel("")

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontname("Times New Roman")
        label.set_fontsize(12)

    # Manuscript axis style
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)
    ax.tick_params(direction="in", length=3.5, width=0.8, colors="black")

    if show_arrows:
        arrow_y = ylim[0] + 0.55
        text_y = arrow_y + 0.25

        down_arrow = FancyArrowPatch(
            posA=(-fc_cutoff - 0.05, arrow_y),
            posB=(xlim[0] + 0.12, arrow_y),
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=1.2,
            color="#3F5FD0",
            zorder=3,
        )
        ax.add_patch(down_arrow)
        ax.text(
            (xlim[0] - fc_cutoff) / 2,
            text_y,
            f"Down {down_n}",
            ha="center",
            va="bottom",
            fontsize=11,
            fontname="Times New Roman",
            color="black",
        )

        up_arrow = FancyArrowPatch(
            posA=(fc_cutoff + 0.05, arrow_y),
            posB=(xlim[1] - 0.12, arrow_y),
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=1.2,
            color="#C0002B",
            zorder=3,
        )
        ax.add_patch(up_arrow)
        ax.text(
            (xlim[1] + fc_cutoff) / 2,
            text_y,
            f"Up {up_n}",
            ha="center",
            va="bottom",
            fontsize=11,
            fontname="Times New Roman",
            color="black",
        )

    cbar = plt.colorbar(sc, ax=ax, fraction=0.052, pad=0.045)
    if show_colorbar_label:
        cbar.set_label(r"log$_2$FC", fontsize=16, fontname="Times New Roman", rotation=270, labelpad=18)
    else:
        cbar.set_label("")
    cbar.ax.tick_params(labelsize=10, length=3.0, width=0.8, direction="in")
    for label in cbar.ax.get_yticklabels():
        label.set_fontname("Times New Roman")
    cbar.outline.set_linewidth(0.8)

    plt.tight_layout()

    output_prefix = Path(output_prefix)
    tif_path = output_prefix.with_suffix(".tif")
    pdf_path = output_prefix.with_suffix(".pdf")
    png_path = output_prefix.with_suffix(".png")

    fig.savefig(tif_path, dpi=dpi, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return {
        "up_n": up_n,
        "down_n": down_n,
        "tif": tif_path,
        "pdf": pdf_path,
        "png": png_path,
    }


def main():
    parser = argparse.ArgumentParser(description="Draw manuscript-style volcano plot from processed volcano-plot data.")
    parser.add_argument("-i", "--input", required=True, help="Input file: .xlsx, .xls, .csv, .tsv, .txt")
    parser.add_argument("-o", "--output-prefix", default="Volcano_plot", help="Output prefix without extension")
    parser.add_argument("--sheet", default="0", help="Excel sheet name or index. Default: 0")
    parser.add_argument("--feature-col", default=None)
    parser.add_argument("--log2fc-col", default=None)
    parser.add_argument("--pvalue-col", default=None)
    parser.add_argument("--fc-cutoff", type=float, default=1.0)
    parser.add_argument("--p-cutoff", type=float, default=0.05)
    parser.add_argument("--xlim", default="-4.0,6.6")
    parser.add_argument("--ylim", default="0,14.0")
    parser.add_argument("--color-vmin", type=float, default=-3.5)
    parser.add_argument("--color-vmax", type=float, default=4.0)
    parser.add_argument("--point-size", type=float, default=42)
    parser.add_argument("--no-labels", action="store_true")
    parser.add_argument("--no-arrows", action="store_true")
    parser.add_argument("--no-colorbar-label", action="store_true")

    args = parser.parse_args()

    sheet = args.sheet
    try:
        sheet = int(sheet)
    except Exception:
        pass

    df = read_table(args.input, sheet_name=sheet)
    volcano = calculate_from_stats_table(df, args.log2fc_col, args.pvalue_col, args.feature_col)
    volcano = add_regulation(volcano, fc_cutoff=args.fc_cutoff, p_cutoff=args.p_cutoff)

    output_prefix = Path(args.output_prefix)
    result_csv = output_prefix.with_name(output_prefix.name + "_volcano_result.csv")
    volcano.to_csv(result_csv, index=False, encoding="utf-8-sig")

    info = plot_volcano(
        volcano,
        output_prefix=output_prefix,
        fc_cutoff=args.fc_cutoff,
        p_cutoff=args.p_cutoff,
        xlim=parse_range(args.xlim),
        ylim=parse_range(args.ylim),
        color_vmin=args.color_vmin,
        color_vmax=args.color_vmax,
        point_size=args.point_size,
        show_labels=not args.no_labels,
        show_arrows=not args.no_arrows,
        show_colorbar_label=not args.no_colorbar_label,
    )

    print("=" * 60)
    print("Volcano plot finished.")
    print(f"Actual Up: {info['up_n']}")
    print(f"Actual Down: {info['down_n']}")
    print(f"Result table: {result_csv}")
    print(f"TIF: {info['tif']}")
    print(f"PDF: {info['pdf']}")
    print(f"PNG: {info['png']}")
    print("=" * 60)


if __name__ == "__main__":
    main()
