#!/usr/bin/env python3
"""
QC plotter for MiXCR align_stats.csv / assemble_stats.csv (output of
log_parser.py). Merges the two on sample name, derives QC percentages
(the parser strips "(xx.x%)" annotations down to raw counts, so
percentages are recomputed here from Total sequencing reads), and
produces a set of cohort-level QC plots plus a merged metrics table.

USAGE:
    python3 qc_plotter.py --align_stats align_stats.csv \
                           --assemble_stats assemble_stats.csv \
                           --out_dir qc_plots/
"""

import argparse
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

CHAINS = ["TRA", "TRB", "TRG", "TRD"]

# === Plot style: clean, minimal, shared conventions across figures ===
plt.rcParams.update({
    "figure.dpi": 150,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": "#333333",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
})
BAR_COLOR = "#4C72B0"
THRESHOLD_COLOR = "#C44E52"
CHAIN_COLORS = {"TRA": "#4C72B0", "TRB": "#55A868", "TRG": "#C44E52", "TRD": "#8172B2"}


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("--align_stats", required=True, help="path to align_stats.csv")
    parser.add_argument("--assemble_stats", required=True, help="path to assemble_stats.csv")
    parser.add_argument("--out_dir", required=True, help="directory to write QC plots + merged table")
    parser.add_argument("--align_rate_threshold", type=float, default=80.0,
                         help="dashed threshold line for alignment rate plot (percent)")
    parser.add_argument("--reads_used_threshold", type=float, default=50.0,
                         help="dashed threshold line for reads-used-in-clonotypes plot (percent)")
    return parser.parse_args()


def sample_from_align_input(val):
    """'SAMPLE_barcode_R1.fastq,SAMPLE_barcode_R2.fastq' -> 'SAMPLE'"""
    first = str(val).split(",")[0].strip()
    return re.sub(r"_barcode_R1\.fastq$", "", first)


def sample_from_assemble_input(val):
    """'extended_SAMPLE.vdjca' -> 'SAMPLE'"""
    base = str(val).strip()
    base = re.sub(r"^extended_", "", base)
    base = re.sub(r"\.vdjca$", "", base)
    return base


def load_data(align_path, assemble_path):
    align = pd.read_csv(align_path)
    assemble = pd.read_csv(assemble_path)

    align.columns = [c.strip() for c in align.columns]
    assemble.columns = [c.strip() for c in assemble.columns]

    align["Sample"] = align["Input file(s)"].apply(sample_from_align_input)
    assemble["Sample"] = assemble["Input file(s)"].apply(sample_from_assemble_input)

    merged = pd.merge(align, assemble, on="Sample", suffixes=("_align", "_assemble"))

    def col(df, name, default=np.nan):
        """Fetch a column if present, else a NaN-filled Series of the right length."""
        if name in df.columns:
            return pd.to_numeric(df[name], errors="coerce")
        print(f"WARNING: column '{name}' not found in merged stats — filling with NaN")
        return pd.Series(default, index=df.index)

    # Derived QC percentages (recomputed since the parser dropped the % text)
    merged["Alignment rate (%)"] = (
        col(merged, "Successfully aligned reads") / col(merged, "Total sequencing reads") * 100
    )
    merged["Reads used in clonotypes (%)"] = (
        col(merged, "Reads used in clonotypes, percent of total") / col(merged, "Total sequencing reads") * 100
    )

    # Chain fractions from the final assembled clonotypes
    chain_cols = [f"{c} chains_assemble" for c in CHAINS]
    for cc in chain_cols:
        if cc not in merged.columns:
            merged[cc] = 0
    merged["Total TCR chain reads"] = merged[chain_cols].apply(pd.to_numeric, errors="coerce").sum(axis=1)
    for c in CHAINS:
        merged[f"{c} fraction (%)"] = (
            pd.to_numeric(merged[f"{c} chains_assemble"], errors="coerce")
            / merged["Total TCR chain reads"].replace(0, np.nan) * 100
        )

    merged = merged.sort_values("Sample").reset_index(drop=True)
    return merged


def _bar(ax, samples, values, title, ylabel, threshold=None, color=BAR_COLOR):
    x = np.arange(len(samples))
    ax.bar(x, values, color=color, width=0.7)
    if threshold is not None:
        ax.axhline(threshold, color=THRESHOLD_COLOR, linestyle="--", linewidth=1,
                    label=f"threshold = {threshold:g}")
        ax.legend(frameon=False, loc="upper right")
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=90, fontsize=6)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(-0.6, len(samples) - 0.4)


def plot_sequencing_depth(df, out_dir):
    fig, ax = plt.subplots(figsize=(max(8, len(df) * 0.25), 4))
    _bar(ax, df["Sample"], df["Total sequencing reads"],
         "Sequencing depth per sample", "Total sequencing reads")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "01_sequencing_depth.png"))
    plt.close(fig)


def plot_alignment_rate(df, out_dir, threshold):
    fig, ax = plt.subplots(figsize=(max(8, len(df) * 0.25), 4))
    _bar(ax, df["Sample"], df["Alignment rate (%)"],
         "Alignment rate per sample", "Successfully aligned reads (%)", threshold=threshold)
    ax.set_ylim(0, 100)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "02_alignment_rate.png"))
    plt.close(fig)


def plot_reads_used_in_clonotypes(df, out_dir, threshold):
    fig, ax = plt.subplots(figsize=(max(8, len(df) * 0.25), 4))
    _bar(ax, df["Sample"], df["Reads used in clonotypes (%)"],
         "Reads used in final clonotypes per sample", "Reads used in clonotypes (%)", threshold=threshold)
    ax.set_ylim(0, 100)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "03_reads_used_in_clonotypes.png"))
    plt.close(fig)


def plot_clonotype_count(df, out_dir):
    fig, ax = plt.subplots(figsize=(max(8, len(df) * 0.25), 4))
    _bar(ax, df["Sample"], df["Final clonotype count"],
         "Final clonotype count per sample", "Clonotype count", color="#55A868")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "04_clonotype_count.png"))
    plt.close(fig)


def plot_chain_distribution(df, out_dir):
    fig, ax = plt.subplots(figsize=(max(8, len(df) * 0.25), 4))
    x = np.arange(len(df))
    bottom = np.zeros(len(df))
    for c in CHAINS:
        vals = df[f"{c} fraction (%)"].fillna(0).values
        ax.bar(x, vals, bottom=bottom, color=CHAIN_COLORS[c], width=0.7, label=c)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(df["Sample"], rotation=90, fontsize=6)
    ax.set_ylabel("Chain fraction of clonotype reads (%)")
    ax.set_title("TCR chain distribution per sample")
    ax.set_ylim(0, 100)
    ax.set_xlim(-0.6, len(df) - 0.4)
    ax.legend(frameon=False, ncol=len(CHAINS), loc="upper center", bbox_to_anchor=(0.5, -0.35))
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "05_chain_distribution.png"))
    plt.close(fig)


def plot_clone_size_vs_count(df, out_dir):
    if "Average number of reads per clonotype" not in df.columns:
        print("WARNING: 'Average number of reads per clonotype' not found — skipping clone-size plot")
        return
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(df["Final clonotype count"], df["Average number of reads per clonotype"],
               color=BAR_COLOR, edgecolor="white", s=40)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Final clonotype count")
    ax.set_ylabel("Average reads per clonotype")
    ax.set_title("Clonality vs. clone-level sequencing depth")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "06_clone_size_vs_count.png"))
    plt.close(fig)


def main():
    args = get_options()
    os.makedirs(args.out_dir, exist_ok=True)

    df = load_data(args.align_stats, args.assemble_stats)

    plot_sequencing_depth(df, args.out_dir)
    plot_alignment_rate(df, args.out_dir, args.align_rate_threshold)
    plot_reads_used_in_clonotypes(df, args.out_dir, args.reads_used_threshold)
    plot_clonotype_count(df, args.out_dir)
    plot_chain_distribution(df, args.out_dir)
    plot_clone_size_vs_count(df, args.out_dir)

    metrics_cols = ["Sample", "Total sequencing reads", "Successfully aligned reads",
                     "Alignment rate (%)", "Final clonotype count",
                     "Reads used in clonotypes (%)", "Average number of reads per clonotype"] + \
                    [f"{c} fraction (%)" for c in CHAINS]
    df.reindex(columns=metrics_cols).to_csv(
        os.path.join(args.out_dir, "qc_summary_metrics.tsv"), sep="\t", index=False)

    print(f"Wrote 6 QC plots and qc_summary_metrics.tsv to {args.out_dir}")


if __name__ == "__main__":
    main()
