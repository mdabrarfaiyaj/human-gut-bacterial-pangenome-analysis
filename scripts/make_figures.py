#!/usr/bin/env python3
"""
make_figures.py

Generates the three figures required by the project:

  1. core_accessory_barplot.png - core/accessory/unique counts plus the
     gene-frequency spectrum that shows whether the pangenome is open.
  2. amr_heatmap.png            - AMR gene presence/absence, coloured by
     drug class, ordered by prevalence.
  3. phylogeny_amr.png          - core-genome tree with an aligned
     AMR presence grid and per-strain gene counts.

Usage:
    python3 make_figures.py \
        --gene-presence-absence results/panaroo/gene_presence_absence.csv \
        --amr-dir results/amr \
        --tree results/phylogeny/core_gene_alignment.aln.treefile \
        --outdir figures
"""
import argparse
import csv
import glob
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

# ---------------------------------------------------------------- styling ---

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.titleweight": "semibold",
    "axes.labelsize": 9,
    "axes.edgecolor": "#4a4a4a",
    "axes.linewidth": 0.8,
    "xtick.color": "#4a4a4a",
    "ytick.color": "#4a4a4a",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.bbox": "tight",
})

INK = "#2b2b2b"
MUTED = "#8a8a8a"
GRID = "#e6e6e6"

CATEGORY_COLOURS = {
    "Core": "#2f5d8a",
    "Accessory": "#5fa8a0",
    "Unique": "#d8b95f",
}

# One colour per drug class, muted but distinguishable.
CLASS_COLOURS = {
    "TETRACYCLINE": "#3d6f9e",
    "MACROLIDE": "#c2703d",
    "LINCOSAMIDE": "#7a9b5e",
    "LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN": "#8a6fa8",
    "BETA-LACTAM": "#b4574f",
}
CLASS_FALLBACK = "#6f6f6f"

PANAROO_META = {
    "Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
    "Avg sequences per isolate", "Genome Fragment", "Order within Fragment",
    "Accessory Fragment", "Accessory Order with Fragment", "QC",
    "Min group size nuc", "Max group size nuc", "Avg group size nuc",
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--gene-presence-absence", required=True)
    p.add_argument("--amr-dir", required=True)
    p.add_argument("--tree", required=False)
    p.add_argument("--outdir", required=True)
    p.add_argument("--core-threshold", type=float, default=0.98)
    return p.parse_args()


def short(name):
    """Bvulgatus_003 -> 003, to keep tick labels compact."""
    return name.split("_")[-1]


# --------------------------------------------------------------- loading ---

def load_presence(gpa_csv):
    df = pd.read_csv(gpa_csv, low_memory=False)
    sample_cols = [c for c in df.columns if c not in PANAROO_META]
    presence = df[sample_cols].notna() & (df[sample_cols] != "")
    return presence, sample_cols


def load_amr(amr_dir):
    """Return (matrix DataFrame genes x samples, dict gene -> drug class)."""
    files = sorted(glob.glob(os.path.join(amr_dir, "*.amr.tsv")))
    samples = [os.path.basename(f).replace(".amr.tsv", "") for f in files]
    hits = defaultdict(set)
    gene_class = {}

    for path, sample in zip(files, samples):
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                gene = row.get("Gene symbol") or row.get("gene_symbol")
                if not gene:
                    continue
                hits[gene].add(sample)
                drug = row.get("Class") or row.get("class") or ""
                gene_class[gene] = drug.strip().upper()

    genes = sorted(hits, key=lambda g: (-len(hits[g]), g))
    mat = pd.DataFrame(0, index=genes, columns=samples, dtype=int)
    for gene, carriers in hits.items():
        for s in carriers:
            mat.loc[gene, s] = 1
    return mat, gene_class


# -------------------------------------------------------------- figure 1 ---

def figure1(presence, sample_cols, outdir, core_threshold):
    n = len(sample_cols)
    counts = presence.sum(axis=1)

    core = int((counts >= core_threshold * n).sum())
    unique = int((counts == 1).sum())
    accessory = int(len(counts) - core - unique)
    total = core + accessory + unique

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 4.1),
                                   gridspec_kw={"width_ratios": [1, 1.35]})

    # (a) category counts
    cats = ["Core", "Accessory", "Unique"]
    vals = [core, accessory, unique]
    bars = ax1.bar(cats, vals, color=[CATEGORY_COLOURS[c] for c in cats],
                   width=0.62, edgecolor="white", linewidth=0.8)
    for bar, v in zip(bars, vals):
        ax1.text(bar.get_x() + bar.get_width() / 2, v + total * 0.012,
                 f"{v:,}", ha="center", va="bottom", fontsize=9.5,
                 color=INK, fontweight="semibold")
        ax1.text(bar.get_x() + bar.get_width() / 2, v / 2,
                 f"{100 * v / total:.0f}%", ha="center", va="center",
                 fontsize=9, color="white", fontweight="semibold")

    ax1.set_ylabel("Gene clusters")
    ax1.set_title(f"a  Pangenome composition (n = {n})", loc="left")
    ax1.set_ylim(0, max(vals) * 1.16)
    ax1.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax1.set_axisbelow(True)
    for side in ("top", "right"):
        ax1.spines[side].set_visible(False)

    # (b) gene frequency spectrum
    spectrum = [int((counts == k).sum()) for k in range(1, n + 1)]
    colours = []
    for k in range(1, n + 1):
        if k >= core_threshold * n:
            colours.append(CATEGORY_COLOURS["Core"])
        elif k == 1:
            colours.append(CATEGORY_COLOURS["Unique"])
        else:
            colours.append(CATEGORY_COLOURS["Accessory"])

    ax2.bar(range(1, n + 1), spectrum, color=colours, width=0.72,
            edgecolor="white", linewidth=0.6)
    ax2.set_xlabel("Number of strains sharing a gene")
    ax2.set_ylabel("Gene clusters")
    ax2.set_title("b  Gene frequency spectrum", loc="left")
    ax2.set_xticks(range(1, n + 1))
    ax2.set_xticklabels([str(k) for k in range(1, n + 1)], fontsize=7.5)
    ax2.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax2.set_axisbelow(True)
    for side in ("top", "right"):
        ax2.spines[side].set_visible(False)

    handles = [Rectangle((0, 0), 1, 1, color=CATEGORY_COLOURS[c]) for c in cats]
    ax2.legend(handles, cats, frameon=False, fontsize=8,
               loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.02))

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "core_accessory_barplot.png"), dpi=300)
    plt.close(fig)

    return {"core": core, "accessory": accessory, "unique": unique,
            "n_samples": n, "total": total}


# -------------------------------------------------------------- figure 2 ---

def figure2(mat, gene_class, outdir):
    if mat.empty:
        fig, ax = plt.subplots(figsize=(7, 1.8))
        ax.text(0.5, 0.5, "No AMR genes detected", ha="center", va="center",
                fontsize=11, color=INK)
        ax.axis("off")
        fig.savefig(os.path.join(outdir, "amr_heatmap.png"), dpi=300)
        plt.close(fig)
        return

    genes = list(mat.index)
    samples = list(mat.columns)
    n_g, n_s = len(genes), len(samples)

    fig, ax = plt.subplots(figsize=(max(7.5, 0.46 * n_s + 3.6),
                                    max(3.2, 0.36 * n_g + 1.9)))

    for i, gene in enumerate(genes):
        colour = CLASS_COLOURS.get(gene_class.get(gene, ""), CLASS_FALLBACK)
        for j, sample in enumerate(samples):
            present = mat.loc[gene, sample] == 1
            ax.add_patch(Rectangle(
                (j, i), 1, 1,
                facecolor=colour if present else "#f5f5f5",
                edgecolor="white", linewidth=1.4))

    ax.set_xlim(0, n_s)
    ax.set_ylim(0, n_g)
    ax.invert_yaxis()
    ax.set_aspect("equal")

    ax.set_xticks(np.arange(n_s) + 0.5)
    ax.set_xticklabels([short(s) for s in samples], fontsize=8)
    ax.set_yticks(np.arange(n_g) + 0.5)
    ax.set_yticklabels([f"{g}  ({int(mat.loc[g].sum())}/{n_s})" for g in genes],
                       fontsize=8.5)
    ax.set_xlabel("Strain (Bvulgatus_ prefix omitted)", labelpad=8)
    ax.set_title("AMR gene presence by strain (AMRFinderPlus)", loc="left", pad=10)

    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)

    present_classes = []
    for g in genes:
        c = gene_class.get(g, "")
        if c and c not in present_classes:
            present_classes.append(c)

    def pretty(c):
        return "MLS" if c == "LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN" else c.title()

    handles = [Rectangle((0, 0), 1, 1,
                         facecolor=CLASS_COLOURS.get(c, CLASS_FALLBACK))
               for c in present_classes]
    ax.legend(handles, [pretty(c) for c in present_classes],
              frameon=False, fontsize=8, title="Drug class",
              title_fontsize=8.5, loc="upper left",
              bbox_to_anchor=(1.01, 1.0), handlelength=1.1)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "amr_heatmap.png"), dpi=300)
    plt.close(fig)


# -------------------------------------------------------------- figure 3 ---

def tree_layout(tree):
    """Return tip order, and dicts of x (root distance) and y per clade."""
    tree.ladderize()
    tips = tree.get_terminals()
    ys = {tip: float(i) for i, tip in enumerate(tips)}

    def assign_y(clade):
        if clade.is_terminal():
            return ys[clade]
        vals = [assign_y(c) for c in clade.clades]
        y = (min(vals) + max(vals)) / 2.0
        ys[clade] = y
        return y

    assign_y(tree.root)

    xs = {}

    def assign_x(clade, acc):
        length = clade.branch_length or 0.0
        xs[clade] = acc + length
        for c in clade.clades:
            assign_x(c, xs[clade])

    assign_x(tree.root, 0.0)
    return tips, xs, ys


def figure3(tree_path, mat, gene_class, outdir):
    from Bio import Phylo

    tree = Phylo.read(tree_path, "newick")
    tips, xs, ys = tree_layout(tree)
    tip_names = [t.name for t in tips]
    n = len(tip_names)

    # Compress any grossly long branch so the rest of the tree stays legible.
    lengths = [c.branch_length or 0.0 for c in tree.find_clades()
               if c.branch_length]
    cap = float(np.percentile(lengths, 90)) * 4 if lengths else 0.0
    compressed = {}
    if cap > 0:
        for clade in tree.find_clades():
            bl = clade.branch_length or 0.0
            if bl > cap:
                compressed[clade] = bl
        if compressed:
            for clade in compressed:
                clade.branch_length = cap
            xs = {}

            def assign_x(clade, acc):
                xs[clade] = acc + (clade.branch_length or 0.0)
                for c in clade.clades:
                    assign_x(c, xs[clade])
            assign_x(tree.root, 0.0)

    amr_genes = list(mat.index) if not mat.empty else []
    grid_w = max(3.0, 0.30 * len(amr_genes) + 1.2)
    fig, (ax_t, ax_g) = plt.subplots(
        1, 2, figsize=(6.6 + grid_w, max(4.2, 0.36 * n + 1.8)),
        gridspec_kw={"width_ratios": [3.1, grid_w / 2.3]})

    # --- tree panel
    for clade in tree.find_clades(order="level"):
        x_end = xs[clade]
        x_start = x_end - (clade.branch_length or 0.0)
        y = ys[clade]
        ax_t.plot([x_start, x_end], [y, y], color=INK, linewidth=1.1,
                  solid_capstyle="round", zorder=2)
        if clade in compressed:
            xm = (x_start + x_end) / 2
            for off in (-0.012 * max(xs.values()), 0.012 * max(xs.values())):
                ax_t.plot([xm + off - 0.004 * max(xs.values()), xm + off + 0.004 * max(xs.values())],
                          [y - 0.22, y + 0.22], color="white", linewidth=2.4, zorder=3)
                ax_t.plot([xm + off - 0.004 * max(xs.values()), xm + off + 0.004 * max(xs.values())],
                          [y - 0.22, y + 0.22], color=MUTED, linewidth=0.9, zorder=4)
        if not clade.is_terminal():
            child_ys = [ys[c] for c in clade.clades]
            ax_t.plot([x_end, x_end], [min(child_ys), max(child_ys)],
                      color=INK, linewidth=1.1, zorder=2)

    x_max = max(xs.values())

    # bootstrap labels: only where support is not maximal, so the tree stays clean
    for clade in tree.get_nonterminals():
        conf = None
        if clade.confidence is not None:
            conf = clade.confidence
        elif clade.name:
            try:
                conf = float(clade.name)
            except ValueError:
                conf = None
        if conf is not None and conf < 100:
            ax_t.text(xs[clade] - x_max * 0.006, ys[clade] - 0.24, f"{int(conf)}",
                      fontsize=7, color=MUTED, ha="right", va="center")

    for tip in tips:
        ax_t.plot([xs[tip], x_max * 1.035], [ys[tip], ys[tip]],
                  color=GRID, linewidth=0.7, zorder=1)
        ax_t.text(x_max * 1.055, ys[tip], tip.name, fontsize=8.5,
                  va="center", ha="left", color=INK)

    ax_t.set_ylim(n - 0.4, -0.6)
    ax_t.set_xlim(-x_max * 0.02, x_max * 1.42)
    ax_t.set_yticks([])
    ax_t.set_xlabel("Substitutions per site"
                    + ("   ( //  long branch compressed )" if compressed else ""))
    ax_t.set_title("a  Core-genome ML phylogeny (IQ-TREE, 1000 bootstraps)",
                   loc="left", pad=10)
    for side in ("top", "right", "left"):
        ax_t.spines[side].set_visible(False)

    # --- AMR grid panel, rows aligned to tree tips
    if amr_genes:
        for j, gene in enumerate(amr_genes):
            colour = CLASS_COLOURS.get(gene_class.get(gene, ""), CLASS_FALLBACK)
            for i, name in enumerate(tip_names):
                present = name in mat.columns and mat.loc[gene, name] == 1
                ax_g.add_patch(Rectangle(
                    (j, i - 0.5), 1, 1,
                    facecolor=colour if present else "#f5f5f5",
                    edgecolor="white", linewidth=1.2))

        counts = [int(mat[name].sum()) if name in mat.columns else 0
                  for name in tip_names]
        for i, c in enumerate(counts):
            ax_g.text(len(amr_genes) + 0.55, i, str(c), fontsize=8.5,
                      va="center", ha="center", color=INK)
        ax_g.text(len(amr_genes) + 0.55, -1.05, "n", fontsize=8.5,
                  va="center", ha="center", color=MUTED, fontweight="semibold")

        ax_g.set_xlim(0, len(amr_genes) + 1.2)
        ax_g.set_ylim(n - 0.4, -0.6)
        ax_g.set_xticks(np.arange(len(amr_genes)) + 0.5)
        ax_g.set_xticklabels(amr_genes, rotation=90, fontsize=8)
        ax_g.set_yticks([])
        ax_g.set_title("b  AMR genes", loc="left", pad=10)
        for side in ("top", "right", "left", "bottom"):
            ax_g.spines[side].set_visible(False)
        ax_g.tick_params(length=0)
    else:
        ax_g.text(0.5, 0.5, "No AMR genes detected", ha="center", va="center",
                  fontsize=10, color=INK)
        ax_g.axis("off")

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "phylogeny_amr.png"), dpi=300)
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    presence, sample_cols = load_presence(args.gene_presence_absence)
    stats = figure1(presence, sample_cols, args.outdir, args.core_threshold)
    print("Figure 1 written:", stats)

    mat, gene_class = load_amr(args.amr_dir)
    figure2(mat, gene_class, args.outdir)
    print(f"Figure 2 written: {mat.shape[0]} AMR genes x {mat.shape[1]} strains")

    if not args.tree or not os.path.exists(args.tree):
        raise FileNotFoundError(f"Tree file not found: {args.tree}")
    figure3(args.tree, mat, gene_class, args.outdir)
    print("Figure 3 written")

    with open(os.path.join(args.outdir, "figure_stats.txt"), "w") as f:
        f.write(f"core_genes={stats['core']}\n")
        f.write(f"accessory_genes={stats['accessory']}\n")
        f.write(f"unique_genes={stats['unique']}\n")
        f.write(f"n_samples={stats['n_samples']}\n")
        f.write(f"total_genes={stats['total']}\n")

    print("Done ->", args.outdir)


if __name__ == "__main__":
    main()
