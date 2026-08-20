#!/usr/bin/env python3
"""
Identify proteins significantly altered by Alprazolam (Alp+CD vs CD) and
REVERSED by the knockdown (Alp+KD vs Alp+CD), then visualise the top 50 as a
z-scored abundance heatmap across the four conditions, annotated with each
protein's Reactome pathway.

Top 50 = ranked by combined reversal strength (see RANK_BY below).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns


# ---------------------------------------------------------------------
# Global font settings  <<< CHANGE FONT SIZES HERE >>>
# (same block as Volcano_plot.py / pca_samples_advanced.py — keep in sync)
# ---------------------------------------------------------------------

plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size":        10,   # baseline (fallback for anything not set below)
    "axes.labelsize":   14,   # x/y axis labels
    "axes.titlesize":   14,   # per-panel titles
    "figure.titlesize": 14,   # fig.suptitle
    "xtick.labelsize":  12,   # numbers on the x axis
    "ytick.labelsize":  12,   # numbers on the y axis
    "legend.fontsize":  12,   # legend text
    "legend.title_fontsize": 12,
    # Keep text as real text (not outlines) in SVG/PDF so it stays editable
    # in Inkscape / Word / PowerPoint. Set to "path" if fonts look wrong on
    # another computer.
    "svg.fonttype":     "none",
    "pdf.fonttype":     42,   # TrueType in the PDF (editable text, not bitmaps)
})

# ---------------------------------------------------------------------
# Heatmap / row-label font sizes
# ---------------------------------------------------------------------
# These figures pack TOP_N rows of small cells, so their text does NOT follow
# the global sizes above — at 12 pt the padj2 annotations and gene symbols
# would overlap. Tune them separately here.

HEATMAP_ANNOT_FONTSIZE   = 8      # padj2 numbers inside the log2FC cells
HEATMAP_XTICK_FONTSIZE   = 10     # comparison labels under the log2FC panel
ZSCORE_XTICK_FONTSIZE    = 10     # condition labels under the z-score panel
CBAR_TICK_FONTSIZE       = 12     # colour-bar tick numbers
COLUMN_HEADER_FONTSIZE   = 12     # "Gene" / "Pathway (source)" column headers
GENE_LABEL_FONTSIZE      = 10     # gene symbols, one per row
PATHWAY_LABEL_FONTSIZE   = 10     # pathway text, one per row

# ---------------------------------------------------------------------
# Long pathway names
# ---------------------------------------------------------------------
# Pathway terms are often far too long for one line, so they are wrapped onto
# several lines (same idea as the GSEA heatmap). Unlike the GSEA heatmap the
# row height here is fixed (one row per protein), so the number of lines is
# capped — anything longer is truncated with an ellipsis rather than running
# into the neighbouring rows.

PATHWAY_WRAP_WIDTH  = 35     # characters per line before wrapping
PATHWAY_MAX_LINES   = 2      # max lines per row (extra text is truncated)
PATHWAY_LINESPACING = 0.95   # line spacing multiplier (<1 tightens the block)


def wrap_pathway(text: str,
                 width: int = PATHWAY_WRAP_WIDTH,
                 max_lines: int = PATHWAY_MAX_LINES) -> str:
    """
    Wrap a long pathway label onto at most `max_lines` lines.

    The trailing source tag ('[R]' / '[GO]') is kept attached to the last line
    so it is never orphaned or cut off by the truncation.
    """
    import textwrap

    text = str(text).strip()
    if not text:
        return ""

    # Detach a trailing source tag so wrapping cannot separate it from the term
    tag = ""
    for t in (TAG_REACTOME, TAG_GO):
        if t and text.endswith(t):
            tag = t
            text = text[: -len(t)].strip()
            break

    lines = textwrap.wrap(text, width=width) or [text]

    if len(lines) > max_lines:
        lines = lines[:max_lines]
        lines[-1] = lines[-1].rstrip(" ,;") + "…"

    if tag:
        # Always attach the tag to the last line. Letting it overflow the wrap
        # width by a few characters looks better than a line containing only
        # "[R]".
        lines[-1] = f"{lines[-1]} {tag}"

    return "\n".join(lines)


# ============================================================
# USER PARAMETERS
# ============================================================
INPUT_FILE = "regulation_table_NAc.tsv"

# Reactome map: 'gene' column holds the GENE SYMBOL, 'pathway' the term.
# Primary annotation source (high specificity, ~60% coverage).
ANNOTATION_FILE = "reactome_map_mouse.tsv"
ANNOTATION_KEY = "gene"             # column in the map that is the gene symbol
ANNOTATION_COL = "pathway"          # column in the map that is the pathway term

# GO Biological Process FALLBACK map (broad coverage). Used only for proteins
# that Reactome does not annotate. Set to None to disable the fallback.
# Build it from the MSigDB mouse GO-BP GMT with go_map.py.
GO_FALLBACK_FILE = "go_map_mouse.tsv"   # or None
GO_KEY = "gene"
GO_COL = "pathway"

# Source tags appended to each annotation so you can see where it came from.
TAG_REACTOME = "[R]"
TAG_GO = "[GO]"

ALPHA = 0.05                        # padj2 significance threshold (both contrasts)
TOP_N = 50                          # top proteins (reverse-regulated) to plot
OUTPUT_PREFIX = "alprazolam_direction_switch"

# How to rank the reverse-regulated proteins for the "top N":
#   "reversal"  -> |log2FC_alp| + |log2FC_kd|  (strongest in BOTH directions)
#   "alp"       -> |log2FC_alp_vs_ctrl| only   (original behaviour)
RANK_BY = "reversal"

# Pretty condition names for display (maps raw group -> label)
GROUP_RENAME = {
    "Control":       "CD",
    "Control_KD":    "KD",
    "Alprazolam":    "Alp+CD",
    "Alprazolam_KD": "Alp+KD",
}
# Column order in the heatmap (raw names; missing ones are skipped)
GROUP_ORDER = ["Control", "Control_KD", "Alprazolam", "Alprazolam_KD"]

# Contrast names in the file
C_ALP_VS_CTRL = "Alprazolam vs Control"
C_KD_VS_ALP   = "Alprazolam_KD vs Alprazolam"

# Pretty labels for comparisons (raw comparison string -> short label)
COMPARISON_RENAME = {
    "Alprazolam vs Control":        "Alp+CD vs CD",
    "Alprazolam_KD vs Alprazolam":  "Alp+KD vs Alp+CD",
    "Alprazolam_KD vs Control_KD":  "Alp+KD vs KD",
    "Control_KD vs Control":        "KD vs CD",
}

# The two figures produced per run:
#   Figure 1 -> the two reversal comparisons
#   Figure 2 -> all four comparisons, in table order
FIG1_COMPARISONS = [C_ALP_VS_CTRL, C_KD_VS_ALP]
FIG2_COMPARISONS = [
    "Control_KD vs Control",
    "Alprazolam vs Control",
    "Alprazolam_KD vs Control_KD",
    "Alprazolam_KD vs Alprazolam",
]
# ============================================================


def clean_cols(df):
    df.columns = [str(c).strip() for c in df.columns]
    return df


def main():
    # ========================================================
    # LOAD MAIN DATA + CLEAN IDENTIFIERS
    # ========================================================
    df = clean_cols(pd.read_csv(INPUT_FILE, sep="\t"))

    # identifier is "SYMBOL~UNIPROT"; keep the UniProt part as the unique id,
    # but also keep the gene symbol (already in display_name) for annotation.
    df["identifier_raw"] = df["identifier"].astype(str)
    df["uniprot"] = df["identifier_raw"].str.split("~").str[-1]
    df["gene_symbol"] = df["display_name"].astype(str).str.strip()
    # unique protein key = UniProt (stable, 1:1); fall back to raw if empty
    df["protein_key"] = df["uniprot"].where(df["uniprot"].str.len() > 0,
                                            df["identifier_raw"])

    # ========================================================
    # LOAD ANNOTATION MAP (Reactome): gene_symbol -> "; "-joined pathways
    # ========================================================
    def load_gene_to_terms(path, key, col):
        """Load a gene->'; '-joined-terms dict from a map file."""
        m = clean_cols(pd.read_csv(path, sep="\t"))
        if key not in m.columns or col not in m.columns:
            raise ValueError(f"{path} must contain '{key}' and '{col}'.")
        m[key] = m[key].astype(str).str.strip()
        m[col] = m[col].astype(str).str.strip()
        m = m[(m[col] != "") & (m[col].str.lower() != "nan")]
        return (
            m.groupby(key)[col]
            .apply(lambda s: "; ".join(sorted(set(v for v in s if v))))
            .to_dict()
        )

    annot_collapsed = load_gene_to_terms(ANNOTATION_FILE, ANNOTATION_KEY, ANNOTATION_COL)

    # GO fallback (optional)
    go_collapsed = {}
    if GO_FALLBACK_FILE:
        try:
            go_collapsed = load_gene_to_terms(GO_FALLBACK_FILE, GO_KEY, GO_COL)
            print(f"Loaded GO fallback: {len(go_collapsed)} genes annotated")
        except FileNotFoundError:
            print(f"  [WARN] GO fallback file '{GO_FALLBACK_FILE}' not found; "
                  f"continuing with Reactome only.")

    def annotate_gene(sym):
        """Reactome first (tagged [R]); GO only if Reactome has nothing ([GO])."""
        r = annot_collapsed.get(sym, "")
        if r:
            return r, TAG_REACTOME
        g = go_collapsed.get(sym, "")
        if g:
            return g, TAG_GO
        return "", ""

    # ========================================================
    # EXTRACT CONTRASTS + MERGE ON PROTEIN
    # ========================================================
    alp = df[df["comparison"] == C_ALP_VS_CTRL].copy()
    kd  = df[df["comparison"] == C_KD_VS_ALP].copy()
    if alp.empty or kd.empty:
        raise SystemExit("One of the required contrasts is missing from the data.")

    keep = ["protein_key", "gene_symbol", "log2FC", "padj2"]
    merged = alp[keep].merge(
        kd[["protein_key", "log2FC", "padj2"]],
        on="protein_key", suffixes=("_alp", "_kd"),
    )

    # ========================================================
    # SIGNIFICANT IN BOTH + REVERSE-REGULATED
    # ========================================================
    sig = merged[(merged["padj2_alp"] < ALPHA) & (merged["padj2_kd"] < ALPHA)].copy()

    reverse = sig[
        ((sig["log2FC_alp"] > 0) & (sig["log2FC_kd"] < 0)) |
        ((sig["log2FC_alp"] < 0) & (sig["log2FC_kd"] > 0))
    ].copy()
    reverse["direction_pattern"] = np.where(
        reverse["log2FC_alp"] > 0, "up_in_alp_down_in_kd", "down_in_alp_up_in_kd"
    )

    # ranking score
    if RANK_BY == "reversal":
        reverse["rank_score"] = reverse["log2FC_alp"].abs() + reverse["log2FC_kd"].abs()
    else:
        reverse["rank_score"] = reverse["log2FC_alp"].abs()
    reverse = reverse.sort_values("rank_score", ascending=False)

    print(f"Significant in both contrasts: {len(sig)}")
    print(f"Reverse-regulated: {len(reverse)} "
          f"(up→down: {(reverse['direction_pattern']=='up_in_alp_down_in_kd').sum()}, "
          f"down→up: {(reverse['direction_pattern']=='down_in_alp_up_in_kd').sum()})")

    top = reverse.head(TOP_N).copy() if TOP_N else reverse.copy()
    print(f"Plotting top {len(top)} proteins")

    # attach pathway annotation: Reactome first, GO fallback, with source tag
    ann_pairs = top["gene_symbol"].map(annotate_gene)
    top["pathway"] = [p for p, _ in ann_pairs]
    top["annotation_source"] = [src for _, src in ann_pairs]

    # save tables
    top.to_csv(f"{OUTPUT_PREFIX}_top{TOP_N}_reverse_regulated.tsv", sep="\t", index=False)
    reverse.to_csv(f"{OUTPUT_PREFIX}_all_reverse_regulated.tsv", sep="\t", index=False)

    n_r = (top["annotation_source"] == TAG_REACTOME).sum()
    n_g = (top["annotation_source"] == TAG_GO).sum()
    n_none = (top["annotation_source"] == "").sum()
    print(f"Annotation source — Reactome: {n_r}, GO fallback: {n_g}, none: {n_none}")

    # ========================================================
    # EXPRESSION MATRIX (4 CONDITIONS), z-scored per protein
    # ========================================================
    expr_rows = []
    for _, r in df.iterrows():
        expr_rows.append({"protein_key": r["protein_key"], "group": r["group1"],
                          "mean": r["mean(group1)"]})
        expr_rows.append({"protein_key": r["protein_key"], "group": r["group2"],
                          "mean": r["mean(group2)"]})
    expr_wide = pd.DataFrame(expr_rows).pivot_table(
        index="protein_key", columns="group", values="mean", aggfunc="mean"
    )

    groups = [g for g in GROUP_ORDER if g in expr_wide.columns]
    order = [k for k in top["protein_key"] if k in expr_wide.index]
    if not order:
        print("No top proteins have expression values; tables saved, skipping heatmap.")
        return

    heat = expr_wide.loc[order, groups]
    z = heat.sub(heat.mean(axis=1), axis=0).div(
        heat.std(axis=1).replace(0, np.nan), axis=0).fillna(0.0)

    # labels
    sym_map = top.set_index("protein_key")["gene_symbol"].to_dict()
    path_map = top.set_index("protein_key")["pathway"].to_dict()
    src_map = top.set_index("protein_key")["annotation_source"].to_dict()
    row_syms = [sym_map.get(k, k) for k in order]

    def first_pathway(s):
        if not isinstance(s, str) or not s.strip():
            return ""
        return s.split(";")[0].strip()

    row_paths = []
    for k in order:
        term = first_pathway(path_map.get(k, ""))
        src = src_map.get(k, "")
        row_paths.append(f"{term} {src}".strip() if term else "")
    col_labels = [GROUP_RENAME.get(g, g) for g in groups]

    # ========================================================
    # Build per-comparison log2FC and padj2 matrices (proteins x comparisons)
    # ========================================================
    # long -> wide: one column per comparison, indexed by protein_key
    lfc_by_comp = (df.pivot_table(index="protein_key", columns="comparison",
                                  values="log2FC", aggfunc="mean"))
    padj_by_comp = (df.pivot_table(index="protein_key", columns="comparison",
                                   values="padj2", aggfunc="mean"))

    # ========================================================
    # PLOT: 3 panels (log2FC heatmap | z-score heatmap | pathway column)
    # ========================================================
    sns.set_style("white")

    def make_three_panel(comparisons, out_png, title_suffix):
        comps = [c for c in comparisons if c in lfc_by_comp.columns]
        if not comps:
            print(f"  [WARN] none of the requested comparisons present; skipping {out_png}")
            return
        lfc = lfc_by_comp.loc[order, comps]
        padj = padj_by_comp.loc[order, comps]
        comp_labels = [COMPARISON_RENAME.get(c, c) for c in comps]

        # padj2 numbers as cell annotations; exact 0 shown as "< 0.0001"
        def fmt_padj(v):
            if pd.isna(v):
                return ""
            if v == 0:
                return "< 0.0001"
            return f"{v:.2g}"
        annot = (padj.applymap(fmt_padj).values if hasattr(padj, "applymap")
                 else padj.map(fmt_padj))

        n = len(order)
        n_comp = len(comps)
        fig_h = max(6.0, 0.32 * n + 2.0)
        # width scales a bit with number of comparison columns
        fig_w = 12 + 1.1 * n_comp
        fig = plt.figure(figsize=(fig_w, fig_h))
        # 4 columns: log2FC | gene labels | z-score | pathway text.
        # Gene symbols sit in the gap BETWEEN the two heatmaps so each row is
        # labelled right where both panels are read. The log2FC colour bar is a
        # small inset on the LEFT of the log2FC panel.
        gs = gridspec.GridSpec(
            1, 4,
            width_ratios=[0.55 * n_comp + 0.8, 0.10, 2, 2],
            wspace=0.45)

        ax_lfc = fig.add_subplot(gs[0])
        ax_gene = fig.add_subplot(gs[1], sharey=ax_lfc)   # gene labels in the gap
        ax_z = fig.add_subplot(gs[2], sharey=ax_lfc)
        ax_ann = fig.add_subplot(gs[3], sharey=ax_lfc)

        # small colour-bar axis to the LEFT of the log2FC panel, matched in
        # height/width to the z-score colour bar (shrink 0.4-ish of panel height)
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cax_lfc = inset_axes(
            ax_lfc, width="4%", height="40%", loc="center left",
            bbox_to_anchor=(-0.10, 0.0, 1, 1), bbox_transform=ax_lfc.transAxes,
            borderpad=0)

        # --- panel 1: log2FC heatmap with padj2 numbers, colour bar on LEFT ---
        vmax = np.nanmax(np.abs(lfc.values)) if np.isfinite(lfc.values).any() else 1.0
        hm1 = sns.heatmap(
            lfc, ax=ax_lfc, cmap="RdBu_r", center=0,
            vmin=-vmax, vmax=vmax,
            xticklabels=comp_labels, yticklabels=False,
            annot=annot, fmt="", annot_kws={"fontsize": HEATMAP_ANNOT_FONTSIZE},
            linewidths=0.4, linecolor="white",
            cbar_ax=cax_lfc,
            cbar_kws={"label": "log2FC"},
        )
        cax_lfc.yaxis.set_ticks_position("left")
        cax_lfc.yaxis.set_label_position("left")
        cax_lfc.tick_params(labelsize=CBAR_TICK_FONTSIZE)
        ax_lfc.tick_params(axis="y", rotation=0, labelsize=GENE_LABEL_FONTSIZE, length=0)
        # horizontal, smaller comparison labels
        ax_lfc.tick_params(axis="x", rotation=0, labelsize=HEATMAP_XTICK_FONTSIZE)
        for t in ax_lfc.get_xticklabels():
            t.set_ha("center")
        ax_lfc.set_xlabel(""); ax_lfc.set_ylabel("")
        # Title size comes from rcParams["axes.titlesize"]
        ax_lfc.set_title("log2FC per comparison\n(padj2 in cells)",
                         fontweight="bold")

        # --- panel 2: z-score heatmap (4 conditions) ---
        # matched small inset colour bar on the RIGHT of this panel, same
        # height/width as the log2FC one so the two legends are equal size.
        cax_z = inset_axes(
            ax_z, width="4%", height="40%", loc="center right",
            bbox_to_anchor=(0.10, 0.0, 1, 1), bbox_transform=ax_z.transAxes,
            borderpad=0)
        hm2 = sns.heatmap(
            z, ax=ax_z, cmap="coolwarm", center=0,
            xticklabels=col_labels, yticklabels=False,
            linewidths=0.4, linecolor="white",
            cbar_ax=cax_z,
            cbar_kws={"label": "Abundance (z-score)"},
        )
        cax_z.tick_params(labelsize=CBAR_TICK_FONTSIZE)
        ax_z.tick_params(axis="x", rotation=0, labelsize=ZSCORE_XTICK_FONTSIZE)
        ax_z.set_xlabel(""); ax_z.set_ylabel("")
        ax_z.set_title("Abundance z-score\n(4 conditions)",
                       fontweight="bold")

        # --- gene labels in the gap between the two heatmaps ---
        ax_gene.set_xlim(0, 1)
        ax_gene.set_ylim(n, 0)
        ax_gene.axis("off")
        ax_gene.text(0.5, -0.6, "Gene", fontsize=COLUMN_HEADER_FONTSIZE, fontweight="bold",
                     va="bottom", ha="center")
        for i, sym in enumerate(row_syms):
            ax_gene.text(0.5, i + 0.5, sym, va="center", ha="center",
                         fontsize=GENE_LABEL_FONTSIZE, fontweight="bold")

        # --- panel: pathway annotation text ---
        ax_ann.set_xlim(0, 1)
        ax_ann.set_ylim(n, 0)
        ax_ann.axis("off")
        px = 0.05
        ax_ann.text(px, -0.6, "Pathway (source)", fontsize=COLUMN_HEADER_FONTSIZE, fontweight="bold",
                    va="bottom")
        for i, txt in enumerate(row_paths):
            ax_ann.text(px, i + 0.5,
                        wrap_pathway(txt) if txt else "—", va="center",
                        fontsize=PATHWAY_LABEL_FONTSIZE,
                        linespacing=PATHWAY_LINESPACING,
                        color="black" if txt else "0.6")

        # Size comes from rcParams["figure.titlesize"]
        fig.suptitle(f"Top reverse-regulated proteins — {title_suffix}",
                     fontweight="bold")
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  {out_png}")

    def make_lfc_gene_pathway(comparisons, out_png, title_suffix):
        """Like make_three_panel but WITHOUT the z-score heatmap:
        log2FC heatmap | gene labels | pathway text."""
        comps = [c for c in comparisons if c in lfc_by_comp.columns]
        if not comps:
            print(f"  [WARN] no comparisons present; skipping {out_png}")
            return
        lfc = lfc_by_comp.loc[order, comps]
        padj = padj_by_comp.loc[order, comps]
        comp_labels = [COMPARISON_RENAME.get(c, c) for c in comps]

        def fmt_padj(v):
            if pd.isna(v):
                return ""
            if v == 0:
                return "< 0.0001"
            return f"{v:.2g}"
        annot = (padj.applymap(fmt_padj).values if hasattr(padj, "applymap")
                 else padj.map(fmt_padj))

        n = len(order)
        n_comp = len(comps)
        fig_h = max(6.0, 0.32 * n + 2.0)
        fig_w = 9 + 1.1 * n_comp
        fig = plt.figure(figsize=(fig_w, fig_h))
        # 3 columns: log2FC | gene labels | pathway text
        gs = gridspec.GridSpec(
            1, 3,
            width_ratios=[0.55 * n_comp + 0.8, 0.25, 2.2],
            wspace=0.18)
        ax_lfc = fig.add_subplot(gs[0])
        ax_gene = fig.add_subplot(gs[1], sharey=ax_lfc)
        ax_ann = fig.add_subplot(gs[2], sharey=ax_lfc)

        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cax_lfc = inset_axes(
            ax_lfc, width="4%", height="40%", loc="center left",
            bbox_to_anchor=(-0.10, 0.0, 1, 1), bbox_transform=ax_lfc.transAxes,
            borderpad=0)

        vmax = np.nanmax(np.abs(lfc.values)) if np.isfinite(lfc.values).any() else 1.0
        sns.heatmap(
            lfc, ax=ax_lfc, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
            xticklabels=comp_labels, yticklabels=False,
            annot=annot, fmt="", annot_kws={"fontsize": HEATMAP_ANNOT_FONTSIZE},
            linewidths=0.4, linecolor="white",
            cbar_ax=cax_lfc, cbar_kws={"label": "log2FC"})
        cax_lfc.yaxis.set_ticks_position("left")
        cax_lfc.yaxis.set_label_position("left")
        cax_lfc.tick_params(labelsize=CBAR_TICK_FONTSIZE)
        ax_lfc.tick_params(axis="x", rotation=0, labelsize=HEATMAP_XTICK_FONTSIZE)
        for t in ax_lfc.get_xticklabels():
            t.set_ha("center")
        ax_lfc.set_xlabel(""); ax_lfc.set_ylabel("")
        ax_lfc.set_title("log2FC per comparison\n(padj2 in cells)",
                         fontweight="bold")

        ax_gene.set_xlim(0, 1); ax_gene.set_ylim(n, 0); ax_gene.axis("off")
        ax_gene.text(0.5, -0.6, "Gene", fontsize=COLUMN_HEADER_FONTSIZE, fontweight="bold",
                     va="bottom", ha="center")
        for i, sym in enumerate(row_syms):
            ax_gene.text(0.5, i + 0.5, sym, va="center", ha="center",
                         fontsize=GENE_LABEL_FONTSIZE, fontweight="bold")

        ax_ann.set_xlim(0, 1); ax_ann.set_ylim(n, 0); ax_ann.axis("off")
        px = 0.05
        ax_ann.text(px, -0.6, "Pathway (source)", fontsize=COLUMN_HEADER_FONTSIZE, fontweight="bold",
                    va="bottom")
        for i, txt in enumerate(row_paths):
            ax_ann.text(px, i + 0.5,
                        wrap_pathway(txt) if txt else "—", va="center",
                        fontsize=PATHWAY_LABEL_FONTSIZE,
                        linespacing=PATHWAY_LINESPACING,
                        color="black" if txt else "0.6")

        # Size comes from rcParams["figure.titlesize"]
        fig.suptitle(f"Top reverse-regulated proteins — {title_suffix}",
                     fontweight="bold")
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  {out_png}")

    print("\nSaved:")
    print(f"  {OUTPUT_PREFIX}_top{TOP_N}_reverse_regulated.tsv")
    print(f"  {OUTPUT_PREFIX}_all_reverse_regulated.tsv")
    make_three_panel(
        FIG1_COMPARISONS,
        f"{OUTPUT_PREFIX}_heatmap_2comparisons.png",
        "reversal comparisons (Alp+CD vs CD, Alp+KD vs Alp+CD)")
    make_three_panel(
        FIG2_COMPARISONS,
        f"{OUTPUT_PREFIX}_heatmap_4comparisons.png",
        "all four comparisons")
    make_lfc_gene_pathway(
        FIG1_COMPARISONS,
        f"{OUTPUT_PREFIX}_heatmap_2comparisons_noZ.png",
        "reversal comparisons (log2FC only)")


if __name__ == "__main__":
    main()
