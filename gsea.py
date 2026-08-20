#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gsea_optimized.py

Optimized standalone pre-ranked GSEA with keyword filtering.

Key improvements:
- Vectorized operations for better performance
- Better memory management
- Enhanced error handling and validation
- Progress indicators for long-running operations
- Parallel processing option for multiple comparisons
- Improved code documentation

Outputs
-------
- gsea_results.csv                 (all tested pathways across comparisons)
- gsea_leading_edge.csv            (leading-edge members for each pathway/comparison)
- gsea_results_filtered.csv        (optional: filtered view used for plots)
- gsea_nes_barplot__<comparison>.png
- gsea_enrichment__<pathway>__<comparison>.png
- gsea_nes_heatmap.png

Usage example
-------------
python gsea_optimized.py \
  --reg-table regulation_table_Amygdala.tsv \
  --pathway-map reactome_map_mouse.tsv \
  --id-col display_name --stat-col log2FC \
  --comparisons "Alprazolam_KD vs Alprazolam, Alprazolam vs Control" \
  --outdir GSEA_out_Amygdala \
  --nperm 500 --min-size 10 --topn 30 --seed 1 \
  --hierarchy reactome_hierarchy_mouse.tsv \
  --ancestors "R-MMU-112316,R-MMU-1428517"
"""

from __future__ import annotations

import argparse
import os
import re
import hashlib
import warnings
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    import seaborn as sns
except ImportError:
    sns = None
    warnings.warn("Seaborn not available. Heatmap will be skipped.")

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    warnings.warn("tqdm not available. Progress bars will be disabled.")


# ---------------------------------------------------------------------
# Global font settings  <<< CHANGE FONT SIZES HERE >>>
# (same block as Volcano_plot.py / pca_samples_advanced.py — keep in sync)
# ---------------------------------------------------------------------

plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size":        10,   # baseline (fallback for anything not set below)
    "axes.labelsize":   14,   # x/y axis labels
    "axes.titlesize":   14,   # per-axes title
    "figure.titlesize": 14,   # fig.suptitle (used by the heatmap)
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
# Heatmap-specific font sizes
# ---------------------------------------------------------------------
# The NES heatmap packs many small cells, so its text does NOT follow the
# global sizes above — 12 pt annotations would overlap and the wrapped pathway
# names would not fit the reserved left margin. Tune these separately.

HEATMAP_ANNOT_FONTSIZE  = 8    # p-value / FDR numbers inside the cells
HEATMAP_XTICK_FONTSIZE  = 10    # comparison labels under the heatmap
HEATMAP_YTICK_FONTSIZE  = 10    # wrapped pathway names on the left
CBAR_TICK_FONTSIZE      = 12    # colour-bar tick numbers


# -----------------------------
# Constants
# -----------------------------

COMPARISON_LABELS = {
    "Alprazolam vs Control": "Alp+CD vs CD",
    "Alprazolam_KD vs Alprazolam": "Alp+KD vs Alp+CD",
    "Alprazolam_KD vs Control_KD": "Alp+KD vs KD",
    "Control_KD vs Control": "KD vs CD",
}


# -----------------------------
# Utilities
# -----------------------------

def pretty_comp(comp: str) -> str:
    """Return abbreviated comparison label for plotting."""
    return COMPARISON_LABELS.get(comp, comp)


def ensure_outdir(outdir: str) -> str:
    """Create output directory if it doesn't exist."""
    Path(outdir).mkdir(parents=True, exist_ok=True)
    return outdir


def read_table(path: str) -> pd.DataFrame:
    """Read TSV or CSV file with automatic format detection."""
    path_lower = str(path).lower()
    try:
        if path_lower.endswith(('.tsv', '.txt')):
            return pd.read_csv(path, sep="\t")
        elif path_lower.endswith('.csv'):
            return pd.read_csv(path)
        else:
            # Try TSV first, then CSV
            try:
                return pd.read_csv(path, sep="\t")
            except Exception:
                return pd.read_csv(path)
    except Exception as e:
        raise ValueError(f"Failed to read file {path}: {str(e)}")


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """
    Benjamini-Hochberg FDR correction.
    
    Optimized version with better handling of edge cases.
    """
    p = np.asarray(pvals, dtype=float)
    n = p.size
    
    if n == 0:
        return np.array([])
    
    q = np.full(n, np.nan, dtype=float)
    ok = np.isfinite(p)
    
    if not np.any(ok):
        return q
    
    p_ok = p[ok]
    n_ok = p_ok.size
    
    # Sort and calculate FDR
    idx = np.argsort(p_ok)
    p_sorted = p_ok[idx]
    ranks = np.arange(1, n_ok + 1, dtype=float)
    
    # Vectorized FDR calculation
    q_sorted = p_sorted * (n_ok / ranks)
    
    # Ensure monotonicity (reverse cumulative minimum)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)
    
    # Map back to original order
    q_ok = np.empty_like(p_ok)
    q_ok[idx] = q_sorted
    q[ok] = q_ok
    
    return q


def stable_hash_int(s: str) -> int:
    """Generate stable integer hash from string for reproducible randomization."""
    h = hashlib.sha256(s.encode("utf-8")).hexdigest()[:16]
    return int(h, 16)


def normalize_id(x: str) -> str:
    """Normalize gene/protein identifiers to uppercase."""
    return str(x).strip().upper()


def sanitize_filename(s: str, max_len: int = 180) -> str:
    """Sanitize string for use in filenames."""
    s = re.sub(r"[^\w\-\.]+", "_", str(s).strip())
    return s[:max_len] if len(s) > max_len else s


def star(p: float) -> str:
    """Convert p-value to significance stars."""
    if not np.isfinite(p):
        return ""
    if p <= 0.001:
        return "***"
    if p <= 0.01:
        return "**"
    if p <= 0.05:
        return "*"
    return ""


# -----------------------------
# GSEA core - Optimized
# -----------------------------

def compute_es_from_hits(
    hit_positions: np.ndarray, 
    hit_weights: np.ndarray, 
    N: int
) -> float:
    """
    Compute enrichment score from hit positions and weights.
    
    Optimized with vectorized operations.
    """
    k = int(hit_positions.size)
    if k == 0 or k == N:
        return np.nan

    w = np.abs(hit_weights).astype(float)
    NR = w.sum()
    
    if not np.isfinite(NR) or NR <= 0:
        return np.nan

    # Vectorized running sum calculation
    p_hit = np.cumsum(w) / NR
    p_miss = (hit_positions - np.arange(k, dtype=float)) / float(N - k)
    running = p_hit - p_miss

    # Find max deviation
    max_es = float(np.nanmax(running))
    min_es = float(np.nanmin(running))
    
    return max_es if abs(max_es) >= abs(min_es) else min_es


def preranked_gsea_one_comparison(
    ranked: pd.Series,
    pathways: Dict[str, List[str]],
    nperm: int = 1000,
    min_size: int = 10,
    max_size: int = 500,
    seed: Optional[int] = 1,
    weight_p: float = 1.0,
    show_progress: bool = True,
) -> Tuple[pd.DataFrame, Dict[str, np.ndarray]]:
    """
    Perform pre-ranked GSEA for one comparison.
    
    Optimized with better memory management and optional progress bars.
    """
    rng = np.random.default_rng(None if seed is None else int(seed))

    # Prepare ranked list
    ranked = ranked.dropna()
    ranked.index = ranked.index.map(normalize_id)
    ranked = ranked.groupby(level=0).mean()
    ranked = ranked.sort_values(ascending=False)

    ids = ranked.index.to_numpy()
    stats = ranked.to_numpy(dtype=float)
    N = stats.size
    
    if N < 2:
        return pd.DataFrame(), {"ranked_ids": ids, "ranked_stats": stats}

    # Pre-compute weights
    w_all = np.power(np.abs(stats), weight_p)
    id_to_idx = {gid: i for i, gid in enumerate(ids)}

    # Filter pathways by size and overlap
    p_members_idx: Dict[str, np.ndarray] = {}
    p_members_id: Dict[str, np.ndarray] = {}

    for p, members in pathways.items():
        m_norm = [normalize_id(x) for x in members if str(x).strip()]
        idxs = [id_to_idx[m] for m in m_norm if m in id_to_idx]
        
        if not idxs:
            continue
            
        idxs = np.array(sorted(set(idxs)), dtype=int)
        k = idxs.size
        
        if k < min_size or k > max_size:
            continue
            
        p_members_idx[p] = idxs
        p_members_id[p] = ids[idxs]

    if not p_members_idx:
        return pd.DataFrame(), {"ranked_ids": ids, "ranked_stats": stats}

    print(f"  → Testing {len(p_members_idx)} pathways (size {min_size}-{max_size})")

    # Compute observed ES and leading edge
    obs = []
    for p, idxs in p_members_idx.items():
        hit_pos = idxs.astype(float)
        hit_w = w_all[idxs]
        es = compute_es_from_hits(hit_pos, hit_w, N)

        # Leading edge calculation
        k = idxs.size
        w = np.abs(hit_w)
        NR = w.sum()
        
        if not np.isfinite(es) or NR <= 0:
            le = ""
        else:
            p_hit = np.cumsum(w) / NR
            p_miss = (hit_pos - np.arange(k, dtype=float)) / float(N - k)
            running = p_hit - p_miss
            peak_i = int(np.nanargmax(running)) if es >= 0 else int(np.nanargmin(running))
            le = ";".join(p_members_id[p][: peak_i + 1].tolist())

        obs.append((p, float(es), int(k), le))

    obs_df = pd.DataFrame(obs, columns=["pathway", "ES", "size", "leading_edge"])

    # Permutation testing - Optimized
    print(f"  → Running {nperm} permutations...")
    null_es: Dict[str, np.ndarray] = {
        p: np.empty(nperm, dtype=float) for p in obs_df["pathway"]
    }
    
    base_indices = np.arange(N, dtype=int)
    
    # Progress bar for permutations
    perm_iter = range(nperm)
    if HAS_TQDM and show_progress:
        perm_iter = tqdm(perm_iter, desc="  Permutations", leave=False)
    
    for b in perm_iter:
        # Generate permutation
        perm = rng.permutation(base_indices)
        
        # Compute inverse permutation efficiently
        invperm = np.empty_like(perm)
        invperm[perm] = base_indices
        w_perm = w_all[perm]

        # Compute null ES for each pathway
        for p, idxs in p_members_idx.items():
            pos = invperm[idxs].astype(int)
            order = np.argsort(pos)
            pos_sorted = pos[order].astype(float)
            hit_w = w_perm[pos[order]]
            null_es[p][b] = compute_es_from_hits(pos_sorted, hit_w, N)

    # Calculate NES and p-values
    rows = []
    for _, rr in obs_df.iterrows():
        p = rr["pathway"]
        es = float(rr["ES"])
        k = int(rr["size"])
        le = rr["leading_edge"]

        null = null_es[p]
        null = null[np.isfinite(null)]
        
        if null.size == 0 or not np.isfinite(es):
            pval = np.nan
            nes = np.nan
        else:
            if es >= 0:
                pval = (np.sum(null >= es) + 1.0) / (null.size + 1.0)
                pos_null = null[null >= 0]
                denom = np.mean(pos_null) if pos_null.size > 0 else np.nan
                nes = es / denom if (np.isfinite(denom) and denom != 0) else np.nan
            else:
                pval = (np.sum(null <= es) + 1.0) / (null.size + 1.0)
                neg_null = null[null <= 0]
                denom = np.mean(np.abs(neg_null)) if neg_null.size > 0 else np.nan
                nes = es / denom if (np.isfinite(denom) and denom != 0) else np.nan

        rows.append((p, es, nes, pval, k, le))

    # Create results dataframe
    res = pd.DataFrame(
        rows, 
        columns=["pathway", "ES", "NES", "pval", "size", "leading_edge"]
    )
    res["FDR_q"] = bh_fdr(res["pval"].to_numpy())
    res = res.sort_values(
        ["FDR_q", "pval", "NES", "pathway"], 
        ascending=[True, True, False, True], 
        kind="mergesort"
    ).reset_index(drop=True)

    debug = {
        "ranked_ids": ids, 
        "ranked_stats": stats, 
        "ranked_weights": w_all
    }
    
    return res, debug


# -----------------------------
# IO: pathways and ranked stats
# -----------------------------

def load_pathway_map(
    path: str, 
    pathway_col: Optional[str] = None, 
    member_col: Optional[str] = None
) -> Dict[str, List[str]]:
    """Load pathway-to-member mapping from file."""
    df = read_table(path)
    
    if df.shape[1] < 2:
        raise ValueError("Pathway map must have at least two columns: pathway, member.")
    
    if pathway_col is None:
        pathway_col = df.columns[0]
    if member_col is None:
        member_col = df.columns[1]
    
    if pathway_col not in df.columns or member_col not in df.columns:
        raise ValueError(
            f"Pathway map missing required columns: {pathway_col}, {member_col}\n"
            f"Available columns: {df.columns.tolist()}"
        )
    
    df = df[[pathway_col, member_col]].dropna()
    
    # Group by pathway
    pathways: Dict[str, List[str]] = {}
    for p, sub in df.groupby(pathway_col):
        pathways[str(p)] = sub[member_col].astype(str).tolist()
    
    return pathways


def load_pathway_id_map(
    path: str,
    pathway_col: Optional[str] = None,
    id_col: str = "reactome_id",
) -> Dict[str, str]:
    """Load a pathway-name -> Reactome-ID lookup from the map file.

    Returns {} if the map has no id column (e.g. an older map without IDs), in
    which case hierarchy filtering simply won't be available.
    """
    df = read_table(path)
    if pathway_col is None:
        pathway_col = df.columns[0]
    if id_col not in df.columns or pathway_col not in df.columns:
        return {}
    sub = df[[pathway_col, id_col]].dropna()
    out: Dict[str, str] = {}
    for p, rid in zip(sub[pathway_col].astype(str), sub[id_col].astype(str)):
        rid = rid.strip()
        if rid and p not in out:
            out[p] = rid
    return out


def load_ranked_from_reg_table(
    path: str,
    comparison_col: str = "comparison",
    id_col: str = "gene",
    stat_col: str = "log2FC",
    comparisons: Optional[List[str]] = None,
) -> Dict[str, pd.Series]:
    """Load ranked gene lists from regulation table."""
    df = read_table(path)
    
    # Validate required columns
    for col in (comparison_col, id_col, stat_col):
        if col not in df.columns:
            raise ValueError(
                f"Reg table missing required column '{col}'\n"
                f"Available columns: {df.columns.tolist()}"
            )
    
    df = df[[comparison_col, id_col, stat_col]].copy()
    df[id_col] = df[id_col].astype(str).map(normalize_id)
    df[stat_col] = pd.to_numeric(df[stat_col], errors="coerce")
    df = df.dropna(subset=[comparison_col, id_col, stat_col])

    # Filter by comparisons if specified
    if comparisons:
        keep = set([c.strip() for c in comparisons if c.strip()])
        df = df[df[comparison_col].astype(str).isin(keep)]

    # Group by comparison
    ranked_by_comp: Dict[str, pd.Series] = {}
    for comp, sub in df.groupby(comparison_col):
        s = sub.groupby(id_col)[stat_col].mean()
        ranked_by_comp[str(comp)] = s
    
    return ranked_by_comp


# -----------------------------
# Plotting helpers + filtering
# -----------------------------

def load_hierarchy(path: Optional[str]) -> Optional[pd.DataFrame]:
    """Load the Reactome hierarchy file produced by build_hierarchy.py.

    Expected columns: reactome_id, display_name, parent_id, ancestor_ids,
    ancestor_names. Returns None if no path is given or the file is missing.
    """
    if not path:
        return None
    if not os.path.exists(path):
        print(f"  [WARN] hierarchy file not found: {path}")
        return None
    h = read_table(path)
    needed = {"reactome_id", "ancestor_ids"}
    if not needed.issubset(h.columns):
        print(f"  [WARN] hierarchy file missing columns {needed}; got {list(h.columns)}")
        return None
    return h.fillna("")


def resolve_ancestor_ids(hierarchy: pd.DataFrame, ancestors: List[str]) -> set:
    """Resolve a list of ancestor IDs or names to a set of Reactome IDs."""
    if hierarchy is None:
        return set()
    name_to_id = {
        str(n).strip().lower(): str(i)
        for i, n in zip(hierarchy["reactome_id"], hierarchy["display_name"])
    }
    wanted = set()
    for a in ancestors:
        a = a.strip()
        if not a:
            continue
        if a.upper().startswith("R-"):
            wanted.add(a)
        else:
            sid = name_to_id.get(a.lower())
            if sid:
                wanted.add(sid)
            else:
                print(f"  [WARN] ancestor '{a}' not found in hierarchy by name")
    return wanted


def pathways_under(hierarchy: pd.DataFrame, ancestors: List[str]) -> set:
    """Return the set of Reactome IDs that ARE, or descend from, any ancestor."""
    wanted = resolve_ancestor_ids(hierarchy, ancestors)
    if not wanted:
        return set()
    keep = set(wanted)  # ancestors themselves count
    for rid, anc in zip(hierarchy["reactome_id"], hierarchy["ancestor_ids"]):
        anc_ids = set(str(anc).split(";")) if anc else set()
        if anc_ids & wanted:
            keep.add(str(rid))
    return keep


def filter_results_for_plots(
    df: pd.DataFrame,
    hierarchy: Optional[pd.DataFrame] = None,
    ancestors: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Filter GSEA result rows to pathways under the given ancestor(s).

    Matching is by the `reactome_id` column threaded through from the map. If no
    ancestors are requested, the input is returned unchanged.
    """
    if not ancestors or hierarchy is None:
        return df
    if "reactome_id" not in df.columns:
        print("  [WARN] results have no 'reactome_id' column; cannot filter by "
              "hierarchy. Is your map missing Reactome IDs?")
        return df
    keep_ids = pathways_under(hierarchy, ancestors)
    if not keep_ids:
        print("  [WARN] no pathways resolved under the requested ancestors")
    return df[df["reactome_id"].astype(str).isin(keep_ids)]


def make_nes_barplot(
    res: pd.DataFrame, 
    outdir: str, 
    comparison: str, 
    top_n: int = 20
) -> None:
    """Create horizontal barplot of NES values."""
    if res is None or res.empty:
        return
    
    sub = res.sort_values(
        ["FDR_q", "pval", "NES", "pathway"], 
        ascending=[True, True, False, True], 
        kind="mergesort"
    ).head(top_n).copy()
    
    # Reverse for better visualization
    sub = sub.iloc[::-1]

    fig, ax = plt.subplots(figsize=(9, max(3.5, 0.35 * len(sub))))
    
    # Color bars by direction
    colors = ['red' if x < 0 else 'blue' for x in sub["NES"]]
    ax.barh(sub["pathway"], sub["NES"], color=colors, alpha=0.7)
    ax.axvline(0, color='black', linewidth=1, linestyle='--')
    
    # Font sizes come from the rcParams block at the top of this file
    ax.set_title(f"GSEA NES — {comparison}", fontweight='bold')
    ax.set_xlabel("NES")
    ax.set_ylabel("")
    
    plt.tight_layout()
    
    out_png = os.path.join(outdir, f"gsea_nes_barplot__{sanitize_filename(comparison)}.png")
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  [OK] Wrote: {out_png}")


def plot_enrichment_curve(
    ranked_ids: np.ndarray,
    ranked_stats: np.ndarray,
    pathway_name: str,
    pathway_members: List[str],
    comparison: str,
    outdir: str,
    weight_p: float = 1.0,
    title_extra: str = "",
) -> None:
    """Plot enrichment curve for a single pathway."""
    ids = np.array([normalize_id(x) for x in ranked_ids], dtype=object)
    stat = np.asarray(ranked_stats, dtype=float)
    N = stat.size
    w_all = np.power(np.abs(stat), weight_p)

    members = set(normalize_id(x) for x in pathway_members)
    hit_mask = np.array([gid in members for gid in ids], dtype=bool)
    k = int(hit_mask.sum())
    
    if k == 0 or k == N:
        return

    hit_w = w_all[hit_mask]
    NR = hit_w.sum()
    
    if not np.isfinite(NR) or NR <= 0:
        return

    # Calculate running enrichment score
    P_hit = np.zeros(N, dtype=float)
    P_miss = np.zeros(N, dtype=float)
    P_hit[hit_mask] = hit_w / NR
    P_miss[~hit_mask] = 1.0 / float(N - k)
    running = np.cumsum(P_hit - P_miss)

    es_pos = float(np.max(running))
    es_neg = float(np.min(running))
    es = es_pos if abs(es_pos) >= abs(es_neg) else es_neg

    # Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 5), 
                                    gridspec_kw={'height_ratios': [3, 1]})
    
    # Enrichment curve
    ax1.plot(running, linewidth=2, color='green' if es > 0 else 'red')
    ax1.axhline(0, color='black', linewidth=1, linestyle='--')
    ax1.set_ylabel("Running enrichment score")
    ax1.set_title(
        f"Enrichment curve — {pathway_name}\n{comparison}  ES={es:.3g} {title_extra}".strip()
    )
    ax1.grid(True, alpha=0.3)
    
    # Hit positions
    hit_positions = np.where(hit_mask)[0]
    ax2.vlines(hit_positions, 0, 1, colors='black', linewidths=0.5)
    ax2.set_xlim(0, N)
    ax2.set_ylim(0, 1)
    ax2.set_xlabel("Rank in ordered list")
    ax2.set_ylabel("Hits")
    ax2.set_yticks([])
    
    plt.tight_layout()
    
    out_png = os.path.join(
        outdir, 
        f"gsea_enrichment__{sanitize_filename(pathway_name)}__{sanitize_filename(comparison)}.png"
    )
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()


def make_nes_heatmap(
    all_results: pd.DataFrame,
    outdir: str,
    top_n: int = 30,
    annotate: str = "fdr",
    sig_threshold: float = 0.1,
    p_fmt: str = ".2g",
    display_cutoff: Optional[float] = None,
    display_metric: str = "FDR_q",
) -> None:
    """
    Create heatmap of NES values across comparisons.
    
    Improved version with better handling of missing values and aesthetics.

    display_cutoff: if set, only pathways whose `display_metric` is <= this value
        in AT LEAST ONE comparison are shown. Applied before the top_n cap.
    display_metric: which column the cutoff applies to ('FDR_q' or 'pval').
    """
    if all_results is None or all_results.empty:
        print("  [WARN] No results for heatmap")
        return
    
    if sns is None:
        print("  [WARN] seaborn not available; skipping heatmap.")
        return

    df = all_results.copy()
    df["FDR_q"] = pd.to_numeric(df["FDR_q"], errors="coerce")
    df["pval"] = pd.to_numeric(df["pval"], errors="coerce")
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")

    # Optional significance cutoff: keep only pathways passing the threshold in
    # at least one comparison (before the top_n selection below).
    if display_cutoff is not None:
        if display_metric not in df.columns:
            print(f"  [WARN] display_metric '{display_metric}' not in results; "
                  f"skipping cutoff.")
        else:
            passing = (
                df.loc[df[display_metric] <= display_cutoff, "pathway"]
                .dropna()
                .unique()
            )
            n_before = df["pathway"].nunique()
            df = df[df["pathway"].isin(passing)]
            print(f"  Heatmap cutoff {display_metric} <= {display_cutoff}: "
                  f"{df['pathway'].nunique()}/{n_before} pathways kept")
            if df.empty:
                print("  [WARN] No pathways pass the heatmap display cutoff; "
                      "skipping heatmap.")
                return

    # Sort and select top pathways
    df_sorted = df.sort_values(
        by=["comparison", "FDR_q", "pval", "NES", "pathway"],
        ascending=[True, True, True, False, True],
        kind="mergesort",
    )

    top_paths = (
        df_sorted.groupby("comparison", sort=True, group_keys=False)
        .head(top_n)["pathway"]
        .dropna()
        .unique()
    )
    top_paths = sorted(top_paths)

    # Create NES matrix
    heat_raw = (
        df[df["pathway"].isin(top_paths)]
        .pivot(index="pathway", columns="comparison", values="NES")
        .reindex(index=top_paths)
    ).sort_index(axis=1)

    # Apply pretty labels
    col_map = {c: pretty_comp(c) for c in heat_raw.columns}
    heat = heat_raw.rename(columns=col_map)

    # Prepare annotations
    annot_tbl = None
    if annotate != "none":
        if annotate == "stars":
            q_raw = (
                df[df["pathway"].isin(top_paths)]
                .pivot(index="pathway", columns="comparison", values="FDR_q")
                .reindex(index=top_paths)
            ).sort_index(axis=1)

            q = q_raw.rename(columns=col_map).reindex(
                index=heat.index, columns=heat.columns
            )
            q = q.apply(pd.to_numeric, errors="coerce")
            
            annot_tbl = q.map(
                lambda x: star(x) if (pd.notna(x) and x <= sig_threshold) else ""
            )

        elif annotate in ("fdr", "pval"):
            metric = "FDR_q" if annotate == "fdr" else "pval"
            vals_raw = (
                df[df["pathway"].isin(top_paths)]
                .pivot(index="pathway", columns="comparison", values=metric)
                .reindex(index=top_paths)
            ).sort_index(axis=1)

            vals = vals_raw.rename(columns=col_map).reindex(
                index=heat.index, columns=heat.columns
            )
            vals = vals.apply(pd.to_numeric, errors="coerce")
            vals_sig = vals.where(vals <= sig_threshold)
            
            annot_tbl = vals_sig.map(
                lambda x: format(float(x), p_fmt) if pd.notna(x) else ""
            )

    # Wrap long pathway (y-axis) names onto multiple lines so they don't run off
    # the figure. The index holds full pathway NAMES (never Reactome IDs), so the
    # heatmap always shows readable names.
    import textwrap
    wrap_width = 45
    wrapped_index = {
        lbl: "\n".join(textwrap.wrap(str(lbl), width=wrap_width)) or str(lbl)
        for lbl in heat.index
    }
    heat = heat.rename(index=wrapped_index)
    if annot_tbl is not None:
        annot_tbl = annot_tbl.rename(index=wrapped_index)

    # ------------------------------------------------------------------
    # Fixed cell size: each cell is a constant CELL_W_IN x CELL_H_IN inches,
    # regardless of how many comparisons or pathways there are. The figure
    # grows to fit; the cells never stretch or shrink. Set the two equal for
    # squares, or different for rectangles (e.g. wide-and-short cells).
    # ------------------------------------------------------------------
    n_rows = heat.shape[0]
    n_cols = heat.shape[1]

    CELL_W_IN = 1.25        # width of each cell (inches)
    CELL_H_IN = 0.25        # height of each cell (inches)
    LABEL_W_IN = 5.5        # left margin reserved for (wrapped) pathway names
    LABEL_H_IN = 1.6        # bottom margin reserved for comparison labels
    CBAR_W_IN = 1.0         # right margin for the colour bar
    TITLE_H_IN = 1.0        # top margin for the title

    grid_w = n_cols * CELL_W_IN
    grid_h = n_rows * CELL_H_IN
    fig_w = LABEL_W_IN + grid_w + CBAR_W_IN
    fig_h = TITLE_H_IN + grid_h + LABEL_H_IN

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    sns.heatmap(
        heat,
        cmap="vlag",
        center=0,
        cbar_kws={"label": "NES"},
        annot=annot_tbl if annot_tbl is not None else False,
        fmt="",
        annot_kws={"fontsize": HEATMAP_ANNOT_FONTSIZE},
        linewidths=0.5,
        linecolor='gray',
        ax=ax
    )

    # Pin the axes rectangle so the cell area is exactly grid_w x grid_h inches
    # and the label/title margins are constant across different-sized heatmaps.
    left = LABEL_W_IN / fig_w
    bottom = LABEL_H_IN / fig_h
    width = grid_w / fig_w
    height = grid_h / fig_h
    ax.set_position([left, bottom, width, height])

    # Move the colour bar into the reserved right margin (set_position above only
    # moves the main axes, so the cbar must be repositioned to match). Make it
    # the SAME height as the grid, and tick every integer (so ±1, ±2, ... show).
    if ax.collections and ax.collections[0].colorbar is not None:
        cbar = ax.collections[0].colorbar
        cbar_ax = cbar.ax
        cbar_left = left + width + 0.3 / fig_w
        cbar_w = 0.22 / fig_w
        cbar_ax.set_position([cbar_left, bottom, cbar_w, height])

        # Integer ticks spanning the colour scale, always including +/-1.
        import numpy as _np
        vmin, vmax = ax.collections[0].get_clim()
        lo = int(_np.floor(vmin))
        hi = int(_np.ceil(vmax))
        ticks = list(range(lo, hi + 1))
        for t in (-1, 1):                 # guarantee +/-1 even if rounding skips them
            if vmin <= t <= vmax and t not in ticks:
                ticks.append(t)
        ticks = sorted(t for t in ticks if vmin <= t <= vmax)
        if ticks:
            cbar.set_ticks(ticks)
            cbar_ax.tick_params(labelsize=CBAR_TICK_FONTSIZE)

    ax.set_xlabel("")
    ax.set_ylabel("")

    # Smaller, non-overlapping x-axis (comparison) labels; stagger every other
    # one downward so neighbouring long names don't collide.
    xlabels = [t.get_text() for t in ax.get_xticklabels()]
    ax.set_xticklabels(xlabels, fontsize=HEATMAP_XTICK_FONTSIZE,
                       rotation=0, ha="center")
    #for i, tick in enumerate(ax.get_xticklabels()):
     #   if i % 2 == 1:
      #      tick.set_y(tick.get_position()[1] - 0.03)

    ax.tick_params(axis="y", labelsize=HEATMAP_YTICK_FONTSIZE)
    
    title = "GSEA NES heatmap (top pathways per contrast)"
    if annotate != "none":
        title += f"\nAnnot: {annotate} (≤ {sig_threshold})"
    # Center the title over the whole figure rather than just the grid.
    fig.suptitle(title, fontweight='bold',
                 y=1 - (TITLE_H_IN / fig_h) * 0.35)   # size: rcParams["figure.titlesize"]
    
    out_png = os.path.join(outdir, "gsea_nes_heatmap.png")
    # NOTE: no tight_layout()/bbox_inches='tight' here — those would rescale the
    # carefully fixed cell geometry. Save at the exact figure size instead.
    plt.savefig(out_png, dpi=300)
    plt.close()
    
    print(f"  [OK] Wrote: {out_png}")


# -----------------------------
# Main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Optimized pre-ranked GSEA with keyword filtering.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python gsea_optimized.py \\
    --reg-table regulation_table.tsv \\
    --pathway-map reactome_map.tsv \\
    --id-col gene --stat-col log2FC \\
    --outdir GSEA_output \\
    --nperm 1000 --topn 30
        """
    )
    
    # Input/output
    ap.add_argument("--reg-table", required=True, 
                    help="Regulation/ranking table (TSV/CSV)")
    ap.add_argument("--pathway-map", required=True, 
                    help="Pathway map (TSV/CSV) with pathway + member columns")
    ap.add_argument("--outdir", required=True, 
                    help="Output directory")
    
    # Column specifications
    ap.add_argument("--comparison-col", default="comparison", 
                    help="Column name for contrast labels in reg table")
    ap.add_argument("--id-col", default="gene", 
                    help="Column name for IDs in reg table (gene/protein)")
    ap.add_argument("--stat-col", default="log2FC", 
                    help="Column name for ranking statistic in reg table")
    ap.add_argument("--comparisons", default="", 
                    help="Comma-separated list of comparisons to run (optional)")
    ap.add_argument("--pathway-col", default="", 
                    help="Column name for pathway in map (optional; default first col)")
    ap.add_argument("--member-col", default="", 
                    help="Column name for member ID in map (optional; default second col)")

    # GSEA parameters
    ap.add_argument("--nperm", type=int, default=1000, 
                    help="Number of permutations (default 1000)")
    ap.add_argument("--min-size", type=int, default=10, 
                    help="Min overlap size for pathway (default 10)")
    ap.add_argument("--max-size", type=int, default=500, 
                    help="Max overlap size for pathway (default 500)")
    ap.add_argument("--topn", type=int, default=30, 
                    help="Top pathways per contrast for plots (default 30)")
    ap.add_argument("--seed", type=int, default=1, 
                    help="Base random seed; set -1 for non-deterministic")
    
    # Plotting parameters
    ap.add_argument("--heatmap-annot", default="fdr", 
                    choices=["fdr", "pval", "stars", "none"], 
                    help="Heatmap annotations")
    ap.add_argument("--sig-threshold", type=float, default=0.05, 
                    help="Annotation threshold for heatmap")
    ap.add_argument("--p-fmt", default=".2g", 
                    help="Format for p/q annotation")
    ap.add_argument("--heatmap-cutoff", type=float, default=None,
                    help="If set, the heatmap shows only pathways meeting this "
                         "significance threshold (on --heatmap-cutoff-metric) in "
                         "at least one comparison. E.g. 0.1")
    ap.add_argument("--heatmap-cutoff-metric", default="FDR_q",
                    choices=["FDR_q", "pval"],
                    help="Metric the --heatmap-cutoff applies to (default FDR_q)")
    
    # Filtering (by Reactome hierarchy)
    ap.add_argument("--ancestors", default="",
                    help="Comma-separated Reactome ancestor IDs or names. Plots/"
                         "filtered output keep only pathways that ARE or descend "
                         "from these. E.g. 'R-MMU-112316' (Neuronal System) or "
                         "'R-MMU-112316,R-MMU-1428517'. Empty = keep all.")
    ap.add_argument("--hierarchy", default="",
                    help="Path to reactome_hierarchy_mouse.tsv (from "
                         "build_hierarchy.py). Required when --ancestors is set.")
    
    # Other options
    ap.add_argument("--no-progress", action="store_true",
                    help="Disable progress bars")

    args = ap.parse_args()
    
    # Setup
    outdir = ensure_outdir(args.outdir)
    show_progress = not args.no_progress and HAS_TQDM
    
    # Parse hierarchy filtering options
    ancestors = [a.strip() for a in args.ancestors.split(",") if a.strip()]
    hierarchy = load_hierarchy(args.hierarchy) if ancestors else None
    if ancestors and hierarchy is None:
        print("  [WARN] --ancestors set but no usable --hierarchy; "
              "filtering disabled (all pathways kept).")
        ancestors = []

    comps = [c.strip() for c in args.comparisons.split(",") if c.strip()] if args.comparisons else None
    pathway_col = args.pathway_col.strip() or None
    member_col = args.member_col.strip() or None

    # Load pathway map
    print("\n" + "="*60)
    print("LOADING PATHWAY MAP")
    print("="*60)
    pathways = load_pathway_map(args.pathway_map, pathway_col=pathway_col, member_col=member_col)
    print(f"✓ Loaded {len(pathways)} pathways")

    # Pathway-name -> Reactome ID lookup (for hierarchy filtering); {} if absent
    pathway_ids = load_pathway_id_map(args.pathway_map, pathway_col=pathway_col)
    if ancestors and not pathway_ids:
        print("  [WARN] map has no 'reactome_id' column; hierarchy filtering "
              "will keep all pathways. Rebuild the map with reactome_map.py.")

    # Load ranked lists
    print("\n" + "="*60)
    print("LOADING RANKED LISTS")
    print("="*60)
    ranked_by_comp = load_ranked_from_reg_table(
        args.reg_table,
        comparison_col=args.comparison_col,
        id_col=args.id_col,
        stat_col=args.stat_col,
        comparisons=comps,
    )
    
    if not ranked_by_comp:
        raise SystemExit("✗ No comparisons found to run.")
    
    print(f"✓ Found {len(ranked_by_comp)} comparison(s): {list(ranked_by_comp.keys())}")

    # Run GSEA
    print("\n" + "="*60)
    print("RUNNING GSEA")
    print("="*60)
    
    all_rows = []
    leading_edge_rows = []

    for i, (comp, ranked) in enumerate(ranked_by_comp.items(), 1):
        print(f"\n[{i}/{len(ranked_by_comp)}] Processing: {comp}")
        print(f"  Genes in ranked list: {len(ranked)}")
        
        # Generate reproducible seed for this comparison
        base_seed = None if (args.seed is None or int(args.seed) < 0) else int(args.seed)
        comp_seed = None if base_seed is None else (base_seed + (stable_hash_int(comp) % 1_000_000_000))

        res, debug = preranked_gsea_one_comparison(
            ranked=ranked,
            pathways=pathways,
            nperm=int(args.nperm),
            min_size=int(args.min_size),
            max_size=int(args.max_size),
            seed=comp_seed,
            weight_p=1.0,
            show_progress=show_progress,
        )
        
        if res.empty:
            print(f"  [WARN] No pathway met size filters for '{comp}'")
            continue

        res["comparison"] = comp
        res["reactome_id"] = res["pathway"].map(pathway_ids).fillna("")
        all_rows.append(res)

        # Leading edge
        le = res[["pathway", "comparison", "leading_edge", "NES", "FDR_q", "pval", "size"]].copy()
        leading_edge_rows.append(le)

        # Filter for plotting
        plot_res = filter_results_for_plots(res, hierarchy=hierarchy, ancestors=ancestors)
        
        if plot_res.empty:
            print(f"  [INFO] No pathways match hierarchy filter for plots")
            continue

        # Create barplot
        make_nes_barplot(plot_res, outdir, pretty_comp(comp), top_n=min(int(args.topn), 50))

        # Create enrichment curves for top pathways
        top_plot = plot_res.sort_values(
            ["FDR_q", "pval", "NES", "pathway"],
            ascending=[True, True, False, True],
            kind="mergesort"
        ).head(min(int(args.topn), 20))
        
        print(f"  Creating {len(top_plot)} enrichment curves...")
        for _, row in top_plot.iterrows():
            p = row["pathway"]
            members = pathways.get(p, [])
            extra = f" NES={row['NES']:.3g}  q={row['FDR_q']:.3g}"
            plot_enrichment_curve(
                ranked_ids=debug["ranked_ids"],
                ranked_stats=debug["ranked_stats"],
                pathway_name=p,
                pathway_members=members,
                comparison=pretty_comp(comp),
                outdir=outdir,
                title_extra=extra,
            )

    if not all_rows:
        raise SystemExit("\n✗ No GSEA results produced for any comparison.")

    # Save results
    print("\n" + "="*60)
    print("SAVING RESULTS")
    print("="*60)
    
    all_res = pd.concat(all_rows, ignore_index=True)
    out_csv = os.path.join(outdir, "gsea_results.csv")
    all_res.to_csv(out_csv, index=False)
    print(f"✓ Wrote: {out_csv}")

    if leading_edge_rows:
        le_all = pd.concat(leading_edge_rows, ignore_index=True)
        out_le = os.path.join(outdir, "gsea_leading_edge.csv")
        le_all.to_csv(out_le, index=False)
        print(f"✓ Wrote: {out_le}")

    # Filtered results for plots
    viz_res = filter_results_for_plots(all_res, hierarchy=hierarchy, ancestors=ancestors)
    print(f"\nPathways: all={all_res['pathway'].nunique()}, filtered={viz_res['pathway'].nunique()}")
    
    out_filt = os.path.join(outdir, "gsea_results_filtered.csv")
    viz_res.to_csv(out_filt, index=False)
    print(f"✓ Wrote: {out_filt}")

    # Heatmap across contrasts
    print("\n" + "="*60)
    print("CREATING HEATMAP")
    print("="*60)
    
    make_nes_heatmap(
        all_results=viz_res,
        outdir=outdir,
        top_n=int(args.topn),
        annotate=args.heatmap_annot,
        sig_threshold=float(args.sig_threshold),
        p_fmt=args.p_fmt,
        display_cutoff=args.heatmap_cutoff,
        display_metric=args.heatmap_cutoff_metric,
    )

    print("\n" + "="*60)
    print("✓ DONE!")
    print("="*60)


if __name__ == "__main__":
    main()