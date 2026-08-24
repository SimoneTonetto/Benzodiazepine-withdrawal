#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gsea.py

Standalone pre-ranked Gene Set Enrichment Analysis (GSEA) for quantitative
proteomics, with Reactome pathway support and optional restriction of figures
to a chosen branch of the Reactome hierarchy.

For each contrast, proteins are ranked by the chosen statistic, every pathway
meeting a detected-size threshold is tested by permutation, and results are
written alongside NES barplots, enrichment curves and a cross-contrast heatmap.

Behaviour notes
---------------
- Results, barplots and heatmap rows are ranked by BH-adjusted permutation
  p-value (BH_q) by default; use --sort-metric to rank by FDR_q or pval
  instead. Ties are broken by nominal p, then |NES| descending, then pathway
  name, which matters because permutation p-values are floored at
  1 / (nperm + 1) and top pathways frequently share an adjusted value.
- Both BH_q (Benjamini-Hochberg on nominal permutation p-values) and FDR_q
  (GSEA-style empirical FDR from the pooled null NES) are reported.
- Hierarchy filtering via --ancestors is applied only after testing and
  multiple-testing correction; it changes which pathways are plotted, never
  the statistics or the size of the tested pathway universe.
- Ranking ties and per-contrast seeds are handled deterministically, so runs
  are reproducible unless --seed is set to -1.

Outputs
-------
- gsea_results.csv                 (all tested pathways across comparisons)
- gsea_leading_edge.csv            (leading-edge members for each pathway/comparison)
- gsea_results_filtered.csv        (the subset used for plots)
- gsea_nes_barplot__<comparison>.png
- gsea_enrichment__<pathway>__<comparison>.png
- gsea_nes_heatmap.png

Usage example
-------------
python gsea.py \
  --reg-table regulation_table_Amygdala.tsv \
  --pathway-map reactome_map_mouse.tsv \
  --id-col display_name --stat-col log2FC \
  --comparisons "Alprazolam_KD vs Alprazolam, Alprazolam vs Control" \
  --outdir GSEA_out_Amygdala \
  --nperm 1000 --min-size 10 --topn 30 --seed 1 \
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
# Global font settings
# ---------------------------------------------------------------------

plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size":        10,
    "axes.labelsize":   14,
    "axes.titlesize":   14,
    "figure.titlesize": 14,
    "xtick.labelsize":  12,
    "ytick.labelsize":  12,
    "legend.fontsize":  12,
    "legend.title_fontsize": 12,
    "svg.fonttype":     "none",
    "pdf.fonttype":     42,
})

# ---------------------------------------------------------------------
# Heatmap-specific font sizes
# ---------------------------------------------------------------------
HEATMAP_ANNOT_FONTSIZE  = 8
HEATMAP_XTICK_FONTSIZE  = 10
HEATMAP_YTICK_FONTSIZE  = 10
CBAR_TICK_FONTSIZE      = 12


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
            try:
                return pd.read_csv(path, sep="\t")
            except Exception:
                return pd.read_csv(path)
    except Exception as e:
        raise ValueError(f"Failed to read file {path}: {str(e)}")


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg correction of nominal p-values."""
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
    
    idx = np.argsort(p_ok)
    p_sorted = p_ok[idx]
    ranks = np.arange(1, n_ok + 1, dtype=float)
    
    q_sorted = p_sorted * (n_ok / ranks)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)
    
    q_ok = np.empty_like(p_ok)
    q_ok[idx] = q_sorted
    q[ok] = q_ok
    
    return q


def gsea_fdr_from_nes(
    observed_nes: np.ndarray,
    null_nes_by_pathway: Dict[str, np.ndarray],
    pathway_names: List[str],
) -> np.ndarray:
    """
    Calculate empirical GSEA-style FDR from observed and null NES.

    Positive and negative enrichment tails are treated separately.
    The null tail probability is compared with the corresponding
    observed NES tail probability, following the standard GSEA FDR idea.
    """
    obs = np.asarray(observed_nes, dtype=float)
    q = np.full(obs.shape, np.nan, dtype=float)
    valid = np.isfinite(obs)
    if not np.any(valid):
        return q

    names = [pathway_names[i] for i in np.flatnonzero(valid)]

    null_parts = []
    for p in names:
        x = np.asarray(null_nes_by_pathway.get(p, []), dtype=float)
        x = x[np.isfinite(x)]
        if x.size:
            null_parts.append(x)

    if not null_parts:
        return q

    null_all = np.concatenate(null_parts)
    ov = obs[valid]
    qv = np.full(ov.shape, np.nan, dtype=float)

    # Positive NES tail.
    pos = ov >= 0
    if np.any(pos):
        null_pos = null_all[null_all >= 0]
        if null_pos.size:
            idx = np.flatnonzero(pos)
            order = idx[np.argsort(ov[idx])]
            xs = ov[order]

            null_tail = np.array(
                [(np.sum(null_pos >= x) + 1.0) / (null_pos.size + 1.0) for x in xs]
            )
            obs_tail = np.array(
                [np.sum(ov[pos] >= x) / float(np.sum(pos)) for x in xs]
            )
            vals = null_tail / obs_tail
            vals = np.minimum.accumulate(vals[::-1])[::-1]
            qv[order] = vals

    # Negative NES tail.
    neg = ov < 0
    if np.any(neg):
        null_neg = null_all[null_all <= 0]
        if null_neg.size:
            idx = np.flatnonzero(neg)
            order = idx[np.argsort(ov[idx])[::-1]]
            xs = ov[order]

            null_tail = np.array(
                [(np.sum(null_neg <= x) + 1.0) / (null_neg.size + 1.0) for x in xs]
            )
            obs_tail = np.array(
                [np.sum(ov[neg] <= x) / float(np.sum(neg)) for x in xs]
            )
            vals = null_tail / obs_tail
            vals = np.minimum.accumulate(vals[::-1])[::-1]
            qv[order] = vals

    q[valid] = np.clip(qv, 0.0, 1.0)
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


DEFAULT_SORT_METRIC = "BH_q"
SORT_METRIC_LABELS = {
    "BH_q": "BH-adjusted p",
    "FDR_q": "GSEA FDR q",
    "pval": "nominal p",
}


def resolve_sort_metric(df: pd.DataFrame, metric: str = DEFAULT_SORT_METRIC) -> str:
    """Return a usable sorting metric, falling back if the column is absent."""
    if df is not None and metric in df.columns:
        return metric
    for fallback in ("BH_q", "FDR_q", "pval"):
        if df is not None and fallback in df.columns:
            if fallback != metric:
                print(f"  [WARN] sort metric '{metric}' not available; using '{fallback}'")
            return fallback
    return metric


def sort_by_significance(
    df: pd.DataFrame,
    metric: str = DEFAULT_SORT_METRIC,
    group_cols: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Sort results ascending by the chosen significance metric.

    Ties are broken by nominal p-value, then |NES| (descending), then pathway
    name so the ordering is fully deterministic. Missing values sort last.
    """
    if df is None or df.empty:
        return df

    metric = resolve_sort_metric(df, metric)
    work = df.copy()
    def _num(col: str) -> pd.Series:
        if col in work.columns:
            return pd.to_numeric(work[col], errors="coerce")
        return pd.Series(np.nan, index=work.index, dtype=float)

    work["_metric"] = _num(metric)
    work["_pval"] = _num("pval")
    work["_absnes"] = _num("NES").abs()

    group_cols = list(group_cols or [])
    by = group_cols + ["_metric", "_pval", "_absnes", "pathway"]
    asc = [True] * len(group_cols) + [True, True, False, True]

    work = work.sort_values(by=by, ascending=asc, kind="mergesort", na_position="last")
    return work.drop(columns=["_metric", "_pval", "_absnes"]).reset_index(drop=True)


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
# GSEA core
# -----------------------------

def compute_es_from_hits(
    hit_positions: np.ndarray, 
    hit_weights: np.ndarray, 
    N: int
) -> float:
    """Compute enrichment score from hit positions and weights."""
    k = int(hit_positions.size)
    if k == 0 or k == N:
        return np.nan

    w = np.abs(hit_weights).astype(float)
    NR = w.sum()
    
    if not np.isfinite(NR) or NR <= 0:
        return np.nan

    p_hit = np.cumsum(w) / NR
    p_miss = (hit_positions - np.arange(k, dtype=float)) / float(N - k)
    running = p_hit - p_miss

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
    sort_metric: str = DEFAULT_SORT_METRIC,
) -> Tuple[pd.DataFrame, Dict[str, np.ndarray]]:
    """Perform pre-ranked GSEA for one comparison."""
    rng = np.random.default_rng(None if seed is None else int(seed))

    # Prepare ranked list with duplicate aggregation
    ranked = ranked.dropna()
    ranked.index = ranked.index.map(normalize_id)
    ranked = ranked.groupby(level=0).mean()

    # Create DataFrame to perform multi-key deterministic tie-breaking
    df_ranked = pd.DataFrame({'stat': ranked.values, 'id': ranked.index})
    df_ranked['abs_stat'] = df_ranked['stat'].abs()
    
    # Sort deterministically: primary by stat (descending), secondary by abs(stat), tertiary by ID
    df_ranked = df_ranked.sort_values(
        by=['stat', 'abs_stat', 'id'],
        ascending=[False, False, True],
        kind='mergesort'
    ).reset_index(drop=True)

    ids = df_ranked['id'].to_numpy()
    stats = df_ranked['stat'].to_numpy(dtype=float)
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
            
            # Corrected leading edge calculation for positive vs negative ES
            if es >= 0:
                peak_i = int(np.nanargmax(running))
                le_members = p_members_id[p][: peak_i + 1]
            else:
                peak_i = int(np.nanargmin(running))
                le_members = p_members_id[p][peak_i:]
                
            le = ";".join(le_members.tolist())

        obs.append((p, float(es), int(k), le))

    obs_df = pd.DataFrame(obs, columns=["pathway", "ES", "size", "leading_edge"])

    # Permutation testing
    print(f"  → Running {nperm} permutations...")
    null_es: Dict[str, np.ndarray] = {
        p: np.empty(nperm, dtype=float) for p in obs_df["pathway"]
    }
    
    base_indices = np.arange(N, dtype=int)
    
    perm_iter = range(nperm)
    if HAS_TQDM and show_progress:
        perm_iter = tqdm(perm_iter, desc="  Permutations", leave=False)
    
    for b in perm_iter:
        perm = rng.permutation(base_indices)
        invperm = np.empty_like(perm)
        invperm[perm] = base_indices
        w_perm = w_all[perm]

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

    res = pd.DataFrame(
        rows, 
        columns=["pathway", "ES", "NES", "pval", "size", "leading_edge"]
    )
    res["FDR_q"] = gsea_fdr_from_nes(
        observed_nes=res["NES"].to_numpy(),
        null_nes_by_pathway=null_es,
        pathway_names=res["pathway"].tolist(),
    )
    # Keep BH-adjusted nominal permutation p-values separately.
    res["BH_q"] = bh_fdr(res["pval"].to_numpy())
    res = sort_by_significance(res, metric=sort_metric)

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
    
    pathways: Dict[str, List[str]] = {}
    for p, sub in df.groupby(pathway_col):
        pathways[str(p)] = sub[member_col].astype(str).tolist()
    
    return pathways


def load_pathway_id_map(
    path: str,
    pathway_col: Optional[str] = None,
    id_col: str = "reactome_id",
) -> Dict[str, str]:
    """Load a pathway-name -> Reactome-ID lookup from the map file."""
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

    if comparisons:
        keep = set([c.strip() for c in comparisons if c.strip()])
        df = df[df[comparison_col].astype(str).isin(keep)]

    ranked_by_comp: Dict[str, pd.Series] = {}
    for comp, sub in df.groupby(comparison_col):
        s = sub.groupby(id_col)[stat_col].mean()
        ranked_by_comp[str(comp)] = s
    
    return ranked_by_comp


# -----------------------------
# Plotting helpers + filtering
# -----------------------------

def load_hierarchy(path: Optional[str]) -> Optional[pd.DataFrame]:
    """Load Reactome hierarchy file."""
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
    """Resolve ancestor IDs or names to Reactome IDs."""
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
    """Return Reactome IDs descending from any ancestor."""
    wanted = resolve_ancestor_ids(hierarchy, ancestors)
    if not wanted:
        return set()
    keep = set(wanted)
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
    """Filter GSEA result rows to pathways under given ancestor(s)."""
    if not ancestors or hierarchy is None:
        return df
    if "reactome_id" not in df.columns:
        print("  [WARN] results have no 'reactome_id' column; cannot filter by hierarchy.")
        return df
    keep_ids = pathways_under(hierarchy, ancestors)
    if not keep_ids:
        print("  [WARN] no pathways resolved under the requested ancestors")
    return df[df["reactome_id"].astype(str).isin(keep_ids)]


def make_nes_barplot(
    res: pd.DataFrame, 
    outdir: str, 
    comparison: str, 
    top_n: int = 20,
    sort_metric: str = DEFAULT_SORT_METRIC,
) -> None:
    """Create horizontal barplot of NES values, ranked by significance."""
    if res is None or res.empty:
        return
    
    sub = sort_by_significance(res, metric=sort_metric).head(top_n).copy()
    
    # Reverse so the most significant pathway ends up at the top of the axis.
    sub = sub.iloc[::-1]

    fig, ax = plt.subplots(figsize=(9, max(3.5, 0.35 * len(sub))))
    
    colors = ['red' if x < 0 else 'blue' for x in sub["NES"]]
    ax.barh(sub["pathway"], sub["NES"], color=colors, alpha=0.7)
    ax.axvline(0, color='black', linewidth=1, linestyle='--')
    
    metric_used = resolve_sort_metric(res, sort_metric)
    ax.set_title(
        f"GSEA NES — {comparison}\n"
        f"top {len(sub)} by {SORT_METRIC_LABELS.get(metric_used, metric_used)}",
        fontweight='bold'
    )
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
    """Plot enrichment curve for a single pathway using deterministically sorted ranks."""
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

    P_hit = np.zeros(N, dtype=float)
    P_miss = np.zeros(N, dtype=float)
    P_hit[hit_mask] = hit_w / NR
    P_miss[~hit_mask] = 1.0 / float(N - k)
    running = np.cumsum(P_hit - P_miss)

    es_pos = float(np.max(running))
    es_neg = float(np.min(running))
    es = es_pos if abs(es_pos) >= abs(es_neg) else es_neg

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 5), 
                                    gridspec_kw={'height_ratios': [3, 1]})
    
    ax1.plot(running, linewidth=2, color='green' if es > 0 else 'red')
    ax1.axhline(0, color='black', linewidth=1, linestyle='--')
    ax1.set_ylabel("Running enrichment score")
    ax1.set_title(
        f"Enrichment curve — {pathway_name}\n{comparison}  ES={es:.3g} {title_extra}".strip()
    )
    ax1.grid(True, alpha=0.3)
    
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
    sort_metric: str = DEFAULT_SORT_METRIC,
) -> None:
    """Create heatmap of NES values across comparisons, rows ranked by significance."""
    if all_results is None or all_results.empty:
        print("  [WARN] No results for heatmap")
        return
    
    if sns is None:
        print("  [WARN] seaborn not available; skipping heatmap.")
        return

    df = all_results.copy()
    for col in ("FDR_q", "BH_q", "pval", "NES"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    metric_col = resolve_sort_metric(df, sort_metric)

    if display_cutoff is not None:
        if display_metric not in df.columns:
            print(f"  [WARN] display_metric '{display_metric}' not in results; skipping cutoff.")
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
                print("  [WARN] No pathways pass the heatmap display cutoff; skipping heatmap.")
                return

    # Pick the top_n most significant pathways within each contrast.
    df_sorted = sort_by_significance(df, metric=metric_col, group_cols=["comparison"])

    top_paths = (
        df_sorted.groupby("comparison", sort=True, group_keys=False)
        .head(top_n)["pathway"]
        .dropna()
        .unique()
    )

    # Order the heatmap rows by each pathway's best (smallest) value of the
    # sorting metric across all contrasts, so the most significant sit on top.
    sub = df[df["pathway"].isin(top_paths)]
    row_order = (
        sub.groupby("pathway")
        .agg(
            _metric=(metric_col, "min"),
            _pval=("pval", "min"),
            _absnes=("NES", lambda s: s.abs().max()),
        )
        .reset_index()
        .sort_values(
            by=["_metric", "_pval", "_absnes", "pathway"],
            ascending=[True, True, False, True],
            kind="mergesort",
            na_position="last",
        )
    )
    top_paths = row_order["pathway"].tolist()

    heat_raw = (
        df[df["pathway"].isin(top_paths)]
        .pivot(index="pathway", columns="comparison", values="NES")
        .reindex(index=top_paths)
    ).sort_index(axis=1)

    col_map = {c: pretty_comp(c) for c in heat_raw.columns}
    heat = heat_raw.rename(columns=col_map)

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

        elif annotate in ("fdr", "bh", "pval"):
            metric = (
                "FDR_q" if annotate == "fdr"
                else "BH_q" if annotate == "bh"
                else "pval"
            )
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

    import textwrap
    wrap_width = 45
    wrapped_index = {
        lbl: "\n".join(textwrap.wrap(str(lbl), width=wrap_width)) or str(lbl)
        for lbl in heat.index
    }
    heat = heat.rename(index=wrapped_index)
    if annot_tbl is not None:
        annot_tbl = annot_tbl.rename(index=wrapped_index)

    n_rows = heat.shape[0]
    n_cols = heat.shape[1]

    CELL_W_IN = 1.25
    CELL_H_IN = 0.25
    LABEL_W_IN = 5.5
    LABEL_H_IN = 1.6
    CBAR_W_IN = 1.0
    TITLE_H_IN = 1.0

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

    left = LABEL_W_IN / fig_w
    bottom = LABEL_H_IN / fig_h
    width = grid_w / fig_w
    height = grid_h / fig_h
    ax.set_position([left, bottom, width, height])

    if ax.collections and ax.collections[0].colorbar is not None:
        cbar = ax.collections[0].colorbar
        cbar_ax = cbar.ax
        cbar_left = left + width + 0.3 / fig_w
        cbar_w = 0.22 / fig_w
        cbar_ax.set_position([cbar_left, bottom, cbar_w, height])

        import numpy as _np
        vmin, vmax = ax.collections[0].get_clim()
        lo = int(_np.floor(vmin))
        hi = int(_np.ceil(vmax))
        ticks = list(range(lo, hi + 1))
        for t in (-1, 1):
            if vmin <= t <= vmax and t not in ticks:
                ticks.append(t)
        ticks = sorted(t for t in ticks if vmin <= t <= vmax)
        if ticks:
            cbar.set_ticks(ticks)
            cbar_ax.tick_params(labelsize=CBAR_TICK_FONTSIZE)

    ax.set_xlabel("")
    ax.set_ylabel("")

    xlabels = [t.get_text() for t in ax.get_xticklabels()]
    ax.set_xticklabels(xlabels, fontsize=HEATMAP_XTICK_FONTSIZE, rotation=0, ha="center")
    ax.tick_params(axis="y", labelsize=HEATMAP_YTICK_FONTSIZE)
    
    metric_label = SORT_METRIC_LABELS.get(metric_col, metric_col)
    title = f"GSEA NES heatmap (top pathways per contrast, ranked by {metric_label})"
    if annotate != "none":
        title += f"\nAnnot: {annotate} (≤ {sig_threshold})"
    fig.suptitle(title, fontweight='bold', y=1 - (TITLE_H_IN / fig_h) * 0.35)
    
    out_png = os.path.join(outdir, "gsea_nes_heatmap.png")
    plt.savefig(out_png, dpi=300)
    plt.close()
    
    print(f"  [OK] Wrote: {out_png}")


# -----------------------------
# Main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Pre-ranked GSEA for quantitative proteomics, with optional "
                    "Reactome hierarchy filtering of figures.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python gsea.py \\
    --reg-table regulation_table.tsv \\
    --pathway-map reactome_map.tsv \\
    --id-col gene --stat-col log2FC \\
    --outdir GSEA_output \\
    --nperm 1000 --topn 30

  # Restrict figures to two Reactome branches (statistics unaffected):
  python gsea.py \\
    --reg-table regulation_table.tsv \\
    --pathway-map reactome_map.tsv \\
    --id-col gene --stat-col log2FC \\
    --outdir GSEA_output \\
    --hierarchy reactome_hierarchy_mouse.tsv \\
    --ancestors "R-MMU-112316,R-MMU-1428517"
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
    ap.add_argument("--stat-col", default="T-statistics", 
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
                    help="Min pathway members detected in the ranked list (default 10)")
    ap.add_argument("--max-size", type=int, default=500, 
                    help="Max pathway members detected in the ranked list (default 500)")
    ap.add_argument("--topn", type=int, default=30, 
                    help="Top pathways per contrast carried into figures (default 30)")
    ap.add_argument("--seed", type=int, default=1, 
                    help="Base random seed; set -1 for non-deterministic")
    
    # Plotting parameters
    ap.add_argument("--heatmap-annot", default="bh", 
                    choices=["fdr", "bh", "pval", "stars", "none"], 
                    help="Value printed in heatmap cells (default bh)")
    ap.add_argument("--sig-threshold", type=float, default=0.05, 
                    help="Only annotate heatmap cells at or below this value (default 0.05)")
    ap.add_argument("--p-fmt", default=".2g", 
                    help="Format for p/q annotation")
    ap.add_argument("--heatmap-cutoff", type=float, default=None,
                    help="Drop pathways from the heatmap above this value (default: keep all)")
    ap.add_argument("--heatmap-cutoff-metric", default="BH_q",
                    choices=["BH_q", "FDR_q", "pval"],
                    help="Metric the --heatmap-cutoff applies to (default BH_q)")
    ap.add_argument("--sort-metric", default="BH_q",
                    choices=["BH_q", "FDR_q", "pval"],
                    help="Metric used to rank results, barplots and heatmap rows "
                         "(default BH_q = BH-corrected permutation p-value)")
    
    # Filtering
    ap.add_argument("--ancestors", default="",
                    help="Comma-separated Reactome ancestor IDs or names. Filters figures only, after correction; statistics are unaffected.")
    ap.add_argument("--hierarchy", default="",
                    help="Reactome hierarchy file with reactome_id and ancestor_ids columns; required by --ancestors.")
    
    # Other options
    ap.add_argument("--no-progress", action="store_true",
                    help="Disable progress bars")

    args = ap.parse_args()
    
    outdir = ensure_outdir(args.outdir)
    show_progress = not args.no_progress and HAS_TQDM
    
    ancestors = [a.strip() for a in args.ancestors.split(",") if a.strip()]
    hierarchy = load_hierarchy(args.hierarchy) if ancestors else None
    if ancestors and hierarchy is None:
        print("  [WARN] --ancestors set but no usable --hierarchy; filtering disabled.")
        ancestors = []

    comps = [c.strip() for c in args.comparisons.split(",") if c.strip()] if args.comparisons else None
    pathway_col = args.pathway_col.strip() or None
    member_col = args.member_col.strip() or None

    print("\n" + "="*60)
    print("LOADING PATHWAY MAP")
    print("="*60)
    pathways = load_pathway_map(args.pathway_map, pathway_col=pathway_col, member_col=member_col)
    print(f"✓ Loaded {len(pathways)} pathways")

    pathway_ids = load_pathway_id_map(args.pathway_map, pathway_col=pathway_col)

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

    print("\n" + "="*60)
    print("RUNNING GSEA")
    print("="*60)
    
    all_rows = []
    leading_edge_rows = []

    for i, (comp, ranked) in enumerate(ranked_by_comp.items(), 1):
        print(f"\n[{i}/{len(ranked_by_comp)}] Processing: {comp}")
        print(f"  Genes in ranked list: {len(ranked)}")
        
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
            sort_metric=args.sort_metric,
        )
        
        if res.empty:
            print(f"  [WARN] No pathway met size filters for '{comp}'")
            continue

        res["comparison"] = comp
        res["reactome_id"] = res["pathway"].map(pathway_ids).fillna("")
        all_rows.append(res)

        le = res[["pathway", "comparison", "leading_edge", "NES", "FDR_q", "BH_q", "pval", "size"]].copy()
        leading_edge_rows.append(le)

        plot_res = filter_results_for_plots(res, hierarchy=hierarchy, ancestors=ancestors)
        
        if plot_res.empty:
            print(f"  [INFO] No pathways match hierarchy filter for plots")
            continue

        make_nes_barplot(
            plot_res, outdir, pretty_comp(comp),
            top_n=min(int(args.topn), 50),
            sort_metric=args.sort_metric,
        )

        top_plot = sort_by_significance(
            plot_res, metric=args.sort_metric
        ).head(min(int(args.topn), 20))
        
        print(f"  Creating {len(top_plot)} enrichment curves...")
        for _, row in top_plot.iterrows():
            p = row["pathway"]
            members = pathways.get(p, [])
            extra = (
                f" NES={row['NES']:.3g}"
                f"  BH q={row['BH_q']:.3g}"
                f"  GSEA q={row['FDR_q']:.3g}"
            )
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

    print("\n" + "="*60)
    print("SAVING RESULTS")
    print("="*60)
    
    all_res = pd.concat(all_rows, ignore_index=True)
    all_res = sort_by_significance(
        all_res, metric=args.sort_metric, group_cols=["comparison"]
    )
    out_csv = os.path.join(outdir, "gsea_results.csv")
    all_res.to_csv(out_csv, index=False)
    print(f"✓ Wrote: {out_csv}")

    if leading_edge_rows:
        le_all = pd.concat(leading_edge_rows, ignore_index=True)
        le_all = sort_by_significance(
            le_all, metric=args.sort_metric, group_cols=["comparison"]
        )
        out_le = os.path.join(outdir, "gsea_leading_edge.csv")
        le_all.to_csv(out_le, index=False)
        print(f"✓ Wrote: {out_le}")

    viz_res = filter_results_for_plots(all_res, hierarchy=hierarchy, ancestors=ancestors)
    viz_res = sort_by_significance(
        viz_res, metric=args.sort_metric, group_cols=["comparison"]
    )
    print(f"\nPathways: all={all_res['pathway'].nunique()}, filtered={viz_res['pathway'].nunique()}")
    
    out_filt = os.path.join(outdir, "gsea_results_filtered.csv")
    viz_res.to_csv(out_filt, index=False)
    print(f"✓ Wrote: {out_filt}")

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
        sort_metric=args.sort_metric,
    )

    print("\n" + "="*60)
    print("✓ DONE!")
    print("="*60)


if __name__ == "__main__":
    main()
