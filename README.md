Benzodiazepine withdrawal — proteomics analysis
Code used for the GSEA analyses and heatmap figures in [paper title / DOI].
Two standalone scripts, both reading the same long-format differential-abundance
table:
`gsea.py` — pre-ranked GSEA against Reactome: permutation testing, BH
correction, and NES barplots, enrichment curves and a cross-contrast heatmap.
`TOP.py` — proteins whose direction of change is reversed by the
knockdown, plotted as annotated log2FC / z-score heatmaps.
They share a styling block, so figures from the two sit together in a panel
without restyling.
Install
```bash
git clone https://github.com/SimoneTonetto/Benzodiazepine-withdrawal.git
cd Benzodiazepine-withdrawal
pip install -r requirements.txt
```
Python 3.9+. Needs numpy, pandas, matplotlib and seaborn; tqdm is optional
(progress bars in `gsea.py`).
Input files
Regulation table — long format, one row per protein per contrast.
`gsea.py` needs `comparison`, an ID column and a ranking statistic (all
configurable). `TOP.py` needs those plus `identifier` (`SYMBOL~UNIPROT`),
`padj2`, `group1`, `group2`, `mean(group1)` and `mean(group2)`, and validates
them up front.
```
identifier	display_name	comparison	log2FC	padj2	group1	group2	mean(group1)	mean(group2)
Calm1~P0DP26	Calm1	Alprazolam vs Control	-1.84	0.0031	Alprazolam	Control	21403.5	76812.1
```
Pathway map — one row per pathway–member pair. A single file with columns
`pathway`, `gene`, `reactome_id` serves both scripts: `TOP.py` reads it as-is,
`gsea.py` needs `--pathway-col pathway --member-col gene`. The `reactome_id`
column is used only by `--ancestors`.
```
pathway	gene	reactome_id
Respiratory electron transport	Ndufa4	R-MMU-611105
```
Reactome hierarchy (optional, `gsea.py --ancestors` only) — needs
`reactome_id` and semicolon-separated `ancestor_ids`; a `display_name` column
also lets you name ancestors instead of using IDs.
---
`gsea.py`
```bash
python gsea.py \
  --reg-table regulation_table_Amygdala.tsv \
  --pathway-map reactome_map_mouse.tsv \
  --id-col display_name --stat-col log2FC \
  --outdir GSEA_out_Amygdala \
  --nperm 1000 --topn 30 --heatmap-topn 10 --seed 1
```
Run `python gsea.py --help` for the full flag list. Outputs:
`gsea_results.csv` (every pathway tested), `gsea_leading_edge.csv`,
`gsea_results_filtered.csv` (the subset plotted), plus
`gsea_nes_barplot__<contrast>.png`, `gsea_enrichment__<pathway>__<contrast>.png`
and `gsea_nes_heatmap.png`.
Method. Proteins are ranked per contrast (duplicates averaged, ties broken
deterministically). Pathways are kept when their members detected in that
ranked list number 10–500 — detected size, not annotated size, which matters
at proteomic coverage. An enrichment score is taken as the maximum absolute
deviation of a running sum weighted by |statistic| (exponent p = 1). The null
comes from permuting protein labels with pathway membership fixed; nominal
p-values use `(x + 1) / (nperm + 1)`. ES is normalised by the mean of the
same-signed null (NES), and p-values are BH-corrected within each contrast.
Two q-values are reported. `BH_q` (Benjamini–Hochberg on the permutation
p-values) is the default for ranking and the one to report — it has a formal
guarantee and stays conservative with discrete p-values. `FDR_q` (GSEA-style
empirical FDR from the pooled null NES) is kept for reference; it is a
heuristic, tends anti-conservative, and needs monotonicity enforced. Switch
with `--sort-metric`.
Row counts in the heatmap. `--topn` selects per contrast and the heatmap
plots the union, so `--topn 10` across four contrasts can give up to 40 rows.
`--heatmap-topn N` caps the total, ranking pathways by their best q across
contrasts. A pathway strong in one contrast only will therefore survive the
cut — usually what you want, but it means some rows look near-blank elsewhere.
Ties. Permutation p-values are floored at `1 / (nperm + 1)`, so top
pathways often share an identical q. Ties break on nominal p, then |NES|
descending, then name. Raising `--nperm` is the real fix.
Hierarchy filtering. `--ancestors` restricts figures only, after
correction. It never changes q-values or the size of the tested universe —
filtering before correction would inflate significance.
---
`TOP.py`
```bash
python TOP.py                                    # settings from the file header
python TOP.py --input regulation_table_NAc.tsv \
              --output-prefix switch_NAc --top-n 50
```
Every flag defaults to the matching constant in the USER PARAMETERS block, so a
bare run reproduces the in-file configuration. Outputs two TSVs (the plotted
subset and all reverse-regulated proteins) and three figures: log2FC + genes +
z-score + pathway for the two reversal contrasts, the same for all four, and a
version without the z-score panel.
Selection. A protein is kept when it is significant in both contrasts
(`padj2 < ALPHA` each) and its log2FCs have opposite signs. Survivors
rank by effect size — `--rank-by reversal` uses |log2FC_alp| + |log2FC_kd|,
`alp` uses the drug contrast alone. `padj2` is a gate here, not a sort key.
Annotation. One term per protein: Reactome first (tagged `[R]`), GO-BP only
where Reactome is silent (`[GO]`), so the tag shows each label's specificity.
Where several terms match, the first alphabetically is shown — a display
convenience, not a claim about primary function. Full term lists stay in the
TSVs. Disable the fallback with `--go-fallback none`.
Figure tuning. Font sizes and pathway wrapping (`PATHWAY_WRAP_WIDTH`,
`PATHWAY_MAX_LINES`) are constants at the top of the file. Rows are
fixed-height, so long terms are truncated with an ellipsis, keeping the
`[R]`/`[GO]` tag attached. `GROUP_RENAME`, `COMPARISON_RENAME` and
`GROUP_ORDER` control display names and column order.
---
Interpretation notes
The GSEA null permutes gene labels, not samples. It assumes proteins are
independent; complex and pathway co-members are not, and quantification adds
more correlation. Both q-values are anti-conservative in absolute terms and are
best read as a ranking device than as calibrated error rates. This is shared
with GSEAPreranked and `fgsea` on pre-ranked input — with the underlying
intensity matrix, a sample-permutation null is better calibrated.
`TOP.py` selection is not an interaction test. Requiring significance in
two contrasts with opposite signs is a conjunction of two separate tests, and
carries no calibrated error rate for the reversal itself. Treat the output as a
prioritised shortlist; for a formal claim, fit an interaction term in the
underlying model.
Ranking on log2FC ignores variance, so a large change in a noisy protein
ranks alongside a well-measured one. A moderated t-statistic is preferable
where available.
Performance. The permutation loop is O(nperm × n_pathways) in Python;
`fgsea` is faster on large universes at high permutation counts.
Reproducibility
`gsea.py` is deterministic with a non-negative `--seed`: each contrast derives
its seed from the base seed and a hash of the contrast label, so results are
identical whether a contrast runs alone or alongside others. `TOP.py` has no
randomness. Figures are 300 dpi with editable vector text
(`svg.fonttype = "none"`, `pdf.fonttype = 42`).
Citation
If you use this code, please cite [paper] along with the underlying methods and
resources:
> Subramanian A, et al. Gene set enrichment analysis: a knowledge-based
> approach for interpreting genome-wide expression profiles. *PNAS*.
> 2005;102(43):15545–15550.
> Milacic M, et al. The Reactome Pathway Knowledgebase. *Nucleic Acids Res*.
> 2024;52(D1):D672–D678.
> Ashburner M, et al. Gene Ontology: tool for the unification of biology.
> *Nat Genet*. 2000;25(1):25–29.
