Pre-ranked GSEA for quantitative proteomics
A standalone, dependency-light implementation of pre-ranked Gene Set Enrichment
Analysis (GSEA) for quantitative proteomics data, with Reactome pathway support,
publication-ready figures, and optional restriction of plots to a chosen branch
of the Reactome hierarchy.
The script takes a differential-abundance table with one or more contrasts,
ranks proteins within each contrast, tests every pathway that meets a
size threshold, and writes results tables plus NES barplots, enrichment curves,
and a cross-contrast NES heatmap.
Contents
Installation
Quick start
Input files
Output files
Command-line options
Method
Choosing a significance metric
Reproducibility
Customisation
Limitations
Citation
Installation
Requires Python 3.9 or newer.
```bash
git clone https://github.com/<user>/<repo>.git
cd <repo>
pip install -r requirements.txt
```
Dependencies:
Package	Required	Purpose
numpy	yes	numerics
pandas	yes	tables
matplotlib	yes	all plots
seaborn	no	heatmap only; skipped with a warning if absent
tqdm	no	permutation progress bars
A minimal `requirements.txt`:
```
numpy>=1.24
pandas>=2.0
matplotlib>=3.7
seaborn>=0.13
tqdm>=4.66
```
Quick start
```bash
python gsea.py \
  --reg-table regulation_table_Amygdala.tsv \
  --pathway-map reactome_map_mouse.tsv \
  --id-col display_name \
  --stat-col log2FC \
  --outdir GSEA_out_Amygdala \
  --nperm 1000 \
  --min-size 10 \
  --topn 30 \
  --seed 1
```
Restricting the figures to two Reactome branches, while still testing and
correcting across the full pathway universe:
```bash
python gsea4.py \
  --reg-table regulation_table_Amygdala.tsv \
  --pathway-map reactome_map_mouse.tsv \
  --id-col display_name --stat-col log2FC \
  --comparisons "Alprazolam_KD vs Alprazolam, Alprazolam vs Control" \
  --outdir GSEA_out_Amygdala \
  --nperm 1000 --topn 30 --seed 1 \
  --hierarchy reactome_hierarchy_mouse.tsv \
  --ancestors "R-MMU-112316,R-MMU-1428517"
```
Input files
Regulation table (`--reg-table`)
Long-format TSV or CSV, one row per protein per contrast. Column names are
configurable; the defaults are shown in brackets.
Column	Flag	Default name	Description
contrast label	`--comparison-col`	`comparison`	groups rows into contrasts
protein/gene ID	`--id-col`	`gene`	must match the pathway map's member IDs
ranking statistic	`--stat-col`	`T-statistics`	e.g. log2FC or a moderated t-statistic
```
comparison	display_name	log2FC
Alprazolam vs Control	CALM1	-3.197
Alprazolam vs Control	NDUFA4	-1.093
```
IDs are uppercased before matching, so case differences between the two input
files are tolerated. Proteins appearing more than once within a contrast are
collapsed by averaging. Rows with a missing ID or statistic are dropped.
Pathway map (`--pathway-map`)
Long-format TSV or CSV, one row per pathway–member pair. By default the first
column is the pathway name and the second the member ID; override with
`--pathway-col` and `--member-col`. An optional `reactome_id` column is
required only if you intend to use `--ancestors`.
```
display_name	member	reactome_id
Respiratory electron transport	NDUFA4	R-MMU-611105
Respiratory electron transport	SDHB	R-MMU-611105
```
Reactome hierarchy (`--hierarchy`, optional)
Needed only for `--ancestors`. Requires `reactome_id` and `ancestor_ids`
(semicolon-separated); a `display_name` column additionally lets you pass
ancestors by name rather than ID.
```
reactome_id	display_name	ancestor_ids
R-MMU-611105	Respiratory electron transport	R-MMU-1428517;R-MMU-1430728
```
Output files
File	Contents
`gsea_results.csv`	every pathway tested, in every contrast
`gsea_leading_edge.csv`	leading-edge members per pathway and contrast
`gsea_results_filtered.csv`	the subset used for figures (identical to the above when no hierarchy filter is applied)
`gsea_nes_barplot__<contrast>.png`	top pathways by NES for one contrast
`gsea_enrichment__<pathway>__<contrast>.png`	running-sum enrichment curve with a hit rug
`gsea_nes_heatmap.png`	NES across all contrasts for the union of per-contrast top pathways
Columns in `gsea_results.csv`:
Column	Description
`pathway`	pathway name from the map
`reactome_id`	Reactome ID, if the map supplied one
`comparison`	contrast label
`ES`	enrichment score
`NES`	ES normalised by the mean of the same-signed null
`pval`	nominal permutation p-value, `(x + 1) / (nperm + 1)`
`BH_q`	Benjamini–Hochberg adjusted `pval`, within contrast
`FDR_q`	GSEA-style empirical FDR from the pooled null NES
`size`	pathway members detected in the ranked list
`leading_edge`	semicolon-separated leading-edge members
Rows are sorted by `BH_q` ascending within each contrast, with ties broken by
nominal p, then by |NES| descending, then by pathway name. The tie-break matters:
permutation p-values are floored at `1 / (nperm + 1)`, so the strongest pathways
frequently share an identical adjusted value, and |NES| decides their order.
Command-line options
Required
Flag	Description
`--reg-table`	differential-abundance table
`--pathway-map`	pathway-to-member map
`--outdir`	output directory, created if absent
Columns and contrast selection
Flag	Default	Description
`--comparison-col`	`comparison`	contrast column in the reg table
`--id-col`	`gene`	protein/gene ID column
`--stat-col`	`T-statistics`	ranking statistic column
`--comparisons`	all	comma-separated subset of contrasts to run
`--pathway-col`	first column	pathway column in the map
`--member-col`	second column	member column in the map
GSEA parameters
Flag	Default	Description
`--nperm`	1000	permutations per contrast
`--min-size`	10	minimum members detected in the ranked list
`--max-size`	500	maximum members detected in the ranked list
`--topn`	30	pathways per contrast carried into figures
`--seed`	1	base seed; `-1` for non-deterministic
Ranking and figures
Flag	Default	Description
`--sort-metric`	`BH_q`	metric ranking tables, barplots and heatmap rows (`BH_q`, `FDR_q`, `pval`)
`--heatmap-annot`	`bh`	cell annotation: `bh`, `fdr`, `pval`, `stars`, `none`
`--sig-threshold`	0.05	annotations are printed only at or below this value
`--p-fmt`	`.2g`	number format for annotations
`--heatmap-cutoff`	none	drop pathways from the heatmap above this value
`--heatmap-cutoff-metric`	`BH_q`	metric the cutoff applies to
Filtering and misc
Flag	Description
`--ancestors`	comma-separated Reactome ancestor IDs or names
`--hierarchy`	path to the hierarchy file; required by `--ancestors`
`--no-progress`	disable progress bars
Method
For each contrast independently:
Proteins are ranked in descending order by the chosen statistic. Duplicate
IDs are averaged, and ties are broken deterministically (by |statistic|, then
ID) so that the ordering is reproducible across runs and machines.
Pathways are retained when the number of members detected in that ranked
list falls within `--min-size` and `--max-size`. Note this is detected
size, not annotated pathway size — a distinction that matters in proteomics,
where coverage is well below transcriptome scale.
An enrichment score is computed from a running sum weighted by the absolute
ranking statistic (weighting exponent p = 1), taken as the maximum absolute
deviation from zero. Its sign indicates enrichment among up- or
downregulated proteins.
A null distribution is built by permuting protein labels across the ranked
list `--nperm` times with pathway membership fixed. Nominal p-values use the
`(x + 1) / (nperm + 1)` estimator against the same-signed tail.
ES values are normalised by the mean of the same-signed null distribution for
that pathway, making scores comparable across pathways of different size.
Nominal p-values are adjusted across all pathways tested in the contrast by
Benjamini–Hochberg (`BH_q`). A GSEA-style empirical FDR (`FDR_q`) is also
reported for reference.
Leading-edge members are those contributing to the running sum up to
(positive ES) or from (negative ES) the point of maximum deviation.
Hierarchy filtering via `--ancestors` is applied after testing and
multiple-testing correction. It affects only which pathways are drawn; it never
changes p-values, q-values, or the size of the tested pathway universe. This is
deliberate — filtering before correction would inflate significance.
Choosing a significance metric
Two adjusted values are reported and either can drive ranking via
`--sort-metric`.
`BH_q` is the default and the one to report. Benjamini–Hochberg controls FDR
under independence or positive regression dependence, and is conservative rather
than anti-conservative with discrete p-values such as permutation p-values.
`FDR_q` is the empirical GSEA-style estimate, computed as a ratio of null-tail
to observed-tail probabilities using the null NES pooled across pathways. It is
a heuristic without a formal guarantee, tends to run anti-conservative relative
to BH, and is not naturally monotone (the implementation enforces monotonicity).
Its one practical advantage is resolution: because it pools
`nperm × n_pathways` null values, it can still discriminate among pathways whose
nominal p-values have all hit the `1 / (nperm + 1)` floor. That is a symptom of
too few permutations rather than a virtue of the estimator — raising `--nperm`
is the better fix.
Reproducibility
With `--seed` set to a non-negative integer, results are deterministic. Each
contrast derives its own seed from the base seed and a SHA-256 hash of the
contrast label, so a contrast produces identical output whether run alone or
alongside others. Ranking ties are broken deterministically rather than by
input order. Pass `--seed -1` for non-deterministic behaviour.
Customisation
Contrast labels in figures are abbreviated through the `COMPARISON_LABELS`
dictionary near the top of the script. Labels not listed there are printed
verbatim, so edit this mapping for your own contrast names:
```python
COMPARISON_LABELS = {
    "Alprazolam vs Control": "Alp+CD vs CD",
    "Alprazolam_KD vs Alprazolam": "Alp+KD vs Alp+CD",
}
```
Figure fonts and heatmap cell geometry are set by module-level constants
(`HEATMAP_ANNOT_FONTSIZE`, `CELL_W_IN`, `CELL_H_IN`, and related). Output is
written at 300 dpi with editable text in vector formats.
Limitations
Gene-label permutation. The null is generated by permuting protein labels,
not sample labels, and so assumes proteins are independent. They are not:
complex members and pathway co-members co-vary strongly, and proteomic
quantification adds further correlation. Both `BH_q` and `FDR_q` are therefore
anti-conservative in absolute terms, and are best read as a ranking device
rather than as calibrated error rates. This limitation is shared with
GSEAPreranked and `fgsea` when run on a pre-ranked list. If the underlying
intensity matrix and sample labels are available, a sample-permutation null is
better calibrated.
Permutation floor. Nominal p-values cannot fall below `1 / (nperm + 1)`.
At the default 1000 permutations the floor is roughly 1e-3, so the top of the
list is often tied after correction.
Performance. The permutation loop is O(nperm × n_pathways) with a
Python-level inner loop. Large permutation counts on large pathway universes are
slow; `fgsea`, whose adaptive multilevel scheme estimates small p-values without
brute force, is the faster option at that scale.
