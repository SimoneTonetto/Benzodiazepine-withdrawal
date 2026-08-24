Proteomics pathway and reversal analysis
Two standalone Python scripts for downstream analysis of quantitative
proteomics differential-abundance tables:
Script	Purpose
`gsea.py`	pre-ranked GSEA against Reactome, with permutation testing and publication-ready figures
`TOP.py`	identify and plot individual proteins whose direction of change is reversed by a knockdown
Both read the same long-format regulation table and share a font/styling block,
so their figures can sit side by side in a panel without restyling.
Contents
Installation
The regulation table
`gsea.py` — pre-ranked GSEA
`TOP.py` — reverse-regulated proteins
Pathway map files
Reproducibility
Limitations
Citation
Installation
Requires Python 3.9 or newer.
```bash
git clone https://github.com/<user>/<repo>.git
cd <repo>
pip install -r requirements.txt
```
Package	Required by	Notes
numpy	both	
pandas	both	
matplotlib	both	
seaborn	`TOP.py`; `gsea.py` heatmap	`gsea.py` skips its heatmap with a warning if absent
tqdm	`gsea.py` only	optional; permutation progress bars
The regulation table
Both scripts consume a long-format TSV with one row per protein per contrast.
`gsea.py` needs only three columns and lets you name them; `TOP.py` needs the
full set below, because it also reconstructs per-condition abundances.
Column	Used by	Description
`comparison`	both	contrast label, e.g. `Alprazolam vs Control`
`display_name`	both	gene symbol
`identifier`	`TOP.py`	`SYMBOL~UNIPROT`; the UniProt part becomes the unique key
`log2FC`	both	log2 fold change
`padj2`	`TOP.py`	adjusted p-value for that contrast
`group1`, `group2`	`TOP.py`	condition names behind the contrast
`mean(group1)`, `mean(group2)`	`TOP.py`	group mean abundances
```
identifier	display_name	comparison	log2FC	padj2	group1	group2	mean(group1)	mean(group2)
Calm1~P0DP26	Calm1	Alprazolam vs Control	-1.84	0.0031	Alprazolam	Control	21403.5	76812.1
```
`TOP.py` validates these up front and names any that are missing, rather than
failing later inside the plotting code.
---
`gsea.py` — pre-ranked GSEA
Ranks proteins within each contrast, tests every pathway meeting a size
threshold by permutation, and writes results tables plus NES barplots,
enrichment curves and a cross-contrast heatmap.
Quick start
```bash
python gsea.py \
  --reg-table regulation_table_Amygdala.tsv \
  --pathway-map reactome_map_mouse.tsv \
  --id-col display_name \
  --stat-col log2FC \
  --outdir GSEA_out_Amygdala \
  --nperm 1000 --min-size 10 --topn 30 --seed 1
```
Restricting the figures to two Reactome branches, while still testing and
correcting across the full pathway universe:
```bash
python gsea.py \
  --reg-table regulation_table_Amygdala.tsv \
  --pathway-map reactome_map_mouse.tsv \
  --id-col display_name --stat-col log2FC \
  --comparisons "Alprazolam_KD vs Alprazolam, Alprazolam vs Control" \
  --outdir GSEA_out_Amygdala \
  --nperm 1000 --topn 30 --seed 1 \
  --hierarchy reactome_hierarchy_mouse.tsv \
  --ancestors "R-MMU-112316,R-MMU-1428517"
```
Outputs
File	Contents
`gsea_results.csv`	every pathway tested, in every contrast
`gsea_leading_edge.csv`	leading-edge members per pathway and contrast
`gsea_results_filtered.csv`	the subset used for figures
`gsea_nes_barplot__<contrast>.png`	top pathways by NES for one contrast
`gsea_enrichment__<pathway>__<contrast>.png`	running-sum enrichment curve with a hit rug
`gsea_nes_heatmap.png`	NES across contrasts for the union of per-contrast top pathways
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
Rows are sorted by `BH_q` ascending within each contrast, ties broken by
nominal p, then |NES| descending, then pathway name. The tie-break matters:
permutation p-values are floored at `1 / (nperm + 1)`, so the strongest
pathways frequently share an identical adjusted value, and |NES| decides
their order.
How many rows the heatmap has
`--topn` selects the most significant pathways within each contrast, and
the heatmap plots the union of those selections. With four contrasts,
`--topn 10` can therefore produce up to 40 rows — the union is only small when
the contrasts agree on which pathways matter.
To fix the total, use `--heatmap-topn`. Rows are first ranked by each pathway's
best (smallest) value of `--sort-metric` across all contrasts, then trimmed:
```bash
# Barplots and enrichment curves keep 30 per contrast;
# the heatmap shows exactly the 10 most significant pathways overall
python gsea.py ... --topn 30 --heatmap-topn 10
```
The run prints what was trimmed (`Heatmap limited to top 10 of 40 pathways`),
and the figure title states the count. Without `--heatmap-topn` the previous
union behaviour is unchanged. Note that `--heatmap-topn` caps rows after
selection, so a pathway ranked first in one contrast but mediocre elsewhere
still survives the cut — the ranking uses each pathway's best contrast, not its
average.
Options
Required
Flag	Description
`--reg-table`	regulation table
`--pathway-map`	pathway-to-member map
`--outdir`	output directory, created if absent
Columns and contrast selection
Flag	Default	Description
`--comparison-col`	`comparison`	contrast column
`--id-col`	`gene`	protein/gene ID column
`--stat-col`	`T-statistics`	ranking statistic column
`--comparisons`	all	comma-separated subset of contrasts
`--pathway-col`	first column	pathway column in the map
`--member-col`	second column	member column in the map
GSEA parameters
Flag	Default	Description
`--nperm`	1000	permutations per contrast
`--min-size`	10	minimum members detected in the ranked list
`--max-size`	500	maximum members detected in the ranked list
`--topn`	30	pathways per contrast carried into figures
`--heatmap-topn`	none	maximum total rows in the heatmap (see note below)
`--seed`	1	base seed; `-1` for non-deterministic
Ranking and figures
Flag	Default	Description
`--sort-metric`	`BH_q`	ranks tables, barplots and heatmap rows (`BH_q`, `FDR_q`, `pval`)
`--heatmap-annot`	`bh`	cell annotation: `bh`, `fdr`, `pval`, `stars`, `none`
`--sig-threshold`	0.05	annotate only at or below this value
`--p-fmt`	`.2g`	number format for annotations
`--heatmap-cutoff`	none	drop pathways above this value
`--heatmap-cutoff-metric`	`BH_q`	metric the cutoff applies to
Filtering and misc
Flag	Description
`--ancestors`	comma-separated Reactome ancestor IDs or names
`--hierarchy`	hierarchy file; required by `--ancestors`
`--no-progress`	disable progress bars
Method
For each contrast independently:
Proteins are ranked in descending order by the chosen statistic. Duplicate
IDs are averaged and ties are broken deterministically (by |statistic|, then
ID) so the ordering is reproducible across runs and machines.
Pathways are retained when the number of members detected in that ranked
list falls within `--min-size` and `--max-size`. This is detected size,
not annotated pathway size — a distinction that matters in proteomics,
where coverage is well below transcriptome scale.
An enrichment score is computed from a running sum weighted by the absolute
ranking statistic (weighting exponent p = 1), taken as the maximum absolute
deviation from zero. Its sign indicates enrichment among up- or
downregulated proteins.
A null distribution is built by permuting protein labels across the ranked
list `--nperm` times with pathway membership fixed. Nominal p-values use the
`(x + 1) / (nperm + 1)` estimator against the same-signed tail.
ES values are normalised by the mean of the same-signed null distribution
for that pathway, making scores comparable across pathways of different size.
Nominal p-values are adjusted across all pathways tested in the contrast by
Benjamini–Hochberg (`BH_q`). A GSEA-style empirical FDR (`FDR_q`) is also
reported for reference.
Leading-edge members are those contributing to the running sum up to
(positive ES) or from (negative ES) the point of maximum deviation.
Hierarchy filtering via `--ancestors` is applied after testing and
multiple-testing correction. It affects only which pathways are drawn; it never
changes p-values, q-values, or the size of the tested pathway universe.
Filtering before correction would inflate significance.
Choosing a significance metric
`BH_q` is the default and the one to report. Benjamini–Hochberg controls FDR
under independence or positive regression dependence, and is conservative
rather than anti-conservative with discrete p-values such as permutation
p-values.
`FDR_q` is the empirical GSEA-style estimate, a ratio of null-tail to
observed-tail probabilities using null NES pooled across pathways. It is a
heuristic without a formal guarantee, tends to run anti-conservative relative
to BH, and is not naturally monotone (the implementation enforces
monotonicity). Its one practical advantage is resolution: pooling
`nperm × n_pathways` null values lets it discriminate among pathways whose
nominal p-values have all hit the `1 / (nperm + 1)` floor. That is a symptom of
too few permutations rather than a virtue of the estimator — raising `--nperm`
is the better fix.
---
`TOP.py` — reverse-regulated proteins
Finds individual proteins that alprazolam changes significantly and the
knockdown pushes back the other way, then plots the strongest with per-protein
pathway annotation.
A protein is kept when it is significant in both contrasts (`padj2 < ALPHA`
in each) and its log2 fold changes have opposite signs. Survivors are
ranked by effect size — `padj2` is a gate here, not a sort key — and the top
`TOP_N` are plotted.
Quick start
```bash
# Use the settings in the USER PARAMETERS block at the top of the file
python TOP.py

# Or override on the command line
python TOP.py \
  --input regulation_table_Amygdala.tsv \
  --output-prefix alprazolam_switch_Amygdala \
  --top-n 50 --alpha 0.05
```
Every command-line default is the corresponding constant in the USER
PARAMETERS block, so running with no arguments reproduces the in-file
configuration exactly. Edit the block for a persistent setup; use flags for
one-off runs across brain regions.
Outputs
All prefixed with `--output-prefix` (default `alprazolam_direction_switch`):
File	Contents
`_top<N>_reverse_regulated.tsv`	the plotted subset, with annotation and source tag
`_all_reverse_regulated.tsv`	every reverse-regulated protein, not just the top N
`_heatmap_2comparisons.png`	log2FC | genes | z-score | pathway, the two reversal contrasts
`_heatmap_4comparisons.png`	as above, all four contrasts
`_heatmap_2comparisons_noZ.png`	log2FC | genes | pathway, no z-score panel
The count in the top-N filename is how many proteins were actually kept, so
`--top-n 0` yields `_top108_...` rather than `_top0_...`.
The log2FC panel prints `padj2` inside each cell; an exact zero is shown as
`< 0.0001`. Abundances are z-scored per protein across the four conditions, so
the z-score panel shows the shape of each protein's response, not its
magnitude.
Options
Flag	Default	Description
`--input`	`regulation_table_NAc.tsv`	regulation table
`--annotation`	`reactome_map_mouse.tsv`	Reactome gene-to-pathway map
`--go-fallback`	`go_map_mouse.tsv`	GO-BP map; `none` disables
`--output-prefix`	`alprazolam_direction_switch`	prefix for all outputs
`--alpha`	0.05	`padj2` threshold, applied to both contrasts
`--top-n`	50	proteins to plot; `0` keeps all
`--rank-by`	`reversal`	`reversal` = |log2FC_alp| + |log2FC_kd|; `alp` = |log2FC_alp| only
`--annotation-key` / `--annotation-col`	`gene` / `pathway`	columns in the Reactome map
`--go-key` / `--go-col`	`gene` / `pathway`	columns in the GO map
Annotation
Each protein gets one displayed pathway term. Reactome is primary and tagged
`[R]`; the GO-BP map is consulted only where Reactome is silent and tagged
`[GO]`, so the tag shows the specificity of each label at a glance. The run
prints how many proteins came from each source.
Where a protein maps to several terms, the one displayed is the first
alphabetically, not the most specific or most significant. This is a display
convenience for a one-line-per-protein figure; the complete term list is kept
in the output TSVs. If a particular figure needs a more meaningful term, pick
it manually rather than relying on this default.
Figure tuning
Layout constants sit near the top of the file: per-element font sizes
(`GENE_LABEL_FONTSIZE`, `PATHWAY_LABEL_FONTSIZE`, `HEATMAP_ANNOT_FONTSIZE` and
others), and pathway wrapping (`PATHWAY_WRAP_WIDTH`, `PATHWAY_MAX_LINES`,
`PATHWAY_LINESPACING`). Because each protein occupies a fixed-height row,
pathway text is capped at `PATHWAY_MAX_LINES` and truncated with an ellipsis;
the `[R]`/`[GO]` tag is always kept attached to the last line.
Display names are set by `GROUP_RENAME` (conditions) and `COMPARISON_RENAME`
(contrasts); `GROUP_ORDER` fixes the column order of the z-score panel. Names
not listed are printed verbatim, so edit these for your own study design.
---
Pathway map files
Both scripts read a gene-to-pathway map, but they expect different column
names by default:
Script	Pathway column	Member column	Configurable via
`gsea.py`	first column	second column	`--pathway-col`, `--member-col`
`TOP.py`	`pathway`	`gene`	`--annotation-col`, `--annotation-key`
A single file with columns `pathway`, `gene`, `reactome_id` serves both — pass
`--pathway-col pathway --member-col gene` to `gsea.py`, and `TOP.py` picks it
up with no flags. Keeping one map for both scripts is worth the small effort,
because it guarantees the pathway-level and protein-level figures are annotated
from the same source.
```
pathway	gene	reactome_id
Respiratory electron transport	Ndufa4	R-MMU-611105
Respiratory electron transport	Sdhb	R-MMU-611105
```
The `reactome_id` column is needed only for `gsea.py --ancestors`.
Reactome hierarchy (`gsea.py --ancestors` only)
Requires `reactome_id` and `ancestor_ids` (semicolon-separated); a
`display_name` column additionally lets you pass ancestors by name rather
than ID.
```
reactome_id	display_name	ancestor_ids
R-MMU-611105	Respiratory electron transport	R-MMU-1428517;R-MMU-1430728
```
Reproducibility
`gsea.py` is deterministic with `--seed` set to a non-negative integer. Each
contrast derives its own seed from the base seed and a SHA-256 hash of the
contrast label, so a contrast produces identical output whether run alone or
alongside others. Ranking ties are broken deterministically rather than by
input order. Pass `--seed -1` for non-deterministic behaviour.
`TOP.py` involves no randomness and is deterministic by construction.
Figures are written at 300 dpi with editable text preserved in vector formats
(`svg.fonttype = "none"`, `pdf.fonttype = 42`), so labels remain selectable in
Inkscape, Illustrator or PowerPoint.
Limitations
Gene-label permutation (`gsea.py`). The null permutes protein labels, not
sample labels, and so assumes proteins are independent. They are not: complex
members and pathway co-members co-vary strongly, and proteomic quantification
adds further correlation. Both `BH_q` and `FDR_q` are therefore
anti-conservative in absolute terms and are best read as a ranking device
rather than as calibrated error rates. This limitation is shared with
GSEAPreranked and `fgsea` when run on a pre-ranked list. If the underlying
intensity matrix and sample labels are available, a sample-permutation null is
better calibrated.
Permutation floor (`gsea.py`). Nominal p-values cannot fall below
`1 / (nperm + 1)`. At the default 1000 permutations the floor is roughly 1e-3,
so the top of the list is often tied after correction.
Performance (`gsea.py`). The permutation loop is O(nperm × n_pathways)
with a Python-level inner loop. Large permutation counts on large pathway
universes are slow; `fgsea`, whose adaptive multilevel scheme estimates small
p-values without brute force, is the faster option at that scale.
Ranking statistic (`gsea.py`). Ranking on log2 fold change ignores
measurement variance, so a large change in a noisy protein ranks alongside a
well-measured one. A moderated t-statistic is generally the better choice where
one is available.
Selection is not a formal interaction test (`TOP.py`). Requiring
significance in two contrasts and opposite signs is a conjunction of two
separate tests, not a test of the interaction between drug and knockdown. It
carries no calibrated error rate for the reversal itself, and the double
threshold is conservative in an uncontrolled way. Treat the output as a
prioritised shortlist for follow-up rather than a set of statistically
certified reversals; if a formal claim is needed, fit an interaction term in
the underlying model.
One term per protein (`TOP.py`). The displayed annotation is the first
alphabetically among possibly many, and mixes two ontologies of differing
specificity. The `[R]`/`[GO]` tags make the source visible, but the labels are
orientation for the reader, not a claim about the protein's primary function.
Citation
If you use `gsea.py`, please cite the original GSEA methodology:
> Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: a
> knowledge-based approach for interpreting genome-wide expression profiles.
> *Proc Natl Acad Sci USA*. 2005;102(43):15545–15550.
Where Reactome annotations are used:
> Milacic M, Beavers D, Conley P, et al. The Reactome Pathway Knowledgebase.
> *Nucleic Acids Res*. 2024;52(D1):D672–D678.
Where the GO Biological Process fallback is used:
> Ashburner M, Ball CA, Blake JA, et al. Gene Ontology: tool for the
> unification of biology. *Nat Genet*. 2000;25(1):25–29.
