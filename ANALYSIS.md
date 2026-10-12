# Corrected statistical workflow: methodology and findings

This document accompanies `metabolomics_pipeline.py` (the restructured
analysis script) and `results/` (the four figures + the exported
significant-metabolite tables). It records the methodological decisions
made when replacing the original from-scratch notebook, so the reasoning
is kept with the code rather than only in chat history.

## 1. Experimental design

3 groups (control, CMVS, CSDS) x 2 tissues (scW, BAT), n = 6 mice/group/tissue
(36 real injections; `BAT_Aq_CMVS6` is excluded as a pre-existing,
PCA-justified outlier, so BAT/CMVS runs on n = 5 vs n = 6). The two stress
protocols are never compared to each other — only each to its own control,
in each tissue. That is 4 pre-specified pairwise questions, not one
omnibus 3-group question.

## 2. Welch's t-test vs ANOVA

"There are three groups" and "the question is a specific pairwise
comparison" are different facts and do not imply the same statistics.
ANOVA (or Welch's ANOVA) answers "does group membership explain any
variance in this metabolite at all", which is not what was asked here —
CMVS vs CSDS is never a question of interest, and an omnibus test would
spend power on that irrelevant contrast and would still need post-hoc
tests to get back to the two answers actually wanted. Welch's t-test
(unequal variance, since nothing here justifies assuming equal spread at
n = 5-6) run separately for each pre-specified pairwise contrast, with
multiple-testing correction across the metabolite panel, answers the
actual questions directly. The three-group design does not by itself
force ANOVA; needing to control error across an unplanned scan of *all*
group differences would.

## 3. What changed from the original script, and why

| Change | Reason |
|---|---|
| Log2-transform peak areas before every statistic | LC-MS peak areas are right-skewed with multiplicative noise; log2 stabilises variance and makes log2FC and the t-test test the same quantity. Previously both the t-test and PLS-DA ran on raw areas. |
| Benjamini-Hochberg FDR added, per (tissue, comparison) panel | 301 tests at raw p < 0.05 gives ~15 false positives by chance alone; there was no correction before. Correcting the 4 panels separately (not pooled) because they are 4 distinct biological questions. |
| Significance rule is FDR < 0.05, full stop | No VIP > 1 AND-gate. VIP is demonstrably confounded with raw abundance in this dataset (that confound is what the figure was built to show) and is not stable at n = 5-6, so it is not used as a second hard significance filter. |
| PLS-DA label-permutation test added (1000 permutations) | The model had never been checked against permuted labels; with 301 variables and ~12 samples, high R2Y is easy to achieve on pure noise. |
| Dashed threshold line now marks the FDR boundary, not a fixed uncorrected p = 0.05 | Keeps the plot consistent with the stated significance rule. |
| Y-axis label corrected to "\|correlation with group\| (consistency)" | It always plotted \|r\|, not a p-value; the old label was simply wrong. |
| Panels with zero FDR survivors show an explicit "no metabolite reaches FDR < 0.05" panel | Silently drawing an empty scatter would look like a plotting bug rather than a real (and important) result. |

Kept unchanged: the NIPALS PLS-DA / VIP implementation (its self-consistency
check, VIP(1 comp) proportional to sqrt(SD)x|r|, still passes), the overall
figure layout (abundance x consistency, size+colour = VIP, dotted
equal-VIP curves, labelled archetypes), and the outlier exclusion.

## 4. Results (FDR < 0.05, Welch t-test on log2 peak area, BH-FDR per panel)

| Tissue | Comparison | n tested | raw p<0.05 | FDR<0.05 | PLS-DA R2Y | permutation p |
|---|---|---|---|---|---|---|
| scW | control vs CMVS | 301 | 136 | **75** | 0.91 | 0.112 |
| scW | control vs CSDS | 301 | 103 | **2**  | 0.93 | 0.064 |
| BAT | control vs CMVS | 301 | 36  | **3**  | 0.96 | 0.107 |
| BAT | control vs CSDS | 301 | 26  | **0**  | 0.95 | 0.240 |

All four PLS-DA models achieve high apparent R2Y (0.91-0.96) but **none**
clears the conventional permutation-test bar of p < 0.05 (closest: scW
CSDS, p = 0.064). This agrees with the workbook's own pre-existing
cross-validated Q2 values (`Cross_validation` sheet: BAT models Q2 as low
as 0.056-0.19; scW models better, Q2 0.64-0.80). Read VIP scores in every
panel as exploratory, not confirmatory, given this. Full tables:
`results/significant_metabolites_by_comparison.xlsx`.

## 5. Pathway / enrichment analysis: method

`Corrected_Stress_Adipose_Metabolomics_Analysis.ipynb` adds a pathway
analysis section. This records exactly what was and wasn't computed there,
since the enrichment *numbers* were not produced by this notebook.

### What actually computed the enrichment numbers

Not this notebook. The workbook's `Pathway_analysis` sheet already existed
before this restructuring, with output columns (`Total metabolites in the
pathway`, `Hits`, `p-value`, `-log(p)`, `Adjused p-value (Holm)`, `Adjused
p-value (FDR)`, `Pathway Impact`) that match MetaboAnalyst's pathway-analysis
module signature almost exactly, including its characteristic "Adjused"
typo. That's an inference from the column fingerprint, not a confirmed fact
— there's no metadata in the workbook stating which software produced it.
Standard MetaboAnalyst pathway analysis combines two independent
calculations per pathway:

* **Enrichment** (`Hits`, `p-value`, Holm/FDR-adjusted) — a hypergeometric
  or Fisher's-exact over-representation test: of the metabolites in this
  KEGG pathway, how many landed in the significant-hit list, versus what
  you'd expect by chance from the background of all tested/annotated
  compounds?
* **Topology** (`Pathway Impact`) — a separate calculation, typically
  relative-betweenness centrality of the hit metabolites' nodes within the
  pathway's KEGG network diagram, normalised to 0-1. A metabolite at a
  structurally central node contributes more impact than a peripheral one,
  independent of its p-value.

This ran against a KEGG pathway library (almost certainly the mouse `mmu`
library, given the tissue) that this session has no access to and cannot
verify.

### What hit list fed it

The pathway module needs a "significant metabolites" input list. That's the
workbook's own `Significant_metabolites` sheet, whose header states its
criterion explicitly: **"Adjusted p-values < 0.05 and/or VIP score > 1"** —
an OR rule. That is **not** the same list as this notebook's significance
rule (BH-FDR < 0.05 only, no VIP gate, Section 4 above). The pathway tables
therefore characterise enrichment among a broader, VIP-inclusive candidate
set, not among this notebook's FDR-confirmed hits.

### What this notebook actually did (extraction and plotting only)

1. **Located the block boundaries.** The sheet packs 3 pathway tables (one
   per comparison) stacked vertically with label rows in between, no
   structural separator: row 1 = `"1. scW: CMVS vs control"`, row 54 =
   `"2. scW: CSDS vs control"`, row 107 = `"3. BAT: CMVS vs control"`, sheet
   ends at row 160. `BAT: control vs CSDS` has no block at all.
2. **Parsed each block** (`load_pathway_block` in the notebook): header row
   is one row below the label, data runs from two rows below the label to
   the row before the next label (or end of sheet for the last block), then
   all non-name columns are coerced to numeric. Result: 51 pathway rows per
   comparison (the full KEGG library tested, not just the ones that came
   back significant).
3. **Mapped table columns to plot encodings** (`plot_pathway_enrichment`),
   one bubble per pathway:
   * x = `Pathway Impact` (topology score, as given)
   * y = `-log10(p-value)` (raw enrichment p, log-transformed for
     legibility)
   * marker size = `Hits`, clipped to a 20-400 pixel range
   * marker colour = `Adjused p-value (FDR)` on a fixed 0-1 scale, so
     colour is comparable across all three plots
   * dashed line at raw p = 0.05 for visual reference
   * the 5 smallest-p pathways get text labels, positioned with
     `adjustText` so overlapping labels separate with leader lines instead
     of rendering as mashed-together text (an earlier version of this plot
     had two labels collide on the scW/CSDS panel; fixed by switching to
     `adjustText`).
4. **No recomputation, no KEGG lookup.** This notebook did not build a
   compound-to-pathway membership map or run a hypergeometric test to
   produce these three plots. Outbound access to KEGG's REST API
   (`rest.kegg.jp`) was tested directly in this session, via both `curl`
   and the platform's own web-fetch tool, and confirmed blocked by the
   network egress proxy both times. Hand-writing an approximate
   compound-to-pathway table from memory to work around that was
   deliberately not done — a wrong pathway assignment would look exactly
   like a right one on the page, which is the kind of silent, unverifiable
   error this whole restructuring effort has been about removing.

The notebook does include one from-scratch, generic
`hypergeometric_ora(hit_ids, background_ids, pathway_to_members)` function,
self-tested against `scipy.stats.hypergeom` directly, but it sits unused in
this run because it has no `pathway_to_members` map to run against.
Supplying one (a KEGG export, or enabling network access to KEGG) is what
would let it run directly against this notebook's own FDR-significant hit
lists instead of reporting the workbook's pre-existing, differently-defined
results.

## 6. Known open item

The workbook also contains an already-normalised, already-FDR-corrected
MetaboAnalyst run (`Annotation`/`All_Data` sheets: `P-value`,
`Adj. P-value`, `VIP score` columns) computed on the full ~5,800-feature
untargeted dataset rather than just the 301 annotated compounds, and a
`Volume_normalization` sheet with per-sample tissue weight that this
script does not yet use (peak areas here are not corrected for tissue
weight/extraction volume). The two pipelines agree qualitatively — scW
CMVS has the strongest, most robust signal; BAT nets zero FDR-significant
hits in both pipelines — which is a reassuring cross-check, but the
weight-based normalisation should still be applied before this goes into
a manuscript.

**Update — the workbook's `P-value`/`Adj. P-value` test method is now
confirmed**, not inferred: "Hypothesis test performed by a one-way ANOVA
model with Tukey as post-hoc test. P-values are adjusted by the
Benjamini-Hochberg algorithm." Recomputing one-way ANOVA + Tukey HSD
(Control vs CMVS, Control vs CSDS contrasts) directly from this notebook's
raw peak areas and correlating against the workbook's values:

| Comparison | corr(Tukey p, workbook p) — raw scale | — log2 scale |
|---|---|---|
| scW/CMVS | 0.87 | 0.91 |
| scW/CSDS | 0.83 | 0.91 |
| BAT/CMVS | 0.77 | 0.83 |
| BAT/CSDS | 0.82 | 0.87 |

— better than assuming a plain two-group Welch t-test (~0.80-0.84
correlation), confirming the test family, but still not an exact
reproduction. The remaining gap is most likely the peak-area normalisation
step (tissue weight / extraction volume, `Volume_normalization` sheet)
that this notebook does not apply.

This does **not** change the significance rule used elsewhere in this
repository (Welch's t-test per planned pairwise contrast, Section 2-4
above): the actual questions here are two specific planned contrasts, not
a 3-group omnibus comparison, so Tukey HSD spends power on the
never-asked CMVS-vs-CSDS contrast and assumes equal variance across all
three groups — a stronger assumption than Welch needs at n = 5-6.
ANOVA+Tukey is a legitimate, common alternative; it isn't obviously the
better choice for this specific design. Knowing the exact method just
resolves what these workbook columns represent.

For reference, `results/workbook_anova_tukey_significant_metabolites.xlsx`
holds what "significant" means under the workbook's method alone
(`Adj. P-value < 0.05`, no VIP involved): 25 / 1 / 0 / 0 metabolites for
scW-CMVS / scW-CSDS / BAT-CMVS / BAT-CSDS respectively — more conservative
than this repository's own Welch-based list (75 / 2 / 3 / 0,
`results/all_significant_metabolites.xlsx`), but agreeing on the same
qualitative pattern: scW/CMVS carries by far the strongest signal, BAT has
next to none. Neither list is what feeds the `Pathway_analysis` sheet
(Section 5 above) — that one uses the workbook's own `Significant_metabolites`
sheet, whose criterion is "Adj. P-value < 0.05 **or** VIP > 1", a third,
looser rule.
