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

## 5. Known open item

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
a manuscript. See the chat write-up for the full comparison of the two
p-value sources.
