# Metabolomics Enrichment & Pathway Analysis

Reproducible analysis of the scW / BAT adipose-tissue metabolomics workbook
(Control vs CMVS vs CSDS), built to read every statistical result directly
from the source Excel file — nothing here is hand-typed.

## How to reproduce this

```bash
cd notebooks
jupyter nbconvert --to notebook --execute metabolomics_enrichment_pathway_analysis.ipynb
```

or open `notebooks/metabolomics_enrichment_pathway_analysis.ipynb` in
Jupyter/VS Code and run all cells top to bottom. It reads
`data/original_workbook.xlsx` and regenerates every file under `results/`
and `figures/`.

## Structure

```
metabolomics_enrichment_analysis/
├── notebooks/
│   └── metabolomics_enrichment_pathway_analysis.ipynb   <- the analysis
├── data/
│   └── original_workbook.xlsx                            <- source of truth (never modified)
├── results/
│   ├── significant_metabolites/   full + verification tables, per comparison
│   ├── enrichment/                 (see limitation below)
│   ├── pathway_analysis/           workbook's own pathway tables, extracted
│   ├── annotations/                identifier validation / annotation-quality tables
│   └── summary_tables/             QC, FDR cross-check, comparison summary, final report
├── figures/
│   ├── volcano/     Sig_pvalue_0.05 metabolites only, raw-p and adjusted-p versions, PNG+SVG
│   ├── enrichment/  bar charts, pathways significant at p<0.05 only (workbook's own pathway tables)
│   ├── pathway/     impact-vs-enrichment bubble plots, same significance filter
│   ├── heatmaps/    Sig_pvalue_0.05 metabolites, hierarchically clustered
│   ├── PCA/         one plot per comparison (stress vs its own control), restricted to that
│   │                comparison's Sig_pvalue_0.05 metabolites; a 1-metabolite comparison would get a
│   │                1-D strip plot instead of a degenerate 2-component PCA (not currently triggered)
│   ├── comparisons/ p-value distributions, MA plots (significant-only), cross-comparison
│   │                summary, UpSet-style overlap
│   └── top_metabolites/ top-25 (or fewer, if a comparison has fewer Sig_pvalue_0.05
│                     hits) metabolites per comparison, two views: a bar chart ranked
│                     by raw p-value (bar length = log2FC, coloured by direction), and
│                     a MetaboAnalyst-style VIP score dot plot (ranked by VIP score,
│                     with a Control/stress relative-abundance heat-strip)
└── README.md
```

## What this notebook does, and what it explicitly does not fabricate

- **Significance**: the workbook's own `Sig_pvalue_0.05` sheet (raw p<0.05,
  sourced from `Annotation`'s own `P-value` columns), per explicit
  instruction. Every one of the 167 metabolites it lists across the four
  comparisons was cross-checked against `Annotation` directly — same
  `Lab_ID`, exact p-value match, confirmed below 0.05, 0 mismatches (see
  `results/summary_tables/Sig_pvalue_0.05_verification.csv`). This raw-p
  criterion gates every figure and table in this notebook. The workbook's
  FDR-adjusted p-value (one-way ANOVA + Tukey HSD post-hoc,
  Benjamini-Hochberg FDR — confirmed method, not recomputed) is still
  computed and reported in every table for context, but is **not** used to
  decide what gets plotted; a secondary, clearly-labelled Benjamini-Hochberg
  recomputation from the raw p-value column is also shown as a cross-check
  only. This is a real methodology change from an earlier version of this
  notebook (which gated on FDR<0.05) and reintroduces the false-positive
  risk FDR correction exists to control — most visible for BAT, where 10
  `Sig_pvalue_0.05` metabolites per comparison do not survive FDR
  correction. See Section 13 and Section 31/32 of the notebook for the full
  discussion.
- **Pathway/enrichment analysis**: this environment has **no outbound
  network access to KEGG or MetaboAnalyst** (both confirmed blocked by the
  network egress policy). True over-representation analysis needs a
  compound-to-pathway membership map that isn't available here, and it was
  not fabricated. What's used instead: the workbook's own pre-computed
  `Pathway_analysis` sheet (3 of 4 comparisons; no `BAT: control vs CSDS`
  table exists in the workbook), reported as-is and clearly labelled as
  externally computed. **`Sig_pvalue_0.05` does not close this gap** — it
  gives BAT: control vs CSDS a real metabolite-level hit list (10
  metabolites) for the first time, but there is still no pathway-level
  table for that comparison in the workbook, so it still has no
  pathway/enrichment figure. A real, self-tested hypergeometric ORA
  function (`run_ora()`, in the notebook) is ready to run the moment a
  real pathway-membership map is supplied. Note also that the workbook's
  own `Pathway_analysis`/`Significant_metabolites` sheets use a third,
  separate significance criterion ("adjusted p<0.05 OR VIP>1"), distinct
  from both `Sig_pvalue_0.05` and this notebook's FDR column.
- **Every pathway/enrichment figure starts from the significance level**:
  only pathways clearing raw p < 0.05 are plotted (not the full 51-pathway
  library greyed out for context).
- **Numerator/denominator**: every ratio/log2FC column in the workbook is
  stress-group / control; "increased" always means higher in the stress
  group.

## Key results (see `results/summary_tables/FINAL_REPORT.md` for the full,
regenerated version)

| Comparison | Tested | **Sig_pvalue_0.05** (raw p<0.05, gates every figure) | FDR<0.05 (reported, not gating) | Increased | Decreased |
|---|---|---|---|---|---|
| scW: control vs CMVS | 301 | **82** | 25 | 82 | 0 |
| scW: control vs CSDS | 301 | **65** | 1 | 65 | 0 |
| BAT: control vs CMVS | 301 | **10** | 0 | 7 | 3 |
| BAT: control vs CSDS | 301 | **10** | 0 | 4 | 6 |

All 167 `Sig_pvalue_0.05` metabolites listed across the four comparisons
were cross-checked against `Annotation`'s own `P-value` columns by `Lab_ID`:
0 not found, 0 with a mismatched or non-significant p-value, max
`|sheet_p − Annotation_p|` = 0.0 (see
`results/summary_tables/Sig_pvalue_0.05_verification.csv`). Under FDR
correction, only scW's two comparisons survive; neither BAT comparison
does, despite now each having 10 raw-p-significant metabolites.

## Limitations

Sample size (n=5-6/group), multiple testing (301 simultaneous
tests/comparison; per explicit instruction, raw p<0.05 via `Sig_pvalue_0.05`
— not the FDR-adjusted column — is the criterion gating every figure here,
which reintroduces real false-positive risk FDR correction exists to
control, most visibly for BAT), metabolite annotation coverage (only
54/301, 17.9%, carry a KEGG ID), pathway database access (KEGG/MetaboAnalyst
network-blocked in this environment; `Sig_pvalue_0.05` fixes the
metabolite-level gap for BAT: control vs CSDS but not the separate
pathway-level gap, since no pathway table for that comparison exists in the
workbook), possible batch effects (workbook's `Volume_normalization` sheet
data not applied here), untargeted-metabolomics identification-confidence
limits, and association-vs-causation throughout. Full discussion in the
notebook, Section 32.
