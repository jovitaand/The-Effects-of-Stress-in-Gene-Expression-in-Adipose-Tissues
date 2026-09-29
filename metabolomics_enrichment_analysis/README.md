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
│   ├── volcano/     raw-p and adjusted-p versions, PNG+SVG, x4 comparisons
│   ├── enrichment/  bar charts, pathways significant at p<0.05 only
│   ├── pathway/     impact-vs-enrichment bubble plots, same significance filter
│   ├── heatmaps/    FDR-significant metabolites, hierarchically clustered
│   ├── PCA/         per-tissue, all 3 groups, log2+z-scored
│   └── comparisons/ p-value distributions, MA plots, cross-comparison summary, UpSet-style overlap
└── README.md
```

## What this notebook does, and what it explicitly does not fabricate

- **Significance**: workbook's own `Adj. P-value` columns (one-way ANOVA +
  Tukey HSD post-hoc, Benjamini-Hochberg FDR — confirmed method, not
  recomputed). Rule: FDR < 0.05, no VIP gate. A secondary, clearly-labelled
  Benjamini-Hochberg recomputation from the raw p-value column is shown as
  a cross-check only.
- **Pathway/enrichment analysis**: this environment has **no outbound
  network access to KEGG or MetaboAnalyst** (both confirmed blocked by the
  network egress policy). True over-representation analysis needs a
  compound-to-pathway membership map that isn't available here, and it was
  not fabricated. What's used instead: the workbook's own pre-computed
  `Pathway_analysis` sheet (3 of 4 comparisons; no `BAT: control vs CSDS`
  table exists in the workbook), reported as-is and clearly labelled as
  externally computed. A real, self-tested hypergeometric ORA function
  (`run_ora()`, in the notebook) is ready to run the moment a real
  pathway-membership map is supplied.
- **Every pathway/enrichment figure starts from the significance level**:
  only pathways clearing raw p < 0.05 are plotted (not the full 51-pathway
  library greyed out for context).
- **Numerator/denominator**: every ratio/log2FC column in the workbook is
  stress-group / control; "increased" always means higher in the stress
  group.

## Key results (see `results/summary_tables/FINAL_REPORT.md` for the full,
regenerated version)

| Comparison | Tested | Raw p<0.05 | FDR<0.05 | Increased | Decreased |
|---|---|---|---|---|---|
| scW: control vs CMVS | 301 | 82 | **25** | 25 | 0 |
| scW: control vs CSDS | 301 | 65 | **1** | 1 | 0 |
| BAT: control vs CMVS | 301 | 10 | **0** | 0 | 0 |
| BAT: control vs CSDS | 301 | 10 | **0** | 0 | 0 |

No metabolite is FDR-significant in more than one of the four comparisons.

## Limitations

Sample size (n=5-6/group), multiple testing (301 simultaneous tests/comparison,
FDR is the only significance criterion treated as confirmatory), metabolite
annotation coverage (only 54/301, 17.9%, carry a KEGG ID), pathway database
access (KEGG/MetaboAnalyst network-blocked in this environment), possible
batch effects (workbook's `Volume_normalization` sheet data not applied
here), untargeted-metabolomics identification-confidence limits, and
association-vs-causation throughout. Full discussion in the notebook,
Section 32.
