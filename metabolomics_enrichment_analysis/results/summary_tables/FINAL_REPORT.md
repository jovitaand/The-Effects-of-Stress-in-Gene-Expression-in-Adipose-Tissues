# Metabolomics Enrichment & Pathway Analysis -- Final Report

Generated from `original_workbook.xlsx`, sheet `Annotation` (301 metabolites).


## Statistical summary

| Comparison           |   Total tested |   Raw p < 0.05 |   FDR < 0.05 |   Increased |   Decreased |   KEGG mapped (of tested) |   Pathway mapped (of FDR-significant) |
|:---------------------|---------------:|---------------:|-------------:|------------:|------------:|--------------------------:|--------------------------------------:|
| scW: control vs CMVS |            301 |             82 |           25 |          25 |           0 |                        54 |                                     4 |
| scW: control vs CSDS |            301 |             65 |            1 |           1 |           0 |                        54 |                                     1 |
| BAT: control vs CMVS |            301 |             10 |            0 |           0 |           0 |                        54 |                                     0 |
| BAT: control vs CSDS |            301 |             10 |            0 |           0 |           0 |                        54 |                                     0 |


## Enrichment / pathway summary (workbook's own precomputed values)


**scW: control vs CMVS** -- 5 pathway(s) at raw p<0.05:

| Metabolic Pathway                           |   Hits |   Total metabolites in the pathway |    p-value |   Adjused p-value (FDR) |   Pathway Impact |
|:--------------------------------------------|-------:|-----------------------------------:|-----------:|------------------------:|-----------------:|
| Arginine biosynthesis                       |      7 |                                 14 | 2.3483e-06 |              0.00018786 |          0.50532 |
| Histidine metabolism                        |      7 |                                 16 | 7.1742e-06 |              0.00028697 |          0.67212 |
| Purine metabolism                           |     12 |                                 71 | 0.00022212 |              0.0053388  |          0.23343 |
| One carbon pool by folate                   |      7 |                                 26 | 0.00026694 |              0.0053388  |          0.20568 |
| Alanine, aspartate and glutamate metabolism |      7 |                                 28 | 0.00044049 |              0.0070478  |          0.66908 |

**scW: control vs CSDS** -- 5 pathway(s) at raw p<0.05:

| Metabolic Pathway                           |   Hits |   Total metabolites in the pathway |    p-value |   Adjused p-value (FDR) |   Pathway Impact |
|:--------------------------------------------|-------:|-----------------------------------:|-----------:|------------------------:|-----------------:|
| Purine metabolism                           |     16 |                                 71 | 6.8192e-07 |              5.4553e-05 |          0.35365 |
| Arginine biosynthesis                       |      7 |                                 14 | 3.5351e-06 |              0.0001414  |          0.50532 |
| Histidine metabolism                        |      7 |                                 16 | 1.0736e-05 |              0.00021671 |          0.67212 |
| Alanine, aspartate and glutamate metabolism |      9 |                                 28 | 1.0836e-05 |              0.00021671 |          0.66908 |
| Arginine and proline metabolism             |      9 |                                 36 | 0.00010151 |              0.0016242  |          0.39767 |

**BAT: control vs CMVS** -- 5 pathway(s) at raw p<0.05:

| Metabolic Pathway               |   Hits |   Total metabolites in the pathway |    p-value |   Adjused p-value (FDR) |   Pathway Impact |
|:--------------------------------|-------:|-----------------------------------:|-----------:|------------------------:|-----------------:|
| Purine metabolism               |     14 |                                 71 | 3.7821e-06 |              0.00030257 |          0.30107 |
| Arginine biosynthesis           |      6 |                                 14 | 2.5633e-05 |              0.0010253  |          0.44149 |
| Histidine metabolism            |      6 |                                 16 | 6.314e-05  |              0.0016837  |          0.54917 |
| One carbon pool by folate       |      7 |                                 26 | 0.00016416 |              0.0032831  |          0.2392  |
| Arginine and proline metabolism |      8 |                                 36 | 0.00023558 |              0.0037693  |          0.21279 |


**BAT: control vs CSDS** -- no pathway table exists in the workbook for this comparison.


## Tissue comparison

- scW is associated with a substantially larger, FDR-robust metabolite response to CMVS than BAT (25 vs 0 FDR-significant metabolites).

- Under CSDS, both tissues show a much smaller response (scW: 1, BAT: 0).

- No metabolite is FDR-significant in more than one of the four comparisons (confirmed by the UpSet-style overlap plot).


## Treatment comparison

- Within scW, CMVS produces a much larger significant response than CSDS (25 vs 1).

- Within BAT, neither stressor produces an FDR-significant metabolite at this threshold.


## Limitations

See Section 32 above for the full list (sample size, multiple testing, annotation coverage, pathway database access, possible batch effects, untargeted-metabolomics caveats, association vs causation).
