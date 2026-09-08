# The-Effects-of-Stress-in-Gene-Expression-in-Adipose-Tissues

Comparing the metabolite profiles of two fat depots in mice after two different chronic stress protocols.

---

## 1. What this project is, in plain words

Mice were put through stress. Two kinds of fat were then taken out of them and measured. The question is whether stress changes what is inside that fat.

The measuring was done by mass spectrometry, which reports a **peak area** for every metabolite it detects. A peak area is not a concentration. It is a signal size, so a bigger number means more of that molecule was there, but you cannot compare the number for one metabolite against the number for a different one. You can only compare the same metabolite between samples.

For each metabolite the analysis answers one question: **does it look different between the stressed mice and the control mice?**

That question is answered twice, in two different ways, on purpose:

1. **A t-test**, which looks at one metabolite at a time and gives a p-value.
2. **A PLS-DA model**, which looks at all 301 metabolites together and gives each one a VIP score.

The two answers do not agree, and the figure produced by this project is built to show exactly where and why they disagree. That is the main output.

### Why bother with two methods

The t-test is easy to explain but it ignores the fact that metabolites move together. If ten metabolites in the same pathway all shift a little, the t-test sees ten weak signals and may report none of them.

PLS-DA looks at the whole picture at once and can pick up that kind of coordinated shift. The cost is that its importance score, VIP, is harder to interpret. It has no p-value, no error control, and, as this project shows, it is partly driven by how abundant a metabolite is rather than by how reliably it changes.

Neither method is the correct one. They answer slightly different questions, and reporting only the one that gives the nicer answer would be cherry-picking. Showing both, and showing where they conflict, is the honest version.

---

## 2. The experiment

**Animals and stress protocols**

| Group | What it means |
|---|---|
| `control` | Unstressed mice |
| `CSDS` | Chronic social defeat stress |
| `CMVS` | Chronic multimodal variable stress |

**Tissues**

| Code | Tissue |
|---|---|
| `scW` | Subcutaneous white adipose tissue |
| `BAT` | Brown adipose tissue |

**Sample counts**: 6 mice per group per tissue, so 36 injections in total, plus 19 QC and blank injections that are skipped by the code.

One BAT sample, `BAT_Aq_CMVS6`, is excluded as an outlier. That decision was made from an earlier PCA and is hard-coded in the settings, which means BAT/CMVS runs on 11 samples instead of 12.

**Comparisons run**: four in total, each control against one stress group in one tissue.

- scW: control vs CSDS
- scW: control vs CMVS
- BAT: control vs CSDS
- BAT: control vs CMVS

The two stress groups are never compared against each other.

---

## 3. The data file

`Yana_copy_24April26_scW_BAT_metabolite_profiling_Data_ID.xlsx`, sheet **`Annotation`**.

This sheet holds the 301 metabolites that were actually identified, 141 with high confidence in the annotation and 160 with medium. The raw output of the instrument contains thousands more features that were never given a name, and none of those are used here.

Two details about how the sheet is read:

- The real column headers sit on the **third row**, so the code loads it with `header=2`.
- Every measurement column is named `Area: <tissue>_Aq_<group><number>.raw`, for example `Area: scW_Aq_CSDS3.raw (F12)`. The tissue and the group are read straight out of that name, so **the naming convention matters**. If injections are renamed, the loader stops recognising them.

Columns containing `Blank` or `QC` are skipped automatically.

---

## 4. How the analysis works, step by step

### Step 1: scaling

Peak areas range from around ten thousand to over a billion. Fed in raw, the biggest metabolites would dominate the model simply by being big.

**Pareto scaling** is used to reduce that. Each metabolite is centred on its own mean, then divided by the **square root** of its standard deviation. Dividing by the full standard deviation would make every metabolite equally important, which is the usual choice in other fields. The square root only shrinks the differences, so abundant metabolites still carry more weight than rare ones. This is the standard convention in metabolomics, and it is the direct reason the abundance effect shows up in the figures.

### Step 2: the model

A PLS-DA model is fitted from scratch using the NIPALS algorithm, with 2 components. In plain terms, the model searches for the combination of metabolites that best separates the control mice from the stressed mice, then repeats the search on whatever is left over.

Each metabolite then gets a **VIP score**, which measures how much it contributed to that separation. The common rule of thumb is that VIP above 1 is worth attention. It is only a convention. There is no p-value attached to it.

### Step 3: the univariate statistics

Separately, every metabolite gets a **Welch t-test**. Welch is used rather than the standard t-test because it does not assume the two groups have the same spread, which is safer with six animals per group.

The absolute correlation between a metabolite and the group label, written `|r|`, is also calculated. This is used as a measure of **consistency**: `|r|` near 1 means the metabolite separates the groups almost perfectly, and near 0 means it does not separate them at all.

### Step 4: a check that the model is behaving

The script prints a diagnostic each run:

```
VIP (1 comp) vs sqrt(SD)x|r| : r = 1.0000
```

With Pareto scaling and a single component, the VIP score should be exactly proportional to `sqrt(SD) x |r|`. If that correlation ever comes back as anything other than 1.0000, the model code is broken. It is a self-test, not a result.

---

## 5. The figure

One PNG per comparison, four in total.

**What is plotted**

- **X axis**: mean peak area on a log scale. Right means abundant.
- **Y axis**: `|correlation with group|`. Up means consistent.
- **Point size and colour**: both show VIP. Big and yellow means high VIP, small and purple means low.
- **Dashed line**: the p = 0.05 threshold, converted into the equivalent correlation so it can be drawn on this axis. For 12 samples it lands at `|r| = 0.576`. Only metabolites above it are inside the plotting window.
- **Dotted curves**: levels of equal VIP, at roughly 0.5, 1 and 3.
- **Labelled points**: four example metabolites, each with its category written underneath.

**The four labelled examples**

| Label | Category | What it demonstrates |
|---|---|---|
| A | high abundance + high VIP | The comfortable case. Both methods agree it matters. |
| B | low abundance + high consistency | Tracks the groups tightly, but gets only a modest VIP because it is small. |
| C | high abundance + high VIP, not significant | The model loves it, the t-test does not. |
| D | low VIP + significant | Tiny p-value, but the model rates it as unimportant. |

C and D are the point of the whole figure.

**What to take from it**

The dotted equal-VIP curves slope downwards from left to right. An abundant metabolite reaches VIP = 1 with only a modest correlation, while a rare one needs a very strong correlation to reach the same score. The big yellow points therefore cluster on the right side of the plot rather than at the top.

The printed diagnostics say the same thing in numbers. `|t|` against `sqrt(SD)` comes out near 0, so the t-test barely notices abundance. `|t|` against `|r|` comes out near 1, because the t-test and the correlation measure nearly the same thing.

**The conclusion**: a high VIP is not evidence that a metabolite responds to stress. It can be earned by being abundant. Any metabolite picked out by VIP alone should be checked against its p-value and its effect size before it goes into a results section.

---

## 6. Files in this project

| File | What it is |
|---|---|
| `Yana_copy_24April26_scW_BAT_metabolite_profiling_Data_ID.xlsx` | The data. Sheet `Annotation` is the one used. |
| `panel1_vip_map.py` | The full analysis and plotting script. Run it and it produces all four figures. |
| `Final_The_Effects_of_Stress_in_Gene_Expression_in_Adipose_Tissues_annotated.ipynb` | The same code as a notebook, with a plain-language note before every cell. Start here if you want to understand the method. |
| `figures/panel1_vip_vs_pvalue_<tissue>_<group>.png` | The four output figures. |

A note on the notebook filename. It says gene expression, but this is metabolite profiling and no transcripts are involved anywhere. The name should be corrected before this goes anywhere official.

---

## 7. Running it

**Requirements**: Python 3, plus `numpy`, `pandas`, `scipy`, `matplotlib` and `openpyxl` for reading the Excel file.

```bash
pip install numpy pandas scipy matplotlib openpyxl
```

**Run the script:**

```bash
METAB_XLSX="/path/to/Yana_copy_24April26_scW_BAT_metabolite_profiling_Data_ID.xlsx" \
FIG_OUT="figures" \
python3 panel1_vip_map.py
```

Both paths can also be edited directly in the settings block at the top of the script instead of being passed as environment variables.

**Expected output**: four PNG files in `FIG_OUT`, and a short diagnostic block per comparison printed to the terminal:

```
scW: control vs CSDS
  n = 12
  components = 2
  significance threshold |r| = 0.5760
  VIP (1 comp) vs sqrt(SD)x|r| : r = 1.0000
  |t| vs sqrt(SD) : r = 0.047
  |t| vs |r| : r = 0.970
```

**Settings worth knowing about**, all at the top of the script:

| Setting | Default | What it does |
|---|---|---|
| `SHEET` | `"Annotation"` | Which sheet to read |
| `N_COMPONENTS` | `2` | Number of PLS components |
| `EXCLUDED_SAMPLE` | `BAT_Aq_CMVS6` | Sample dropped from BAT |
| `ABUNDANT_TOP` | `60` | Rank cut-off for calling something abundant |
| `VIP_CURVE_LEVELS` | `[0.5, 1, 3]` | Which equal-VIP curves get drawn |
| `BOTTOM_MARGIN` | `0.04` | Space kept below the significance line so it stays visible |

---

## 8. Known issues

These are open, not fixed. Anyone using this analysis should read this section first.

1. **No multiple-testing correction.** With 301 metabolites tested at p < 0.05, roughly 15 hits are expected from chance alone. Until a Benjamini-Hochberg FDR column is added, no individual metabolite on these plots should be described as a finding.

2. **Small sample size.** Six animals per group, five in one case. Both p-values and VIP scores are unstable at that size.

3. **The model has not been validated.** No cross-validation (Q2) and no permutation test have been run, so there is currently no evidence that the separation the model finds would hold in new data. PLS-DA can separate two groups of pure noise if given enough variables, and 301 variables against 12 samples is exactly that situation.

4. **The y-axis label does not match the plotted data.** The label says p-value on a log scale. The values plotted are `|r|` on a linear scale. One of the two needs to change.

5. **Two different tests are mixed on one plot.** The dashed threshold line is derived from a Pearson correlation test while the p-values in the table come from a Welch t-test. They agree closely but not exactly, so a small number of points sit on the wrong side of the line.

6. **The excluded sample is undocumented.** `BAT_Aq_CMVS6` is dropped in the settings with a one-line comment. The PCA that justified dropping it needs to be written up and kept with the analysis.

7. **Peak areas are not normalised** for sample amount, tissue weight or instrument drift beyond whatever was done upstream in the processing software. If that normalisation did not happen, it should.

---

## 9. Next steps

- Add FDR-adjusted p-values and use them, rather than raw p, for the significance line.
- Cross-validate the PLS-DA model and run a permutation test before making any claim about group separation.
- Fix the axis label, and decide whether the y-axis should show `|r|` or the p-value.
- Add a volcano plot (log2 fold change against p-value) alongside the current figure, since that shows effect size, which neither VIP nor p-value does.
- Check whether the same metabolites come out on top in both tissues, and whether CSDS and CMVS push them in the same direction.
- Map the surviving hits onto pathways using the KEGG IDs already present in the data sheet.
