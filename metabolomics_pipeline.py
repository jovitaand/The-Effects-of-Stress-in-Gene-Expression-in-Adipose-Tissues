"""
Corrected metabolomics analysis and visualization pipeline.

Stress x adipose-tissue metabolomics: Control vs CMVS and Control vs CSDS,
in subcutaneous white fat (scW) and brown fat (BAT). Four comparisons total.

This replaces the earlier from-scratch notebook (see README section
"Known issues" / commit history). The scientific questions and the overall
figure design are unchanged. What changed, and why, is documented in the
section headers below and in METHODOLOGY_CHANGES.md.

Run:
    METAB_XLSX="/path/to/Yana copy_24April26_scW_BAT_metabolite_profiling_Data_ID.xlsx" \
    FIG_OUT="results" \
    python3 metabolomics_pipeline.py
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests

# ==========================================================================
# 1. SETTINGS
# ==========================================================================

EXCEL_PATH = os.environ.get(
    "METAB_XLSX",
    "Yana copy_24April26_scW_BAT_metabolite_profiling_Data_ID.xlsx",
)

OUTPUT_DIR = os.environ.get("FIG_OUT", "results")

SHEET = "Annotation"

N_COMPONENTS = 2

# Sample dropped from BAT. It sat away from everything else in an earlier
# PCA (see Cross_validation sheet: "All" vs "All no outlier" PCA-X models).
# This is a pre-existing decision carried over from the original analysis;
# it is not re-derived here. See Part 7 / Known issues for why this still
# needs a proper write-up with the PCA that justified it.
EXCLUDED_SAMPLE = {
    "scW": None,
    "BAT": "BAT_Aq_CMVS6",
}

# ---- Statistical significance rule -------------------------------------
# FDR (Benjamini-Hochberg) applied SEPARATELY within each of the four
# (tissue x comparison) panels. See Part 4 of the write-up for why this
# unit of correction was chosen over correcting across all four at once,
# or across tissues.
ALPHA_FDR = 0.05

# ---- PLS-DA validation ---------------------------------------------------
N_PERMUTATIONS = 1000
RNG_SEED = 0

# ---- Plot cosmetics (unchanged from the original script) -----------------
VIP_COLOUR_CAP = 4.0
VIP_CURVE_LEVELS = [0.5, 1, 3]
ABUNDANT_TOP = 60
BOTTOM_MARGIN = 0.04
X_PAD_DECADES_LEFT = 0.18
X_PAD_DECADES_RIGHT = 0.35

COMPARISONS = [
    ("scW", "CMVS"),
    ("scW", "CSDS"),
    ("BAT", "CMVS"),
    ("BAT", "CSDS"),
]


# ==========================================================================
# 2. LOAD DATA
# ==========================================================================

def describe_injection(column_name):
    """Return (tissue, group) for an injection column, or None for QC/blanks."""

    if "Blank" in column_name or "QC" in column_name:
        return None

    tissue = "scW" if column_name.startswith("Area: scW") else "BAT"

    if "CMVS" in column_name:
        group = "CMVS"
    elif "CSDS" in column_name:
        group = "CSDS"
    else:
        group = "control"

    return tissue, group


def load_table(excel_path=None, sheet=SHEET):
    """Load the 301-metabolite annotation sheet and its peak-area columns."""

    excel_path = excel_path or EXCEL_PATH

    table = pd.read_excel(excel_path, sheet_name=sheet, header=2)

    area_columns = [c for c in table.columns if str(c).startswith("Area: ")]

    intensities = table[area_columns].apply(pd.to_numeric, errors="coerce")

    return table, intensities


# ==========================================================================
# 3 & 4. DEFINE GROUPS / SEPARATE TISSUES
# ==========================================================================

def select_samples(intensities, tissue, case_group, control_group="control"):
    """
    Columns for one (tissue, case vs control) comparison, QC/blanks and the
    documented outlier removed.
    """

    excluded = EXCLUDED_SAMPLE.get(tissue)
    selected = []

    for column in intensities.columns:
        info = describe_injection(column)
        if info is None:
            continue
        if info[0] != tissue:
            continue
        if info[1] not in (control_group, case_group):
            continue
        if excluded and excluded in column:
            continue
        selected.append(column)

    labels = np.array([describe_injection(c)[1] for c in selected])

    return selected, labels


# ==========================================================================
# 5. PREPROCESS METABOLITE DATA
# ==========================================================================
#
# CHANGE from the original script: peak areas are log2-transformed before
# every downstream statistic (Welch t-test, correlation, PLS-DA). LC-MS
# peak areas are strongly right-skewed and multiplicative in their noise
# (a metabolite twice as abundant tends to have twice the absolute spread),
# which is exactly the situation log transformation is for: it symmetrises
# the distribution, stabilises the variance across metabolites of very
# different abundance, and turns "fold change" into a difference of means
# (so log2FC and the t-test are now testing the same quantity). The
# original script ran the t-test and the PLS-DA on raw areas, so a metabolite
# at 10^9 and one at 10^5 were assumed to have comparable noise structure,
# which is not realistic for this kind of data. This is the single largest
# methodological change in this rewrite; see Part 7 for the full comparison.
#
# Pareto scaling for the PLS-DA step is unchanged in spirit (still
# center-and-divide-by-sqrt(SD)) but now operates on log2 data.

def log2_transform(matrix):
    return np.log2(matrix)


def pareto_scale(matrix):
    """Centre each metabolite, then divide by the square root of its SD."""

    mean = matrix.mean(axis=0)
    sd = matrix.std(axis=0, ddof=1)
    sd = np.where(sd == 0, 1.0, sd)

    return (matrix - mean) / np.sqrt(sd)


def make_dummy_y(labels):
    """Class labels -> centred, unit-variance dummy response matrix."""

    classes = sorted(set(labels))
    response = np.zeros((len(labels), len(classes)))

    for row, label in enumerate(labels):
        response[row, classes.index(label)] = 1.0

    response = response - response.mean(axis=0)
    sd = response.std(axis=0, ddof=1)
    sd = np.where(sd == 0, 1.0, sd)

    return response / sd


# ==========================================================================
# 6. STATISTICAL TESTING (Welch pairwise t-test, planned contrasts)
# ==========================================================================
#
# Each tissue x stress-group panel is a separate, PRE-SPECIFIED pairwise
# question ("does CSDS shift this metabolite relative to control in scW?"),
# not an exploratory scan of a three-group omnibus design. See Part 2/3 of
# the write-up for the full justification of Welch's t-test over one-way
# ANOVA here. Welch's (unequal-variance) t-test is used throughout because
# with n=5-6 per group there is no reason to assume equal within-group
# variance, and Welch's test costs essentially nothing when variances
# happen to be similar.

def welch_test(matrix, is_case):
    return stats.ttest_ind(
        matrix[is_case], matrix[~is_case], axis=0, equal_var=False
    )


# ==========================================================================
# 7. MULTIPLE-TESTING CORRECTION
# ==========================================================================
#
# CHANGE from the original script: this is the fix for "Known issue #1" in
# the README. 301 metabolites x 4 comparisons with no correction means
# roughly 15 features per comparison are expected to cross p < 0.05 by
# chance alone even if stress changed nothing. Benjamini-Hochberg FDR is
# applied WITHIN each of the four (tissue, comparison) panels separately,
# not pooled across all four. See Part 4 for why: the four panels are four
# different biological questions (two tissues, two independent stressors),
# each with its own 301-metabolite test family; pooling would let a strong
# result in one tissue "borrow" significance for an unrelated test in the
# other tissue, and would make the correction depend on which comparisons
# happen to be run together, which is not a property multiple-testing
# correction should have.

def fdr_correct(p_values):
    _, adjusted, _, _ = multipletests(p_values, method="fdr_bh")
    return adjusted


# ==========================================================================
# 8. PLS-DA (NIPALS, unchanged algorithm; new validation wrapped around it)
# ==========================================================================

def extract_pls_component(x_block, y_block, tolerance=1e-12, max_iterations=2000):
    """Extract one PLS2 component using NIPALS. Unchanged from the original."""

    start = int(np.argmax((y_block ** 2).sum(axis=0)))
    y_score = y_block[:, [start]].copy()
    x_weight = x_score = y_loading = None

    for _ in range(max_iterations):
        x_weight = x_block.T @ y_score
        norm = np.linalg.norm(x_weight)
        if norm == 0:
            break
        x_weight = x_weight / norm
        x_score = x_block @ x_weight
        y_loading = (y_block.T @ x_score) / (x_score.T @ x_score)
        updated = (y_block @ y_loading) / (y_loading.T @ y_loading)
        change = np.linalg.norm(updated - y_score)
        scale = max(1.0, np.linalg.norm(y_score))
        y_score = updated
        if change < tolerance * scale:
            break

    x_loading = (x_block.T @ x_score) / (x_score.T @ x_score)

    return x_weight, x_score, x_loading, y_loading


def fit_pls(scaled_matrix, dummy_y, n_components):
    """Fit an n-component PLS2 model. Returns weights/scores/loadings and R2Y."""

    x_residual = scaled_matrix.copy()
    y_residual = dummy_y.copy()

    weights, scores, y_loadings = [], [], []
    ss_y_total = (dummy_y ** 2).sum()

    for _ in range(n_components):
        w, t, p, c = extract_pls_component(x_residual, y_residual)
        x_residual = x_residual - t @ p.T
        y_residual = y_residual - t @ c.T
        weights.append(w.ravel())
        scores.append(t.ravel())
        y_loadings.append(c.ravel())

    weights = np.array(weights).T
    scores = np.array(scores).T
    y_loadings = np.array(y_loadings).T

    r2y = 1.0 - (y_residual ** 2).sum() / ss_y_total

    return weights, scores, y_loadings, r2y


def vip_from_fit(scaled_matrix, weights, scores, y_loadings):
    """VIP scores from an already-fitted PLS model."""

    n_components = weights.shape[1]

    explained = np.array([
        (scores[:, i] @ scores[:, i]) * (y_loadings[:, i] @ y_loadings[:, i])
        for i in range(n_components)
    ])

    normalised = weights / np.linalg.norm(weights, axis=0, keepdims=True)
    contribution = (normalised ** 2) @ explained

    return np.sqrt(scaled_matrix.shape[1] * contribution / explained.sum())


def fit_vip(scaled_matrix, dummy_y, n_components):
    weights, scores, y_loadings, r2y = fit_pls(scaled_matrix, dummy_y, n_components)
    vip = vip_from_fit(scaled_matrix, weights, scores, y_loadings)
    return vip, r2y


# ==========================================================================
# 9. PLS-DA VALIDATION: a label-permutation test
# ==========================================================================
#
# NEW section. The original script fit one PLS-DA model per comparison and
# reported VIP with no check that the model's apparent group separation
# would hold up on unseen data. With 301 (correlated) variables and 11-12
# samples this is exactly the regime where PLS-DA can separate pure noise
# (see README "Known issue #3"). The check added here: shuffle the group
# labels many times, refit R2Y (the same model, same code path, just wrong
# labels) each time, and see how often permuted data reaches an R2Y as high
# as the real one. That gives an empirical p-value for "the model found
# real group structure, not just 301 variables' worth of chances to fit
# noise."
#
# This is read alongside the pre-existing MetaboAnalyst cross-validation
# already recorded in the workbook's Cross_validation sheet (Q2(cum) per
# comparison), which is the authoritative check: it was computed on all
# ~5,800 detected features, not just the 301 annotated ones, so it is more
# powerful than anything that can be re-derived here from the annotated
# subset alone. The permutation test below is a reproducible, from-the-
# same-data confirmation of that same conclusion, not a replacement for it
# (see Part 6 of the write-up for how the two agree).

def permutation_test_r2y(matrix, labels, n_components, n_permutations=N_PERMUTATIONS, seed=RNG_SEED):
    """Empirical p-value: how often a random relabelling reaches this R2Y."""

    scaled = pareto_scale(matrix)
    dummy = make_dummy_y(labels)
    _, _, _, observed_r2y = fit_pls(scaled, dummy, n_components)

    rng = np.random.default_rng(seed)
    permuted_r2y = np.empty(n_permutations)

    for k in range(n_permutations):
        shuffled = rng.permutation(labels)
        dummy_perm = make_dummy_y(shuffled)
        _, _, _, r2y = fit_pls(scaled, dummy_perm, n_components)
        permuted_r2y[k] = r2y

    p_value = (np.sum(permuted_r2y >= observed_r2y) + 1) / (n_permutations + 1)

    return observed_r2y, p_value, permuted_r2y


# ==========================================================================
# 10. RUN ONE COMPARISON: stats + VIP + significance
# ==========================================================================

def analyse(tissue, case_group, control_group="control",
            n_components=N_COMPONENTS, excel_path=None, run_permutation=True):
    """Compute every quantity needed for the tables and the figure."""

    table, intensities = load_table(excel_path)
    selected, labels = select_samples(intensities, tissue, case_group, control_group)

    raw_matrix = intensities[selected].values.T
    log_matrix = log2_transform(raw_matrix)
    is_case = labels == case_group

    # ---- Consistency (|r|) and VIP, computed on log2 data -----------------
    correlation = np.array([
        np.corrcoef(log_matrix[:, j], is_case.astype(float))[0, 1]
        for j in range(log_matrix.shape[1])
    ])

    scaled = pareto_scale(log_matrix)
    dummy = make_dummy_y(labels)
    vip, r2y = fit_vip(scaled, dummy, n_components)

    # ---- Univariate statistics: Welch t-test on log2 data ------------------
    t_stat, p_value = welch_test(log_matrix, is_case)
    adj_p_value = fdr_correct(p_value)

    log2fc = (
        log_matrix[is_case].mean(axis=0) - log_matrix[~is_case].mean(axis=0)
    )

    results = pd.DataFrame({
        "Tissue": tissue,
        "Comparison": f"{control_group} vs {case_group}",
        "Metabolite": table["Metabolite"].astype(str).values,
        "Confidence": table["Confidence in annotation"].astype(str).values,
        "area": raw_matrix.mean(axis=0),
        "SD_log2": log_matrix.std(axis=0, ddof=1),
        "abs_r": np.abs(correlation),
        "VIP": vip,
        "p_value": p_value,
        "adj_p_value": adj_p_value,
        "t_stat": np.abs(t_stat),
        "log2FC": log2fc,
        "Direction": np.where(log2fc > 0, "UP in " + case_group, "DOWN in " + case_group),
        "significant_FDR": adj_p_value < ALPHA_FDR,
    })

    results["area_rank"] = results["area"].rank(ascending=False, method="min").astype(int)
    results["vip_rank"] = results["VIP"].rank(ascending=False, method="min").astype(int)
    results["p_rank"] = results["p_value"].rank(method="min").astype(int)

    validation = {}
    if run_permutation:
        observed_r2y, perm_p, _ = permutation_test_r2y(log_matrix, labels, n_components)
        validation = {"R2Y": observed_r2y, "permutation_p": perm_p}

    return {
        "results": results,
        "tissue": tissue,
        "case_group": case_group,
        "control_group": control_group,
        "n_samples": len(selected),
        "n_components": n_components,
        "validation": validation,
    }


# ==========================================================================
# 11. IDENTIFY SIGNIFICANT METABOLITES
# ==========================================================================

def significant_metabolites(analysis):
    """
    FDR < ALPHA_FDR is the significance rule (see Part 4/7 write-up).
    VIP is reported alongside for context but is NOT used as an additional
    hard filter: in this dataset VIP is demonstrably driven in part by raw
    abundance rather than by group separation alone (that confound is the
    entire point of the figure below), and with n=5-6 per group VIP is not
    stable enough to trust as a second independent significance criterion.
    A metabolite is not excluded from the "significant" table for having a
    low VIP, and a high VIP is never treated as significance on its own.
    """

    results = analysis["results"]

    cols = ["Tissue", "Comparison", "Metabolite", "log2FC", "p_value",
            "adj_p_value", "VIP", "Direction"]

    return results.loc[results["significant_FDR"], cols].sort_values("adj_p_value")


# ==========================================================================
# 12. PREPARE VISUALIZATION DATA / AXIS LIMITS
# ==========================================================================

def significance_correlation(n_samples, alpha=0.05):
    """Absolute correlation corresponding to p = alpha, two-sided."""

    degrees_of_freedom = n_samples - 2
    t_critical = stats.t.ppf(1 - alpha / 2, degrees_of_freedom)

    return t_critical / np.sqrt(t_critical ** 2 + degrees_of_freedom)


def fdr_threshold_as_correlation(results, n_samples):
    """
    The dashed threshold line is now drawn at the |r| that corresponds to
    the LARGEST raw p-value among the FDR-significant metabolites in this
    panel (i.e. the data-dependent boundary FDR actually drew), rather than
    at a fixed p = 0.05. If nothing survives FDR, there is no line to draw:
    the caller must handle that case explicitly (see build_figure). This
    replaces the original fixed p < 0.05 threshold, which is no longer the
    stated significance rule.
    """

    sig = results.loc[results["significant_FDR"]]
    if len(sig) == 0:
        return None

    boundary_p = sig["p_value"].max()
    return significance_correlation(n_samples, alpha=boundary_p)


def compute_limits(results, threshold):
    """Work out the visible window before anything is drawn."""

    visible = results[
        (results["abs_r"] >= threshold)
        & (results["area"] > 0)
        & np.isfinite(results["area"])
    ]

    if len(visible) == 0:
        visible = results[results["area"] > 0]

    top = float(visible["abs_r"].max())
    span = max(top - threshold, 0.05)

    y_min = threshold - BOTTOM_MARGIN * span
    y_max = min(1.0, top + 0.08 * span)

    log_min = np.log10(visible["area"].min())
    log_max = np.log10(visible["area"].max())

    x_min = 10 ** (log_min - X_PAD_DECADES_LEFT)
    x_max = 10 ** (log_max + X_PAD_DECADES_RIGHT)

    return (x_min, x_max), (y_min, y_max)


def axes_fraction(axis, x, y):
    x_min, x_max = axis.get_xlim()
    y_min, y_max = axis.get_ylim()
    fx = (np.log10(x) - np.log10(x_min)) / (np.log10(x_max) - np.log10(x_min))
    fy = (y - y_min) / (y_max - y_min)
    return fx, fy


def pick_archetypes(results, abundant_top=ABUNDANT_TOP):
    """Representative metabolites for annotation, restricted to FDR-significant ones."""

    sig = results[results["significant_FDR"]]

    abundant = sig["area_rank"] <= abundant_top
    high_vip = sig["VIP"] > 1

    def best(subset, column, largest=True):
        if len(subset) == 0:
            return None
        picked = subset.nlargest(1, column) if largest else subset.nsmallest(1, column)
        return picked.iloc[0]

    return {
        "A": (best(sig[high_vip & abundant], "VIP"), "high abundance + high VIP"),
        "B": (best(sig[high_vip & ~abundant], "abs_r"), "low abundance + high consistency"),
        "D": (best(sig[~high_vip], "abs_r", largest=True), "low VIP + FDR-significant"),
    }


# ==========================================================================
# 13. DRAW THE PANEL / EXPORT
# ==========================================================================

def draw_equal_vip_curves(axis, results, levels=VIP_CURVE_LEVELS):
    usable = results[
        (results["area"] > 0) & (results["SD_log2"] > 0)
        & np.isfinite(results["VIP"]) & np.isfinite(results["abs_r"])
    ]
    if len(usable) < 10:
        return

    slope, intercept = np.polyfit(
        np.log10(usable["area"]), np.log10(np.sqrt(usable["SD_log2"])), 1
    )

    x_min, x_max = axis.get_xlim()
    y_min, y_max = axis.get_ylim()
    x_values = np.logspace(np.log10(x_min), np.log10(x_max), 400)
    predicted_sqrt_sd = 10 ** (intercept + slope * np.log10(x_values))

    product = np.sqrt(usable["SD_log2"].values) * usable["abs_r"].values
    good = (product > 0) & (usable["VIP"].values > 0)
    if not good.any():
        return

    constant = np.median(usable["VIP"].values[good] / product[good])
    label_gap = 0.015 * (y_max - y_min)

    for level in levels:
        curve = level / (constant * predicted_sqrt_sd)
        inside = (curve >= y_min) & (curve <= y_max)
        if inside.sum() < 3:
            continue

        axis.plot(x_values[inside], curve[inside], linestyle=":", linewidth=1.0,
                   color="black", alpha=0.6, zorder=1, clip_on=True)

        highest = int(np.argmax(np.where(inside, curve, -np.inf)))
        label_x = x_values[highest]
        label_y = min(curve[highest], y_max - label_gap)
        fx, _ = axes_fraction(axis, label_x, label_y)
        alignment = "left" if fx < 0.06 else ("right" if fx > 0.94 else "center")

        axis.text(label_x, label_y, f"VIP ≈ {level:g}", fontsize=7, alpha=0.85,
                   va="top", ha=alignment, zorder=5, clip_on=True,
                   bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1.0))


def panel_map(axis, analysis, archetypes):
    """
    Same layout as the original Panel 1: x = abundance (log scale),
    y = |correlation with group| ("consistency"), size + colour = VIP.

    CHANGES from the original:
      * y-axis label now says "|correlation with group|" instead of
        "p-value" -- it always plotted |r|, the label was simply wrong
        (README known issue #4).
      * the dashed threshold line marks the FDR < 0.05 boundary for this
        panel, not a fixed uncorrected p < 0.05.
      * only FDR-significant metabolites are plotted, consistent with the
        corrected significance rule.
    """

    results = analysis["results"]
    threshold = fdr_threshold_as_correlation(results, analysis["n_samples"])

    axis.set_xscale("log")

    if threshold is None:
        # Nothing survives FDR correction in this panel. Say so on the
        # figure rather than silently drawing an empty scatter that looks
        # like a plotting bug.
        axis.set_xlim(1e4, 1e9)
        axis.set_ylim(0, 1)
        axis.text(
            0.5, 0.5,
            "No metabolite reaches FDR < 0.05 in this comparison\n"
            f"(smallest adjusted p-value = "
            f"{results['adj_p_value'].min():.3f}, n = {analysis['n_samples']})",
            transform=axis.transAxes, ha="center", va="center", fontsize=12,
            color="dimgrey"
        )
        axis.set_xlabel("Peak area (abundance) - log scale")
        axis.set_ylabel("|correlation with group| (consistency)")
        axis.set_title(
            "No metabolites pass FDR < 0.05 in this comparison — "
            "see write-up for what this does and doesn't mean",
            fontsize=11, pad=8
        )
        return None

    (x_min, x_max), (y_min, y_max) = compute_limits(results, threshold)
    axis.set_xlim(x_min, x_max)
    axis.set_ylim(y_min, y_max)

    plotted = results[results["significant_FDR"]]

    scatter = axis.scatter(
        plotted["area"], plotted["abs_r"],
        s=np.clip(plotted["VIP"] * 22, 4, 260),
        c=np.clip(plotted["VIP"], 0, VIP_COLOUR_CAP),
        cmap="viridis", vmin=0, vmax=VIP_COLOUR_CAP,
        alpha=0.78, edgecolors="none", zorder=3, clip_on=True,
    )

    draw_equal_vip_curves(axis, results)

    axis.axhline(threshold, linestyle="--", color="grey", linewidth=1.2, zorder=2)
    axis.annotate(
        "FDR < 0.05", xy=(0.012, threshold), xycoords=("axes fraction", "data"),
        xytext=(0, 4), textcoords="offset points", fontsize=9, color="dimgrey",
        ha="left", va="bottom", zorder=6,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.2),
    )

    label_box = dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1.5)

    for _key, (row, descriptor) in archetypes.items():
        if row is None:
            continue
        if not (y_min <= row["abs_r"] <= y_max):
            continue
        if not (x_min <= row["area"] <= x_max):
            continue

        fx, fy = axes_fraction(axis, row["area"], row["abs_r"])
        dx, ha = (-18, "right") if fx > 0.55 else (18, "left")
        if fy > 0.80:
            dy, va = -14, "top"
        elif fy < 0.20:
            dy, va = 14, "bottom"
        else:
            dy, va = 11, "bottom"

        name = str(row["Metabolite"])[:26]
        line_gap = 9.5

        if va == "bottom":
            name_offset = (dx, dy + line_gap)
            descriptor_offset = (dx, dy)
            anchored = "descriptor"
        else:
            name_offset = (dx, dy)
            descriptor_offset = (dx, dy - line_gap)
            anchored = "name"

        arrow = dict(arrowstyle="-", lw=0.6, color="grey", shrinkA=2, shrinkB=3)

        axis.annotate(name, xy=(row["area"], row["abs_r"]), textcoords="offset points",
                       xytext=name_offset, fontsize=8, fontweight="bold", ha=ha, va=va,
                       zorder=7, annotation_clip=True, bbox=label_box,
                       arrowprops=arrow if anchored == "name" else None)

        axis.annotate(descriptor, xy=(row["area"], row["abs_r"]), textcoords="offset points",
                       xytext=descriptor_offset, fontsize=7, style="italic", color="#444444",
                       ha=ha, va=va, zorder=7, annotation_clip=True, bbox=label_box,
                       arrowprops=arrow if anchored == "descriptor" else None)

    axis.set_xlabel("Peak area (abundance) - log scale")
    axis.set_ylabel("|correlation with group| (consistency)")
    axis.set_title(
        f"VIP depends on abundance and consistency: only FDR-significant "
        f"metabolites are displayed (BH-FDR < {ALPHA_FDR:g}, "
        f"{len(plotted)} of {len(results)} metabolites)",
        fontsize=11, pad=8,
    )

    return scatter


def build_figure(analysis, save_as=None):
    results = analysis["results"]
    archetypes = pick_archetypes(results)

    figure, axis = plt.subplots(figsize=(10, 7), layout="constrained")
    scatter = panel_map(axis, analysis, archetypes)

    if scatter is not None:
        colour_bar = figure.colorbar(scatter, ax=axis, fraction=0.046, pad=0.02)
        colour_bar.set_label("VIP")

    validation = analysis["validation"]
    perm_note = ""
    if validation:
        perm_note = (
            f"  |  PLS-DA permutation p = {validation['permutation_p']:.3f} "
            f"(R2Y = {validation['R2Y']:.2f})"
        )

    figure.suptitle(
        f"{analysis['tissue']}: {analysis['control_group']} vs "
        f"{analysis['case_group']} ({len(results)} annotated metabolites)"
        f"{perm_note}",
        fontsize=12,
    )

    if save_as:
        path = os.path.join(OUTPUT_DIR, save_as)
        figure.savefig(path, dpi=200, facecolor="white")
        plt.close(figure)
        return path

    return figure


# ==========================================================================
# 14. MAIN
# ==========================================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    all_results = []
    all_significant = []
    validation_rows = []

    for tissue, case in COMPARISONS:
        analysis = analyse(tissue, case)
        results = analysis["results"]
        all_results.append(results)

        sig = significant_metabolites(analysis)
        all_significant.append(sig)

        v = analysis["validation"]
        n_tested = len(results)
        n_sig = int(results["significant_FDR"].sum())

        print(f"\n{tissue}: {analysis['control_group']} vs {case}")
        print(f"  n = {analysis['n_samples']}, components = {analysis['n_components']}")
        print(f"  metabolites tested = {n_tested}")
        print(f"  raw p < 0.05        = {(results['p_value'] < 0.05).sum()}")
        print(f"  FDR (BH) < {ALPHA_FDR:g}      = {n_sig}")
        print(f"  smallest adj. p     = {results['adj_p_value'].min():.4f}")
        print(f"  PLS-DA R2Y = {v['R2Y']:.3f}, permutation p = {v['permutation_p']:.3f}")

        path = build_figure(analysis, save_as=f"panel1_vip_vs_pvalue_{tissue}_{case}.png")
        print(f"  -> {path}")

        validation_rows.append({
            "Tissue": tissue,
            "Comparison": f"control vs {case}",
            "n_samples": analysis["n_samples"],
            "n_tested": n_tested,
            "n_significant_FDR": n_sig,
            "n_plotted": n_sig,
            "criterion": f"BH-FDR < {ALPHA_FDR:g} (Welch t-test on log2 peak area)",
            "PLS_DA_R2Y": round(v["R2Y"], 4),
            "PLS_DA_permutation_p": round(v["permutation_p"], 4),
        })

    combined_results = pd.concat(all_results, ignore_index=True)
    combined_significant = pd.concat(all_significant, ignore_index=True)
    validation_table = pd.DataFrame(validation_rows)

    export_path = os.path.join(OUTPUT_DIR, "significant_metabolites_by_comparison.xlsx")
    with pd.ExcelWriter(export_path) as writer:
        for (tissue, case), sig in zip(COMPARISONS, all_significant):
            sheet_name = f"{tissue}_{case}"[:31]
            sig.to_excel(writer, sheet_name=sheet_name, index=False)
        combined_results.to_excel(writer, sheet_name="all_metabolites_all_tests", index=False)
        validation_table.to_excel(writer, sheet_name="validation_summary", index=False)

    print(f"\nExported: {export_path}")
    print("\nValidation summary:")
    print(validation_table.to_string(index=False))


if __name__ == "__main__":
    main()
