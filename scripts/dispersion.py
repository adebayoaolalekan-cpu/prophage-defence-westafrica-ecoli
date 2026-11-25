"""
Estimate the negative binomial dispersion from the observed pilot counts.

This is the number Section 10.1 of the protocol is held up on. The feasibility
simulation fixed theta rather than estimating it, and condition 4 of Section
8.0 is satisfied for the human versus poultry contrast at a rate ratio of 2.0
but not at 1.5 unless the pilot returns theta of 5 or better.

Theta is estimated three ways, because they can disagree and the disagreement
is informative rather than something to average away.

  1. Method of moments on the raw counts. Var = mu + mu^2 / theta.
  2. Maximum likelihood on the intercept-only negative binomial.
  3. Maximum likelihood with the Section 10.1 covariates included, which is the
     estimate the re-simulation should use, since dispersion conditional on
     the model is what the power calculation needs.

A confidence interval on theta comes from a nonparametric bootstrap over
genomes, resampled within lineage so the interval carries the clustering the
analysis model carries.

Separate estimates are reported for the inference block and the catalogue
block. Pooling them would let the catalogue block's wider source mix inflate
the dispersion that the human versus poultry contrast is judged against.
"""

import json
import math
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import write_tsv

THETA_THRESHOLD = 5.0   # protocol Section 10.1
N_BOOT = 2000
SEED = 20261125


def mom_theta(y):
    """Method of moments. Returns None when the data are not overdispersed."""
    y = np.asarray(y, dtype=float)
    mu, var = y.mean(), y.var(ddof=1) if len(y) > 1 else 0.0
    if mu <= 0 or var <= mu:
        return None
    return float(mu ** 2 / (var - mu))


def nb_mle_theta(y, X=None):
    """
    Profile the negative binomial over alpha = 1/theta. statsmodels needs alpha
    fixed, so the likelihood is evaluated on a grid and then refined. Returns
    (theta, loglik) or (None, None).
    """
    import statsmodels.api as sm

    y = np.asarray(y, dtype=float)
    X = np.ones((len(y), 1)) if X is None else np.asarray(X, dtype=float)

    def ll(alpha):
        try:
            fam = sm.families.NegativeBinomial(alpha=alpha)
            res = sm.GLM(y, X, family=fam).fit()
            return res.llf
        except Exception:
            return -np.inf

    grid = np.exp(np.linspace(math.log(1e-3), math.log(50.0), 120))
    lls = np.array([ll(a) for a in grid])
    if not np.isfinite(lls).any():
        return None, None
    i = int(np.nanargmax(lls))
    lo = grid[max(0, i - 1)]
    hi = grid[min(len(grid) - 1, i + 1)]
    fine = np.exp(np.linspace(math.log(lo), math.log(hi), 60))
    flls = np.array([ll(a) for a in fine])
    j = int(np.nanargmax(flls))
    alpha = float(fine[j])
    return (float("inf") if alpha < 1e-3 else 1.0 / alpha), float(flls[j])


def design(df):
    """
    Section 10.1 covariates: log N50, log contig count, completeness,
    collection year and source. Columns absent from the pilot output are
    dropped and the omission is reported rather than filled in.
    """
    cols, names = [np.ones(len(df))], ["intercept"]
    used, skipped = [], []

    def add(name, vals):
        v = pd.to_numeric(vals, errors="coerce")
        if v.isna().all() or v.nunique(dropna=True) < 2:
            skipped.append(name)
            return
        v = v.fillna(v.median())
        cols.append(v.to_numpy(dtype=float))
        names.append(name)
        used.append(name)

    add("log_n50", np.log(pd.to_numeric(df["contig_n50"], errors="coerce")))
    add("log_contigs", np.log(pd.to_numeric(df["number_of_contigs"],
                                            errors="coerce")))
    add("completeness", df["completeness"])
    add("year", pd.to_numeric(
        df["collection_date"].astype(str).str.slice(0, 4), errors="coerce"))

    src = df["ecological_category"].astype(str)
    levels = sorted(x for x in src.unique() if x and x != "nan")
    if len(levels) > 1:
        for lv in levels[1:]:
            cols.append((src == lv).to_numpy(dtype=float))
            names.append(f"source_{lv}")
            used.append(f"source_{lv}")
    else:
        skipped.append("source")

    return np.column_stack(cols), names, used, skipped


def bootstrap_theta(df, ycol, rng, n=N_BOOT):
    """Cluster bootstrap over lineage. Falls back to genome-level resampling
    when no lineage assignment is available."""
    lin = df["lineage"].astype(str)
    groups = [g.index.to_numpy() for _, g in df.groupby(lin)] \
        if lin.notna().any() and lin.nunique() > 1 else None

    out = []
    for _ in range(n):
        if groups:
            pick = rng.integers(0, len(groups), len(groups))
            idx = np.concatenate([groups[k] for k in pick])
        else:
            idx = rng.integers(0, len(df), len(df))
        y = df.loc[idx, ycol].to_numpy(dtype=float)
        t = mom_theta(y)
        if t is not None and np.isfinite(t):
            out.append(t)
    if len(out) < 50:
        return None, None, len(out)
    a = np.asarray(out)
    return float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5)), len(out)


def describe(df, label, ycol, rng):
    y = pd.to_numeric(df[ycol], errors="coerce").dropna().to_numpy(dtype=float)
    if len(y) == 0:
        return {"set": label, "n": 0, "note": "no genomes"}

    X, names, used, skipped = design(df)
    t_mom = mom_theta(y)
    t_int, ll_int = nb_mle_theta(y)
    t_adj, ll_adj = nb_mle_theta(
        pd.to_numeric(df[ycol], errors="coerce").fillna(0).to_numpy(float), X)
    lo, hi, nb = bootstrap_theta(df.reset_index(drop=True), ycol, rng)

    mu, var = float(y.mean()), float(y.var(ddof=1)) if len(y) > 1 else 0.0
    primary = t_adj if t_adj is not None else (t_int if t_int is not None else t_mom)

    return {
        "set": label,
        "n": int(len(y)),
        "mean_count": round(mu, 4),
        "variance": round(var, 4),
        "variance_to_mean_ratio": round(var / mu, 4) if mu > 0 else None,
        "zero_count_genomes": int((y == 0).sum()),
        "max_count": int(y.max()),
        "theta_method_of_moments": None if t_mom is None else round(t_mom, 4),
        "theta_mle_intercept_only": None if t_int is None else round(t_int, 4),
        "theta_mle_covariate_adjusted": None if t_adj is None else round(t_adj, 4),
        "theta_for_resimulation": None if primary is None else round(primary, 4),
        "theta_bootstrap_ci_low": None if lo is None else round(lo, 4),
        "theta_bootstrap_ci_high": None if hi is None else round(hi, 4),
        "bootstrap_replicates_used": nb,
        "covariates_used": used,
        "covariates_unavailable": skipped,
        "poisson_adequate": bool(t_mom is None),
        "meets_theta_threshold_5": (
            None if primary is None else bool(primary >= THETA_THRESHOLD)),
        "loglik_intercept_only": None if ll_int is None else round(ll_int, 3),
        "loglik_covariate_adjusted": None if ll_adj is None else round(ll_adj, 3),
    }


def main(counts_path, out_json, out_tsv):
    rng = np.random.default_rng(SEED)
    df = pd.read_csv(counts_path, sep="\t", dtype=str)
    for c in ("n_prophage_primary", "n_prophage_intersection"):
        df[c] = pd.to_numeric(df[c], errors="coerce")

    sets = {
        "inference_block_high_contiguity": df[
            (df["block"] == "A_inference")
            & (df["passes_high_contiguity"] == "YES")],
        "inference_block_all": df[df["block"] == "A_inference"],
        "catalogue_block_high_contiguity": df[
            (df["block"] == "B_catalogue")
            & (df["passes_high_contiguity"] == "YES")],
        "catalogue_block_all": df[df["block"] == "B_catalogue"],
        "whole_pilot": df,
    }

    results = []
    for label, sub in sets.items():
        if len(sub) == 0:
            results.append({"set": label, "n": 0, "note": "no genomes"})
            continue
        r = describe(sub.reset_index(drop=True), label,
                     "n_prophage_primary", rng)
        r["measure"] = "n_prophage_primary"
        results.append(r)
        s = describe(sub.reset_index(drop=True),
                     label + "__intersection_sensitivity",
                     "n_prophage_intersection", rng)
        s["measure"] = "n_prophage_intersection"
        results.append(s)

    primary_set = next(
        (r for r in results
         if r["set"] == "inference_block_high_contiguity"
         and r.get("measure") == "n_prophage_primary"), None)

    payload = {
        "seed": SEED,
        "bootstrap_replicates_requested": N_BOOT,
        "theta_threshold_section_10_1": THETA_THRESHOLD,
        "estimates": results,
        "verdict": _verdict(primary_set),
    }
    with open(out_json, "w") as fh:
        json.dump(payload, fh, indent=2)

    flat_fields = ["set", "measure", "n", "mean_count", "variance",
                   "variance_to_mean_ratio", "zero_count_genomes", "max_count",
                   "theta_method_of_moments", "theta_mle_intercept_only",
                   "theta_mle_covariate_adjusted", "theta_for_resimulation",
                   "theta_bootstrap_ci_low", "theta_bootstrap_ci_high",
                   "poisson_adequate", "meets_theta_threshold_5"]
    write_tsv(out_tsv, results, flat_fields)

    sys.stderr.write(json.dumps(payload["verdict"], indent=2) + "\n")


def _verdict(p):
    if p is None or p.get("n", 0) == 0:
        return {"status": "no_estimate",
                "note": "The inference block returned no high-contiguity "
                        "genomes. Section 10.1 cannot be re-run."}
    theta = p.get("theta_for_resimulation")
    lo = p.get("theta_bootstrap_ci_low")
    if theta is None:
        return {"status": "poisson_adequate",
                "note": "The counts show no overdispersion. Re-run Section "
                        "10.1 with a Poisson model, which is more favourable "
                        "than the negative binomial the simulation assumed."}
    if lo is not None and lo >= THETA_THRESHOLD:
        s = "rr_1_5_supported"
        note = ("Theta and the lower bound of its interval both clear 5. "
                "Condition 4 can be claimed at a rate ratio of 1.5. Re-run "
                "Section 10.1 with this theta and the observed lineage "
                "partition to confirm.")
    elif theta >= THETA_THRESHOLD:
        s = "rr_1_5_borderline"
        note = ("Theta clears 5 but the lower bound of its interval does not. "
                "Keep the minimum detectable effect at 2.0 and report the "
                "1.5 result as the sensitivity case.")
    else:
        s = "rr_2_0_only"
        note = ("Theta falls short of 5. The minimum detectable effect stays "
                "at 2.0, as Section 10.1 already states. If theta is very low "
                "the contrast is reported descriptively, per the protocol.")
    return {"status": s, "theta": theta,
            "ci": [lo, p.get("theta_bootstrap_ci_high")], "note": note}


if __name__ == "__main__":
    try:
        sm = snakemake  # noqa: F821
    except NameError:
        sm = None
    if sm is not None:
        main(sm.input.counts, sm.output.json, sm.output.tsv)
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
