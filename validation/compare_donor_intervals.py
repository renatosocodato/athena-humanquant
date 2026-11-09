"""Fixed synthetic comparison of two nominal 95% intervals for donor means."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import stdtrit

SEED = 2026110807
REPLICATES = 2000
DRAWS = 999
TRUTH = 1.0
LAWS = ("normal", "unequal_normal", "heavy_tail", "skewed_treatment")
SIZES = ((2, 2), (4, 4), (4, 8), (8, 8), (16, 16))
LAYOUTS = ("balanced", "unequal_fields")


def generate(rng, n, law, layout, treated):
    shape = (REPLICATES, n)
    donor = rng.normal(size=shape)
    if law == "heavy_tail":
        donor = rng.standard_t(3, size=shape) / np.sqrt(3)
    elif law == "skewed_treatment" and treated:
        donor = rng.exponential(size=shape) - 1
    elif law == "unequal_normal" and treated:
        donor *= 2
    # Normal field and cell noise aggregated at exactly the cell → field → donor hierarchy.
    if layout == "balanced":
        measurement_sd = np.sqrt(0.4**2 / 3 + 0.3**2 / 15)
    else:
        fields = rng.integers(1, 6, size=shape)
        cells = rng.choice([1, 3, 15], size=(*shape, 5))
        active = np.arange(5) < fields[..., None]
        measurement_sd = np.sqrt(0.4**2 / fields + 0.3**2 * (active / cells).sum(axis=-1) / fields**2)
    return donor + rng.normal(size=shape) * measurement_sd + (TRUTH if treated else 0)


def simulate(rng, law, layout, na, nb):
    a, b = generate(rng, na, law, layout, True), generate(rng, nb, law, layout, False)
    estimates = a.mean(1) - b.mean(1)
    va, vb = a.var(1, ddof=1) / na, b.var(1, ddof=1) / nb
    df = (va + vb)**2 / (va**2 / (na - 1) + vb**2 / (nb - 1))
    half = stdtrit(df, .975) * np.sqrt(va + vb)
    welch_hit = np.abs(estimates - TRUTH) <= half
    lo, hi = np.empty(REPLICATES), np.empty(REPLICATES)
    for start in range(0, REPLICATES, 100):
        stop = min(start + 100, REPLICATES)
        idx = np.arange(start, stop)[:, None, None]
        ma = a[idx, rng.integers(na, size=(stop-start, DRAWS, na))].mean(2)
        mb = b[idx, rng.integers(nb, size=(stop-start, DRAWS, nb))].mean(2)
        lo[start:stop], hi[start:stop] = np.quantile(ma-mb, [.025, .975], axis=1)
    boot_hit = (lo <= TRUTH) & (TRUTH <= hi)
    paired = welch_hit.astype(float) - boot_hit
    return {
        "law": law, "layout": layout, "treatment_n": na, "control_n": nb,
        "percentile_coverage": float(boot_hit.mean()), "welch_coverage": float(welch_hit.mean()),
        "coverage_change": float(paired.mean()),
        "paired_mcse": float(paired.std(ddof=1) / np.sqrt(REPLICATES)),
        "percentile_mean_width": float((hi-lo).mean()), "welch_mean_width": float((2*half).mean()),
        "bias": float((estimates-TRUTH).mean()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = np.random.default_rng(SEED)
    results = [simulate(rng, law, layout, *size) for law in LAWS for layout in LAYOUTS for size in SIZES]
    with args.output.open("x") as f:
        json.dump({"seed": SEED, "replicates_per_scenario": REPLICATES, "bootstrap_draws": DRAWS,
                   "true_mean_difference": TRUTH, "results": results}, f, indent=2)


if __name__ == "__main__":
    main()
