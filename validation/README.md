# Finite-donor interval calibration

The implementation matches independent aggregation and an exact two-donor
resampling distribution. That establishes the implemented estimand; it does
not imply that a nominal 95% interval has 95% population coverage.

The design in `CALIBRATION_DESIGN.json` was fixed before the first simulation
outcomes. It uses a condition difference of one, independent normal donor
effects (SD 1), shared field effects (SD 0.4), cell noise (SD 0.3), and independent
measurement counts. Every scenario has 400 independent synthetic experiments,
each with 400 donor-bootstrap draws. The balanced layout has three fields and
five cells per field; the unequal layout has one to five fields and 1/3/15 cells.
No result-dependent change of seed, estimator, exclusion or threshold was made.

| Donors per condition | Balanced coverage | Unequal coverage |
|---:|---:|---:|
| 2 | 66.25% | 66.00% |
| 4 | 86.75% | 83.25% |
| 8 | 91.00% | 91.75% |
| 16 | 93.75% | 94.25% |

At two donors/group, the Wilson intervals for coverage were 61.5–70.7% and
61.2–70.5%; the undercoverage is large relative to Monte Carlo uncertainty.
At sixteen donors/group they were 90.9–95.7% and 91.5–96.1%. These are uncertainty
intervals for the simulation's estimated coverage, not intervals for a biological
effect. Full Monte Carlo standard errors, bias, RMSE and interval widths are in
`CALIBRATION_RESULTS.json`.

The useful conclusion is that the current estimator can produce a meaningful
donor-level point contrast while its percentile limits substantially under-cover
at small donor n in this population. More cells do not supply the missing donor
replication. No universally sufficient donor count is established by these eight
synthetic scenarios; different distributions, variance structures, confounding or
designs can change the result. The function therefore labels the interval as
nominal and reports `coverage_guarantee: false`. It does not silently substitute
another inferential method or manufacture a sample-size threshold.

To reproduce, from an installed source distribution or checkout:

```sh
python validation/calibrate_donor_bootstrap.py --output new_calibration_output
```

The output directory must be new. The script copies the frozen design and writes
the summary plus all synthetic replicate estimates. No network or real human
data are required. The two JSON files are also installed as package resources
under `athena_humanquant/calibration/`.

## Comparison with Welch intervals

A separate fixed simulation evaluated 80,000 experiments across 40 scenarios.
Both methods use the same donor means. The scenarios cross normal, unequal
normal, heavy-tailed and treatment-skewed donor effects with balanced or
unequal field layouts and five donor count pairs.

| Donors (treated/control) | Percentile coverage | Welch coverage | Welch / percentile mean width |
|---|---:|---:|---:|
| 2/2 | 63.7–68.4% | 97.0–98.1% | 6.20–7.10 |
| 4/4 | 83.7–85.8% | 94.8–96.9% | 1.56–1.65 |
| 4/8 | 82.2–87.5% | 93.2–97.0% | 1.43–1.67 |
| 8/8 | 89.6–91.4% | 94.5–96.1% | 1.19–1.22 |
| 16/16 | 91.6–93.6% | 94.9–96.0% | 1.08–1.10 |

Ranges describe the eight simulated populations and layouts at each pair.
The large widths at two donors quantify how little information the design
has about the mean difference. Welch coverage depends on the donor
distribution and variance structure. The percentile interval remains the
default. The optional Welch interval uses independent donor means and
allows different variances. It is undefined when both observed variances
are zero.

Reproduce with `python validation/compare_donor_intervals.py --output comparison.json`
after installing SciPy. The full scenario results and paired Monte Carlo
standard errors are in `INTERVAL_COMPARISON.json`.

Method reference: [Welch, Biometrika 34, 28–35 (1947)](https://doi.org/10.1093/biomet/34.1-2.28).
