# HumanQuant input and estimator

CSV columns: `donor,field,condition,value`. Each row is a measurement/cell.
Field identifiers are nested within donor; reusing `field_1` in different donors
is permitted. Donor, field and condition identifiers are read as strings,
preserving leading zeros. Whitespace is trimmed; empty identifiers are rejected.
Numeric values must be finite. No donor exclusions or private aliases are built
in. The input must contain exactly the two requested conditions; additional
conditions are rejected instead of silently ignored.

Cells average within a field; fields receive equal weight within a donor;
donors receive equal weight within a condition. This deliberately differs from
weighting a donor by its number of cells or weighting a field by its number of
measurements. Confirm that equal field weights match the intended sampling
design. The tool does not remove duplicate rows because multiple identical
measurements can be valid; prevent accidental duplication during acquisition.

The reported difference is treatment minus control. The percentile bootstrap
resamples donors within each condition, with a fixed seed and no cell-level
resampling. At least two independent donors per group are required to execute,
but this minimum is not evidence of adequate precision. With small donor n,
bootstrap distributions are discrete, intervals can be narrow or unstable, and
the fixed synthetic calibration observed substantial undercoverage at small donor
n (about 66% at two donors/group in both tested layouts). The `interval` label is
explicitly nominal and `coverage_guarantee` is false; the complete population,
simulation uncertainty and results are documented in `validation/README.md`.
If both groups have at least three
donors, leave-one-donor-out difference extremes describe influence; these are
not confidence limits. No p-value or significance gate is supplied.

This is an unadjusted descriptive contrast. Shared/paired donors, repeated
longitudinal measurements, batch effects, sex imbalance and other confounding
require an appropriate design/model before causal or disease interpretation.
Each donor must belong to one condition; paired designs are rejected.

CLI options `--draws` (integer >=100) and `--seed` control resampling; the API
accepts the same arguments. The result records the seed, number of draws, and
treatment and control labels. Only donor and condition identifiers appear in the
donor table. Review identifiers before publicly sharing future
analysis outputs; only synthetic identifiers are distributed with this package.

### Optional Welch interval

Install the `welch` extra and pass `--interval welch` to the CLI or
`interval='welch'` to `donor_difference`. The point estimate retains the
same cell, field and donor weights. Welch uses the variability of the
independent donor means in each condition and reports the standard error
and approximate degrees of freedom. The interval is nominal and the
two-zero-variance case raises an error. Bootstrap settings apply to the
percentile interval. See the simulation comparison in `validation/README.md`
for coverage and width across the tested donor designs.
