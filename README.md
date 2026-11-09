# athena-humanquant

Hierarchical donor-level measurement contrasts. Author: Renato Socodato. Version **0.1.0rc3**, prepared for review; not published.

## Install and run

```sh
python -m pip install .
athena-humanquant example --output example_output
athena-humanquant run --input examples/input.csv --output analysis_output --treatment treated --control control --draws 4000 --seed 2026
python -m pip install '.[test]'
python -m pytest tests -q
```

Both examples use synthetic values. Input schemas and assumptions are in
[INPUTS.md](INPUTS.md). Output directories must be new. Input validation and
serialization finish in a temporary staging directory before the requested
output directory is reserved. Existing outputs are refused; invalid inputs do
not leave an empty output directory.

The Python API is `athena_humanquant.core`. No network, cloud account, original grant
tree, external engine or private dataset is required. The core dependency set
is NumPy and pandas. Code and synthetic examples are under MIT; no human data
or private model is redistributed.

## Scientific scope

Measurements average within fields, fields within donors, and donors within conditions. Percentile intervals resample independent donors. This is an unadjusted descriptive contrast, not an LMM or a causal estimate. Small donor n can produce fragile intervals; see INPUTS.md.

When each condition has at least three independent donors, leave-one-donor-out
effect extremes describe fragility. They are not confidence limits and do not validate
the model assumptions. Tests use known numerical answers, invalid inputs,
independence boundaries and synthetic CLI examples. They do not establish
clinical, biomarker or therapeutic validity.

## Independent reference and interval calibration

Independent dictionary aggregation and exact discrete donor-resampling checks
are included in the tests. A fixed-design synthetic calibration found nominal
95% percentile coverage near 66% with two donors/group, 83–87% with four, and
94% with sixteen under the documented normal hierarchical population. These
results expose an inferential limitation; the estimator was not tuned to make
coverage pass. See [the complete design, results and reproduction command](validation/README.md).
The returned interval is explicitly nominal and does not carry a coverage guarantee.

## Relation to the original scientific workflow

This is a general implementation extracted from a scientific use case. It does
not reproduce every original study-specific estimate. No private donor aliases,
project-specific paths, exclusion decisions or source-data trees are shipped.
Future user-supplied identifiers may appear in output tables; review those
outputs before sharing them publicly.

An optional donor-level Welch interval is available through the `welch` extra
and `--interval welch`. The [fixed simulation comparison](validation/README.md)
reports its coverage and the width trade-off with small donor samples.
