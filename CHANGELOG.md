# 0.1.0rc5

- Add an optional Welch interval and a fixed synthetic comparison of interval coverage and width.

# 0.1.0rc4

- Record the bootstrap seed and contrast direction without changing the donor-level estimator or interval method.

# 0.1.0rc3

- Added independent aggregation and exact discrete donor-resampling checks.
- Ran a frozen synthetic coverage study; documented small-donor undercoverage rather than tuning the estimator.
- Marked percentile intervals explicitly nominal and without a coverage guarantee; included reproducible design/results as resources.

# 0.1.0rc2

- Added explicit schema, nonempty identifier and design-boundary validation.
- Added leave-one-independent-unit-out effect ranges as fragility diagnostics.
- Made CLI arguments package-specific and preserved string IDs in CSV ingestion.
- Validated and serialized outputs before reserving a new output directory.
- Documented conditional interval assumptions and the limits of small independent-unit samples.
- Added regression tests for numerical/reference answers, privacy-relevant identifier preservation and invalid-input behavior.

# 0.1.0rc1

Initial review candidate with portable API, explicit input schema and synthetic example.

The rc2 contrast rejects extra unrequested conditions and noninteger bootstrap draw counts. Equal field/donor weights are preserved.
