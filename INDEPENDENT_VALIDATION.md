# Independent implementation and scientific checks

An independent dictionary aggregation and an exact discrete two-donor bootstrap distribution agree with the implemented estimand. A design fixed before simulation outcomes examined 3,200 synthetic experiments. Nominal 95% percentile interval coverage was about 66% at two donors per condition, 83–87% at four, 91–92% at eight and 94% at sixteen, in the tested balanced and unequal measurement layouts.

This supports the donor-level aggregation implementation while exposing substantial interval undercoverage at small donor n in the simulated population. More cells do not replace donor replication. There is no universal sufficient donor count or general 95% coverage guarantee; the API explicitly labels the interval nominal. Frozen design, Monte Carlo uncertainty and reproduction code are in validation/.

These are internal computational checks against independent implementations, public inputs and real engines. No external peer review is claimed. Full execution details are in VALIDATION.md.

## Revision after independent candidate review

Record the bootstrap seed and contrast direction without changing the donor-level estimator or interval method. Hosted execution records for the exact distributions are tracked separately.
