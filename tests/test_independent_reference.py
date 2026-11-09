"""Independent dictionary aggregation and exact donor-resampling reference."""
from collections import defaultdict
from itertools import product
import numpy as np
import pandas as pd

from athena_humanquant.core import donor_difference


def independent_donor_means(frame):
    fields = defaultdict(list)
    for row in frame.itertuples(index=False):
        fields[(row.donor, row.condition, row.field)].append(float(row.value))
    donors = defaultdict(list)
    for (donor, condition, _), values in fields.items():
        donors[(donor, condition)].append(sum(values) / len(values))
    return {key: sum(values) / len(values) for key, values in donors.items()}


def test_unequal_cell_and_field_counts_match_independent_aggregation():
    rows = []
    for i in range(6):
        condition = 'control' if i < 3 else 'treated'
        for field in range(1 + i % 3):
            for cell in range(1 + (i + field) % 5):
                rows.append({'donor': f'd{i}', 'condition': condition,
                             'field': f'f{field}', 'value': 2 * i + field + cell / 10})
    frame = pd.DataFrame(rows)
    reference = independent_donor_means(frame)
    result, donors = donor_difference(frame, 'treated', 'control', draws=100)
    for row in donors.itertuples(index=False):
        np.testing.assert_allclose(row.value, reference[(row.donor, row.condition)])
    a = [v for (d, c), v in reference.items() if c == 'treated']
    b = [v for (d, c), v in reference.items() if c == 'control']
    np.testing.assert_allclose(result['difference'], sum(a) / len(a) - sum(b) / len(b))


def test_small_donor_bootstrap_matches_exact_resampling_distribution():
    # With two donors/group there are 4*4 ordered bootstrap combinations.
    # The fixed example's 2.5/97.5% quantiles are the extrema, not a normal CI.
    frame = pd.DataFrame({'donor': ['a', 'b', 'c', 'd'],
                          'field': ['f'] * 4,
                          'condition': ['control'] * 2 + ['treated'] * 2,
                          'value': [0., 2., 4., 7.]})
    exact = [np.mean(a) - np.mean(b)
             for a in product([4., 7.], repeat=2)
             for b in product([0., 2.], repeat=2)]
    result, _ = donor_difference(frame, 'treated', 'control', draws=20000, seed=20261108)
    # Each enumerated combination has probability 1/16. Use the discrete CDF
    # quantile, not linear interpolation of a sample of only 16 order statistics.
    np.testing.assert_allclose([result['ci_low'], result['ci_high']],
                               np.quantile(exact, [.025, .975], method='inverted_cdf'))


def test_duplicate_cells_preserve_reference_donor_estimand():
    frame = pd.DataFrame({'donor': ['a', 'b', 'c', 'd'],
                          'field': ['f'] * 4,
                          'condition': ['control'] * 2 + ['treated'] * 2,
                          'value': [0., 2., 4., 7.]})
    duplicated = pd.concat([frame, pd.concat([frame.iloc[:1]] * 30)], ignore_index=True)
    assert independent_donor_means(frame) == independent_donor_means(duplicated)
    left, _ = donor_difference(frame, 'treated', 'control', draws=400, seed=43)
    right, _ = donor_difference(duplicated, 'treated', 'control', draws=400, seed=43)
    assert left == right
