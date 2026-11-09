import json

import pandas as pd
import pytest
from scipy.stats import ttest_ind

from athena_humanquant.cli import main
from athena_humanquant.core import donor_difference


def observations():
    return pd.DataFrame([
        {'donor': f'{condition}_{i}', 'condition': condition, 'field': str(field),
         'value': donor_value + offset}
        for condition, means in [('control', [0., 1., 2., 3.]), ('treated', [2., 7., 9.])]
        for i, donor_value in enumerate(means)
        for field, offset in enumerate([-1., 1.])
    ])


def test_cli_welch_agrees_with_independent_scipy_reference(tmp_path):
    source = tmp_path / 'observations.csv'
    output = tmp_path / 'result'
    observations().to_csv(source, index=False)
    assert main(['run', '--input', str(source), '--output', str(output), '--interval', 'welch']) == 0
    result = json.loads((output / 'result.json').read_text())
    reference = ttest_ind([2., 7., 9.], [0., 1., 2., 3.], equal_var=False)
    low, high = reference.confidence_interval()
    assert result['difference'] == pytest.approx(4.5)
    assert result['degrees_of_freedom'] == pytest.approx(reference.df)
    assert [result['ci_low'], result['ci_high']] == pytest.approx([low, high])
    assert 'bootstrap_draws' not in result


def test_welch_does_not_assign_zero_uncertainty_to_constant_samples():
    with pytest.raises(ValueError, match='nonzero observed variance'):
        donor_difference(observations().assign(value=3), 'treated', 'control', interval='welch')
