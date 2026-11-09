"""Equal-field, equal-donor descriptive contrasts of nested measurements."""
import numpy as np
import pandas as pd


def donor_difference(frame, treatment, control, draws=4000, seed=2026, *, interval='percentile'):
    """Average cells within fields, then fields within independent donors.

    The nominal 95% percentile interval resamples donors within group. It
    summarizes the empirical donor distribution; it is not coverage-calibrated
    for small donor samples (see the included fixed-design calibration).
    Set interval='welch' for an approximate normal-theory interval allowing
    unequal donor variances. All rows must belong to the two requested conditions.
    """
    columns = ['donor', 'field', 'condition', 'value']
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError('Missing required columns: ' + ', '.join(missing))
    if interval not in ('percentile', 'welch'):
        raise ValueError('interval must be percentile or welch')
    if (isinstance(draws, (bool, np.bool_)) or not isinstance(draws, (int, np.integer))
            or draws < 100):
        raise ValueError('draws must be an integer of at least 100')
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError('seed must be a nonnegative integer')
    if not isinstance(treatment, str) or not isinstance(control, str):
        raise ValueError('Condition names must be nonempty strings')
    treatment, control = treatment.strip(), control.strip()
    if not treatment or not control or treatment == control or frame[columns].isna().any().any():
        raise ValueError('Complete observations and distinct nonempty conditions required')
    data = frame.copy()
    for column in ['donor', 'field', 'condition']:
        data[column] = data[column].astype(str).str.strip()
        if data[column].eq('').any():
            raise ValueError('Donor, field and condition identifiers must be nonempty')
    if set(data.condition) != {treatment, control}:
        raise ValueError('Input must contain exactly the requested two conditions')
    try:
        data['value'] = pd.to_numeric(data.value, errors='raise').astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError('Measurements must be numeric') from exc
    if not np.isfinite(data.value.to_numpy()).all():
        raise ValueError('Finite measurements required')
    if (data.groupby('donor').condition.nunique() > 1).any():
        raise ValueError('Conflicting donor condition')
    fields = data.groupby(['donor', 'condition', 'field']).value.mean()
    donors = fields.groupby(['donor', 'condition']).mean().reset_index()
    a = donors.loc[donors.condition.eq(treatment), 'value'].to_numpy()
    b = donors.loc[donors.condition.eq(control), 'value'].to_numpy()
    if min(len(a), len(b)) < 2:
        raise ValueError('At least two independent donors per group required')
    difference = a.mean() - b.mean()
    if interval == 'percentile':
        rng = np.random.default_rng(seed)
        boot = np.array([rng.choice(a, len(a), replace=True).mean() -
                         rng.choice(b, len(b), replace=True).mean() for _ in range(draws)])
        low, high = np.quantile(boot, [.025, .975])
        interval_info = {
            'interval': 'donor_percentile_bootstrap_nominal_95_percent',
            'bootstrap_draws': int(draws), 'bootstrap_seed': int(seed),
            'precision_limit': 'Small donor samples give discrete, fragile bootstrap intervals; fixed-design synthetic calibration showed undercoverage. See calibration documentation.',
        }
    else:
        try:
            from scipy.special import stdtrit
        except ImportError as error:
            raise ValueError('Welch requires SciPy. Install the welch extra.') from error
        scale = max(float(np.abs(a).max()), float(np.abs(b).max()))
        if not np.isfinite(scale):
            raise ValueError('Measurement scale exceeds supported numeric range')
        if scale == 0:
            raise ValueError('Welch requires nonzero observed variance in at least one group')
        va = np.var(a / scale, ddof=1) / len(a)
        vb = np.var(b / scale, ddof=1) / len(b)
        total = va + vb
        if total == 0:
            raise ValueError('Welch requires nonzero observed variance in at least one group')
        df = 1 / ((va / total)**2 / (len(a) - 1) + (vb / total)**2 / (len(b) - 1))
        se = np.sqrt(total) * scale
        half_width = stdtrit(df, .975) * se
        low, high = difference - half_width, difference + half_width
        if se == 0 or low == high:
            raise ValueError('Welch interval is below supported numeric resolution')
        interval_info = {
            'interval': 'donor_welch_satterthwaite_nominal_95_percent',
            'standard_error': float(se), 'degrees_of_freedom': float(df),
            'precision_limit': 'Approximate Welch interval for independent donor means. With few donors, intervals can be wide and coverage depends on the population and design.',
        }
    values = [difference, low, high]
    if not np.isfinite(values).all():
        raise ValueError('Measurement scale exceeds supported numeric range')
    result = {
        'difference': float(values[0]), 'ci_low': float(values[1]), 'ci_high': float(values[2]),
        'n_treatment_donors': len(a), 'n_control_donors': len(b),
        'weighting': 'equal fields within donor; equal donors within condition',
        'interval_nominal_level': 0.95, 'coverage_guarantee': False,
        'treatment': treatment, 'control': control,
        **interval_info,
    }
    if min(len(a), len(b)) >= 3:
        leave_out = [np.delete(a, i).mean() - b.mean() for i in range(len(a))]
        leave_out += [a.mean() - np.delete(b, i).mean() for i in range(len(b))]
        result['leave_one_donor_out_min'] = float(min(leave_out))
        result['leave_one_donor_out_max'] = float(max(leave_out))
    return result, donors
