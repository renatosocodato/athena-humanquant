"""Run the fixed synthetic calibration design; no coverage-based tuning or gate."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from athena_humanquant.core import donor_difference


def simulate_frame(rng, n, scenario, design):
    records = []
    for condition in ['control', 'treated']:
        for donor in range(n):
            donor_mean = (design['true_condition_difference'] if condition == 'treated' else 0.)
            donor_mean += rng.normal(0., design['donor_sd'])
            fields = 3 if scenario == 'balanced' else int(rng.integers(1, 6))
            for field in range(fields):
                field_mean = donor_mean + rng.normal(0., design['field_sd'])
                cells = 5 if scenario == 'balanced' else int(rng.choice([1, 3, 15]))
                for value in field_mean + rng.normal(0., design['cell_sd'], size=cells):
                    records.append({'donor': f'{condition}_{donor}', 'field': f'f{field}',
                                    'condition': condition, 'value': float(value)})
    return pd.DataFrame(records)


def wilson(covered, total):
    z = 1.95996398454
    p = covered / total
    denominator = 1. + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    width = z * np.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return [float(center - width), float(center + width)]


def run(design, output):
    if output.exists():
        raise ValueError('Output must be a new directory')
    output.mkdir(parents=True)
    truth = design['true_condition_difference']
    raw = []
    summaries = []
    for scenario_index, scenario in enumerate(design['scenarios']):
        for n in design['donors_per_condition']:
            estimates = []
            widths = []
            covered = 0
            for replicate in range(design['replicates_per_scenario']):
                seed = design['base_seed'] + scenario_index * 1000000 + n * 10000 + replicate
                frame = simulate_frame(np.random.default_rng(seed), n, scenario, design)
                result, _ = donor_difference(frame, 'treated', 'control',
                                             draws=design['bootstrap_draws'], seed=seed + 90000000)
                hit = result['ci_low'] <= truth <= result['ci_high']
                covered += int(hit)
                estimates.append(result['difference'])
                widths.append(result['ci_high'] - result['ci_low'])
                raw.append({'scenario': scenario, 'donors_per_condition': n,
                            'replicate': replicate, 'seed': seed, 'estimate': result['difference'],
                            'ci_low': result['ci_low'], 'ci_high': result['ci_high'],
                            'covers_truth': hit})
            total = len(estimates)
            p = covered / total
            error = np.asarray(estimates) - truth
            summary = {'scenario': scenario, 'donors_per_condition': n,
                       'replicates': total, 'coverage': p,
                       'coverage_monte_carlo_se': float(np.sqrt(p * (1 - p) / total)),
                       'coverage_wilson_95': wilson(covered, total),
                       'mean_estimate': float(np.mean(estimates)), 'bias': float(np.mean(error)),
                       'rmse': float(np.sqrt(np.mean(error ** 2))),
                       'mean_interval_width': float(np.mean(widths))}
            summaries.append(summary)
            print(json.dumps(summary), flush=True)
    (output / 'design.json').write_text(json.dumps(design, indent=2) + '\n')
    (output / 'summary.json').write_text(json.dumps({'data_kind': 'synthetic', 'results': summaries}, indent=2) + '\n')
    pd.DataFrame(raw).to_csv(output / 'replicates.csv', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    design = json.loads((Path(__file__).parent / 'CALIBRATION_DESIGN.json').read_text())
    run(design, args.output)
