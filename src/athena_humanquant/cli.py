DESCRIPTION = 'Hierarchical donor-level measurement contrasts'

"""Offline command-line interface with validated, non-overwriting outputs."""
import argparse
import json
from pathlib import Path
import shutil
import tempfile
import pandas as pd
from . import core


def main(argv=None):
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    commands = parser.add_subparsers(dest='command', required=True)
    example = commands.add_parser('example', help='Run synthetic data, without network')
    example.add_argument('--output', type=Path, required=True)
    run = commands.add_parser('run', help='Analyze a CSV using the documented input schema')
    run.add_argument('--input', type=Path, required=True)
    run.add_argument('--output', type=Path, required=True)
    run.add_argument('--treatment', default='treated')
    run.add_argument('--control', default='control')
    run.add_argument('--draws', type=int, default=4000)
    run.add_argument('--seed', type=int, default=2026)
    run.add_argument('--interval', choices=('percentile', 'welch'), default='percentile')
    args = parser.parse_args(argv)
    output = args.output
    if output.exists():
        parser.error('Output already exists; choose a new directory')
    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.athena-stage-', dir=output.parent) as tmp:
            args.output = Path(tmp)
            result = analyze(args)
            (args.output / 'result.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
            # Reserve exclusively only after all validation and serialization succeed.
            output.mkdir(exist_ok=False)
            for item in args.output.iterdir():
                shutil.move(str(item), output / item.name)
    except (ValueError, KeyError, TypeError, OSError) as exc:
        parser.error(str(exc))
    return 0

def analyze(args):
    if args.command == 'example':
        frame = pd.DataFrame([{'donor': f'synthetic_{i}', 'field': str(j),
            'condition': 'control' if i < 4 else 'treated', 'value': 10 + i + j / 10}
            for i in range(8) for j in range(3)])
        treatment, control, draws, seed = 'treated', 'control', 4000, 2026
    else:
        frame = pd.read_csv(args.input, keep_default_na=False, dtype={'donor': 'string', 'field': 'string', 'condition': 'string'})
        treatment, control, draws, seed = args.treatment, args.control, args.draws, args.seed
    result, donors = core.donor_difference(frame, treatment, control, draws=draws, seed=seed,
                                           interval=getattr(args, 'interval', 'percentile'))
    donors.to_csv(args.output / 'donors.csv', index=False)
    return {'data_kind': 'synthetic' if args.command == 'example' else 'user_supplied', **result}
