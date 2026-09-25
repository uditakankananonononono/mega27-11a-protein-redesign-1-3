"""Command-line interface for reproducible structure-group evaluation and mutation feature audit."""
import argparse
import json
from pathlib import Path
from .benchmark import run
from .features import graph_for_row
import pandas as pd


def main():
    p = argparse.ArgumentParser(description='Protein rescue mutation stability benchmark')
    sub = p.add_subparsers(dest='command', required=True)
    benchmark = sub.add_parser('benchmark', help='Run PDB-held-out GNN and ridge comparison')
    benchmark.add_argument('--csv', type=Path, default=Path('data/raw/S669.csv'))
    benchmark.add_argument('--structures', type=Path, default=Path('data/raw'))
    benchmark.add_argument('--output', type=Path, default=Path('results/s669_group_holdout.json'))
    benchmark.add_argument('--epochs', type=int, default=80)
    inspect = sub.add_parser('inspect', help='Validate and inspect a mutation graph by S669 row name')
    inspect.add_argument('name')
    inspect.add_argument('--csv', type=Path, default=Path('data/raw/S669.csv'))
    inspect.add_argument('--structures', type=Path, default=Path('data/raw'))
    args = p.parse_args()
    if args.command == 'benchmark':
        result = run(args.csv, args.structures, args.output, args.epochs)
        print(json.dumps(dict(audit=result['audit'], split={k:v for k,v in result['split'].items() if isinstance(v,int)}, metrics=result['metrics']), indent=2))
    elif args.command == 'inspect':
        df = pd.read_csv(args.csv)
        rows = df[df.name == args.name]
        if len(rows) != 1:
            p.error('Mutation identifier not found uniquely')
        g = graph_for_row(rows.iloc[0], args.structures)
        if g is None:
            p.error('Structure missing, sequence mismatch, or insufficient aligned residue graph')
        print(json.dumps(dict(name=g['name'], pdb=g['pdb'], residues=len(g['x']), experimental_ddg=g['y'], contact_edges=int((g['a']>0).sum()))))

if __name__ == '__main__':
    main()
