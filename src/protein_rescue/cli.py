"""Command-line interface for reproducible structure-group evaluation and mutation feature audit."""
import argparse
import json
from pathlib import Path
from .benchmark import run
from .features import graph_for_row
from .audit import audit
from .fireprot import analyze
from .targets import map_site
from .epistasis import double_mutant_cycle, literature_controls
from .filters import safe_to_rank
import pandas as pd


def main():
    p = argparse.ArgumentParser(description='Protein rescue mutation stability benchmark')
    sub = p.add_subparsers(dest='command', required=True)
    benchmark = sub.add_parser('benchmark', help='Run PDB-held-out GNN and ridge comparison')
    benchmark.add_argument('--csv', type=Path, default=Path('data/raw/S669.csv'))
    benchmark.add_argument('--structures', type=Path, default=Path('data/raw'))
    benchmark.add_argument('--output', type=Path, default=Path('results/s669_group_holdout.json'))
    benchmark.add_argument('--epochs', type=int, default=80)
    benchmark.add_argument('--graph-cache', type=Path, default=None)
    inspect = sub.add_parser('inspect', help='Validate and inspect a mutation graph by S669 row name')
    inspect.add_argument('name')
    inspect.add_argument('--csv', type=Path, default=Path('data/raw/S669.csv'))
    inspect.add_argument('--structures', type=Path, default=Path('data/raw'))
    evaluation = sub.add_parser('audit', help='Group-aware error audit of a held-out result')
    evaluation.add_argument('--results', type=Path, default=Path('results/s669_group_holdout.json'))
    evaluation.add_argument('--csv', type=Path, default=Path('data/raw/S669.csv'))
    evaluation.add_argument('--output', type=Path, default=Path('results/group_sensitivity.json'))
    fp = sub.add_parser('fireprot-audit', help='Check experimental mutation numbering and sign discordance')
    fp.add_argument('--csv', type=Path, default=Path('data/fireprot_pdb_ddg.csv'))
    fp.add_argument('--structures', type=Path, default=Path('data/raw'))
    fp.add_argument('--output', type=Path, default=Path('results/fireprot_mapping_audit.json'))
    target = sub.add_parser('map-site', help='Map canonical UniProt position to observed target PDB residue')
    target.add_argument('gene', choices=['TP53','SOD1','PTEN'])
    target.add_argument('position', type=int)
    target.add_argument('--structures',type=Path,default=Path('data/raw'))
    target.add_argument('--sequences',type=Path,default=Path('data/targets'))
    ep = sub.add_parser('cycle', help='Calculate experimental double-mutant interaction; not rescue prediction')
    ep.add_argument('single_a',type=float)
    ep.add_argument('single_b',type=float)
    ep.add_argument('combined',type=float)
    ep.add_argument('--se-a',type=float)
    ep.add_argument('--se-b',type=float)
    ep.add_argument('--se-combined',type=float)
    gate = sub.add_parser('rank-check', help='Abstain unless proposed ranking has verified assay and target context')
    gate.add_argument('--context',type=Path,required=True,help='JSON object with structure, assay, training, validation flags')
    args = p.parse_args()
    if args.command == 'benchmark':
        result = run(args.csv, args.structures, args.output, args.epochs, graph_cache=args.graph_cache)
        print(json.dumps(dict(audit=result['audit'], split={k:v for k,v in result['split'].items() if isinstance(v,int)}, metrics=result['metrics']), indent=2))
    elif args.command == 'audit':
        result = audit(args.results, args.csv, args.output)
        print(json.dumps({k:v for k,v in result.items() if k not in ('grouped',)},indent=2))
    elif args.command == 'fireprot-audit':
        result = analyze(args.csv,args.structures,args.output)
        print(json.dumps({k:v for k,v in result.items() if k not in ('per_pdb','opposite_sign_examples')},indent=2))
    elif args.command == 'rank-check':
        print(json.dumps(safe_to_rank(json.loads(args.context.read_text())),indent=2))
    elif args.command == 'cycle':
        print(json.dumps(double_mutant_cycle(args.single_a,args.single_b,args.combined,args.se_a,args.se_b,args.se_combined),indent=2))
    elif args.command == 'map-site':
        print(json.dumps(map_site(args.gene,args.position,args.structures,args.sequences),indent=2))
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
