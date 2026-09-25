import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def test_strata_partition_and_direction():
    x=json.loads((ROOT/'results/error_strata.json').read_text())
    assert sum(s['n'] for s in x['bins'])==x['n_test']==149
    assert x['gnn_better_rows']+x['comparator_better_rows']==149
    assert x['gnn_better_pdb_groups']+x['comparator_better_pdb_groups']==13
    assert x['bins'][-1]['n']==11
    assert x['bins'][-1]['delta_mae']>0
