import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def test_sparse_thresholds_keep_negative_comparison():
    data=json.loads((ROOT/'results/sparse_group_sensitivity.json').read_text())
    thresholds=data['thresholds']
    assert thresholds[0]['retained_proteins']==24
    assert thresholds[0]['retained_variants']==66
    assert thresholds[1]['retained_proteins']==17
    assert thresholds[4]['retained_proteins']==3
    assert thresholds[4]['gnn_better_proteins']==0
    assert thresholds[0]['row_weighted_delta_mae']>0
