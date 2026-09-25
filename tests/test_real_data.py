from pathlib import Path
from protein_rescue.features import load_graphs, graph_for_row
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_real_s669_strict_mapping():
    path = ROOT / 'data/raw/S669.csv'
    if not (ROOT / 'data/raw/1A0F.pdb').exists():
        return  # source files fetched explicitly; hermetic CI doesn't contact RCSB
    df = pd.read_csv(path)
    graph = graph_for_row(df.iloc[0], ROOT / 'data/raw')
    assert graph is not None and graph['pdb'] == '1A0F'
    assert graph['x'].shape[1] == 42 and len(graph['x']) <= 32
    assert abs(graph['y'] - (-1.8)) < 1e-8
