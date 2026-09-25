import numpy as np
import torch
from protein_rescue.features import align_index, AA
from protein_rescue.model import ResidueGNN,collate
from protein_rescue.benchmark import metrics


def test_alignment_handles_missing_segment():
    assert align_index('ABCDEFGH', 'ABCEFGH', 5) == 4
    assert align_index('ABCDEFGH', 'ABCEFGH', 3) is None


def test_graph_batch_and_forward():
    def g(n):
        x=np.zeros((n,42),np.float32); x[0,40]=1
        return dict(x=x,a=np.eye(n,dtype=np.float32),y=-1.)
    x,a,y=collate([g(5),g(7)])
    assert tuple(x.shape)==(2,7,42)
    assert ResidueGNN()(x,a).shape==(2,)


def test_metrics_exact():
    result=metrics(np.array([1.,2.,3.]),np.array([1.,2.,3.]))
    assert result['mae']==0 and result['spearman']==1


def test_aa_alphabet():
    assert len(AA)==20 and len(set(AA))==20


def test_prediction_invariant_to_padding():
    torch.manual_seed(3)
    x=torch.zeros((1,4,42)); x[:,:,0]=1; x[:,0,40]=1
    a=torch.eye(4).unsqueeze(0)
    padded=torch.zeros((1,7,42)); padded[:,:4]=x
    ap=torch.zeros((1,7,7)); ap[:,:4,:4]=a
    model=ResidueGNN()
    assert torch.allclose(model(x,a),model(padded,ap),atol=1e-6)


def test_invalid_mutant_is_rejected_without_index_error(tmp_path):
    import pandas as pd
    from protein_rescue.features import graph_for_row
    row=pd.Series({'name':'rcsb_1ABC_A_A1Z_7_25','wt_seq':'ACDE','mut_seq':'ZCDE','ddG':1.})
    assert graph_for_row(row,tmp_path) is None


def test_mutation_site_from_only_sequence_difference_with_pdb_number_check():
    import pandas as pd
    from pathlib import Path
    from protein_rescue.features import graph_for_row
    root=Path(__file__).resolve().parents[1]
    if not (root/'data/raw/1EKG.pdb').exists():return
    df=pd.read_csv(root/'data/raw/S669.csv')
    row=df[df.name=='rcsb_1EKG_A_L198C_7_25'].iloc[0]
    graph=graph_for_row(row,root/'data/raw')
    assert graph is not None
    assert graph['name']=='rcsb_1EKG_A_L198C_7_25'
