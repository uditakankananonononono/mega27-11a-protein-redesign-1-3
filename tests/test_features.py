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
