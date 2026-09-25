"""Fast leave-one-PDB-out ridge and train-mean baselines, one row per protein."""
from pathlib import Path
import json
import numpy as np
from scipy.stats import spearmanr
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from .features import load_graphs, tabular_features


def run(data:Path,structures:Path,output:Path):
    graphs,audit=load_graphs(data,structures)
    x=np.array([tabular_features(g) for g in graphs]);y=np.array([g['y'] for g in graphs]);groups=np.array([g['pdb'] for g in graphs]);scores=[]
    for group in np.unique(groups):
        train=groups!=group;test=~train
        pred=make_pipeline(StandardScaler(),Ridge(alpha=100.)).fit(x[train],y[train]).predict(x[test])
        baseline=np.full(test.sum(),y[train].mean())
        scores.append(dict(pdb=group,n=int(test.sum()),ridge_mae=float(np.abs(pred-y[test]).mean()),train_mean_mae=float(np.abs(baseline-y[test]).mean()),ridge_spearman=float(spearmanr(y[test],pred).statistic) if test.sum()>2 and np.std(y[test])>0 else None))
    result=dict(audit=audit,n_groups=len(scores),groups=scores,
      row_weighted_ridge_mae=float(np.average([r['ridge_mae'] for r in scores],weights=[r['n'] for r in scores])),
      protein_equal_ridge_mae=float(np.mean([r['ridge_mae'] for r in scores])),
      row_weighted_train_mean_mae=float(np.average([r['train_mean_mae'] for r in scores],weights=[r['n'] for r in scores])),
      protein_equal_train_mean_mae=float(np.mean([r['train_mean_mae'] for r in scores])),
      caveat='Leave-one-structure-out within a heterogeneous 52-PDB subset; protein family homology is not grouped, and hyperparameter alpha=100 was fixed a priori, not tuned on held-out rows.')
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');return result
