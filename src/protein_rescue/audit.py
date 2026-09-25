"""Group-aware error analysis on stored predictions. Does not train or relabel."""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd


def audit(results: Path, csv: Path, output: Path, samples: int = 5000, seed: int = 27):
    report=json.loads(results.read_text())
    predictions=pd.DataFrame(report['predictions'])
    labels=pd.read_csv(csv)[['name', 'GeoDDG-Seq_dir', 'GeoDDG-3D_dir']]
    p=predictions.merge(labels,on='name',validate='one_to_one')
    if len(p)!=len(predictions): raise ValueError('Missing comparator predictions')
    p['gnn_abs_error']=(p.gnn-p.experimental_ddg).abs()
    p['geo_abs_error']=(p['GeoDDG-Seq_dir']-p.experimental_ddg).abs()
    p['delta']=p.gnn_abs_error-p.geo_abs_error
    grouped=p.groupby('pdb').agg(n=('name','count'),gnn_mae=('gnn_abs_error','mean'),geo_mae=('geo_abs_error','mean'),delta=('delta','mean')).reset_index()
    rng=np.random.default_rng(seed)
    groups=grouped.pdb.to_numpy()
    group_rows={k:p[p.pdb==k] for k in groups}
    row_boot=[]
    group_boot=[]
    for _ in range(samples):
        row_boot.append(float(rng.choice(p.delta.to_numpy(),len(p),replace=True).mean()))
        chosen=rng.choice(groups,len(groups),replace=True)
        group_boot.append(float(np.mean([group_rows[g].delta.mean() for g in chosen])))
    q=lambda x:[float(z) for z in np.quantile(x,[.025,.5,.975])]
    result=dict(source_results=str(results), comparator='GeoDDG-Seq_dir', definition='positive delta: GNN absolute error exceeds GeoDDG-Seq absolute error',
        rows=len(p),protein_groups=len(groups),largest_two_row_fraction=float(grouped.nlargest(2,'n').n.sum()/len(p)),
        row_weighted_delta_mae=float(p.delta.mean()),protein_equal_delta_mae=float(grouped.delta.mean()),
        bootstrap_row_95=q(row_boot),bootstrap_protein_equal_95=q(group_boot),
        grouped=grouped.to_dict('records'),seed=seed,bootstrap_samples=samples,
        caveat='Bootstrap conditions on this single fixed held-out split; group count is small. Existing comparator scores may have unknown training overlap with these proteins; no head-to-head train-data equivalence is claimed.')
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result
