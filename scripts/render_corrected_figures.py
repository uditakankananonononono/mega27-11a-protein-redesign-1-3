"""Render corrected full-S669 diagnostics from committed machine-readable results."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('results'); P=Path('paper')
def load(s): return json.loads((R/s).read_text())
def save(name):
    plt.tight_layout(); plt.savefig(P/name, bbox_inches='tight'); plt.close()
b=load('s669_group_holdout.json'); m=b['metrics']; names=['Train mean','Ridge','Residue GNN','GeoDDG-Seq','GeoDDG-3D','DDMut']; keys=['train_mean','ridge','gnn','GeoDDG-Seq_dir','GeoDDG-3D_dir','DDMut_dir']; fig,ax=plt.subplots(figsize=(9,4)); bars=ax.bar(names,[m[k]['mae'] for k in keys],color=['#aaaaaa','#777777','#b85c38','#26796d','#4b79a1','#785eaa']);ax.set(ylabel='MAE (kcal/mol)',title='S669: 66 variants from 24 unseen PDBs'); ax.tick_params(axis='x',rotation=25);ax.bar_label(bars,fmt='%.3f',padding=2);save('figure_benchmark.pdf')
g=load('group_sensitivity.json'); rows=sorted(g['grouped'],key=lambda x:x['n'],reverse=True);fig,ax=plt.subplots(figsize=(10,4));ax.bar([r['pdb'] for r in rows],[r['n'] for r in rows],color='#4b79a1');ax.set(ylabel='Test variants',xlabel='Held-out PDB',title='Group sizes in the fixed test split');ax.tick_params(axis='x',rotation=75);save('figure_groups.pdf')
fig,ax=plt.subplots(figsize=(10,4)); rows=sorted(g['grouped'],key=lambda x:x['delta']);ax.bar([r['pdb'] for r in rows],[r['delta'] for r in rows],color=['#26796d' if r['delta']<0 else '#b85c38' for r in rows]);ax.axhline(0,color='black',lw=.8);ax.set(ylabel='GNN minus GeoDDG-Seq MAE',xlabel='Held-out PDB',title='Descriptive group errors; positive favors comparator');ax.tick_params(axis='x',rotation=75);save('figure_group_delta.pdf')
loo=load('leave_one_pdb_out.json');fig,ax=plt.subplots(figsize=(8,4));ax.scatter([r['n'] for r in loo['groups']],[r['ridge_mae'] for r in loo['groups']],alpha=.65,color='#4b79a1');ax.set_xscale('log');ax.set(xlabel='Variants in held-out PDB (log scale)',ylabel='Ridge MAE',title='94 leave-one-PDB-out ridge evaluations');save('figure_loo.pdf')
strata=load('error_strata.json');fig,ax=plt.subplots(figsize=(8,4));rows=strata['bins'];bars=ax.bar([r['bin'] for r in rows],[r['delta_mae'] for r in rows],color=['#b85c38' if r['delta_mae']>0 else '#26796d' for r in rows]);ax.axhline(0,color='black',lw=.8);ax.bar_label(bars,labels=[f"n={r['n']}" for r in rows]);ax.set(ylabel='GNN minus GeoDDG-Seq MAE',title='Observed-label strata, descriptive only');ax.tick_params(axis='x',rotation=15);save('figure_error_strata.pdf')
sparse=load('sparse_group_sensitivity.json');rows=[r for r in sparse['thresholds'] if r['retained_proteins']];fig,ax=plt.subplots(figsize=(8,4));bars=ax.bar([str(r['min_rows_per_pdb']) for r in rows],[r['row_weighted_delta_mae'] for r in rows],color='#b85c38');ax.bar_label(bars,labels=[f"{r['retained_proteins']} PDBs" for r in rows]);ax.set(xlabel='Minimum variants per PDB',ylabel='GNN minus GeoDDG-Seq MAE',title='Post hoc sparsity check; population changes by threshold');save('figure_sparse_groups.pdf')
