"""Prespecified DDG bins for diagnostic error, not post-hoc performance claims."""
from pathlib import Path
import json
import pandas as pd
import numpy as np


def run(results:Path,source:Path,output:Path):
    report=json.loads(results.read_text());pred=pd.DataFrame(report['predictions'])
    table=pd.read_csv(source)[['name','GeoDDG-Seq_dir']]
    df=pred.merge(table,on='name',validate='one_to_one')
    if len(df)!=len(pred):raise ValueError('Comparator unavailable on some rows')
    df['gnn_error']=(df.gnn-df.experimental_ddg).abs()
    df['comparator_error']=(df['GeoDDG-Seq_dir']-df.experimental_ddg).abs()
    df['delta']=df.gnn_error-df.comparator_error
    bins=[(-float('inf'),-1.5,'ddg <= -1.5'),(-1.5,-.5,'-1.5 < ddg <= -0.5'),(-.5,.5,'-0.5 < ddg <= 0.5'),(.5,float('inf'),'ddg > 0.5')]
    strata=[]
    for low,high,label in bins:
        sample=df[(df.experimental_ddg>low)&(df.experimental_ddg<=high)]
        strata.append(dict(bin=label,n=int(len(sample)),pdb_groups=int(sample.pdb.nunique()),gnn_mae=float(sample.gnn_error.mean()),comparator_mae=float(sample.comparator_error.mean()),delta_mae=float(sample.delta.mean())))
    result=dict(source_results=str(results),comparator='GeoDDG-Seq_dir',bins=strata,n_test=len(df),
      gnn_better_rows=int((df.delta<0).sum()),comparator_better_rows=int((df.delta>0).sum()),
      gnn_better_pdb_groups=int((df.groupby('pdb').delta.mean()<0).sum()),comparator_better_pdb_groups=int((df.groupby('pdb').delta.mean()>0).sum()),
      caveats=['DDG bins are exploratory descriptors based on observed test labels; do not use them as a prospective gate without calibration.','The high-positive bin has only 11 variants from a small number of proteins; group dependence and label imbalance preclude a subgroup discovery claim.'])
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');return result
