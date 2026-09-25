"""Sample-size sensitivity: does the sign of per-protein advantage survive sparse-group removal?"""
from pathlib import Path
import json
import pandas as pd
import numpy as np


def run(audit:Path,output:Path):
    x=json.loads(audit.read_text());g=pd.DataFrame(x['grouped'])
    thresholds=[]
    for min_size in [1,2,3,4,5,10,20,40]:
        subset=g[g.n>=min_size]
        if subset.empty:continue
        thresholds.append(dict(min_rows_per_pdb=min_size,retained_proteins=int(len(subset)),retained_variants=int(subset.n.sum()),equal_protein_delta_mae=float(subset.delta.mean()),row_weighted_delta_mae=float(np.average(subset.delta,weights=subset.n)),gnn_better_proteins=int((subset.delta<0).sum()),comparator_better_proteins=int((subset.delta>0).sum())))
    result=dict(source=str(audit),direction='positive: GNN higher MAE than comparator',thresholds=thresholds,
                caveat='Post-hoc filtering on test-group size changes the target population. This is sensitivity analysis, not an optimized evaluation protocol or an independent test.')
    output.write_text(json.dumps(result,indent=2)+'\n');return result
