"""A cautious uncertainty audit of repeated experimental database measurements."""
from pathlib import Path
import json
import pandas as pd
import numpy as np


def run(source:Path,output:Path):
    df=pd.read_csv(source)
    groups=[]
    for (pdb,mutation),g in df.groupby(['WWPDB','SUBSTITUTION']):
        if len(g)<2: continue
        spans=g.groupby(['PH','EXP_TEMPERATURE'],dropna=False).DDG.agg(['count','min','max'])
        groups.append(dict(pdb=pdb,mutation=mutation,n=int(len(g)),span=float(g.DDG.max()-g.DDG.min()),
            same_ph_temp_max_span=float((spans['max']-spans['min']).max()),
            sign_discordant=bool(g.DDG.min()<0<g.DDG.max()),
            same_recorded_ph_temp_sign_discordant=bool(((spans['min']<0)&(spans['max']>0)).any())))
    repeated=pd.DataFrame(groups)
    result=dict(total_filtered_experiments=len(df),repeated_pdb_mutation_groups=len(repeated),
        repeated_with_any_nonzero_span=int((repeated.span>0).sum()),
        span_quantiles={str(q):float(repeated.span.quantile(q)) for q in [.25,.5,.75,.9,.95]},
        repeated_sign_discordant=int(repeated.sign_discordant.sum()),
        repeated_same_recorded_ph_temp_sign_discordant=int(repeated.same_recorded_ph_temp_sign_discordant.sum()),
        repeated_group_detail=groups,
        caveat='PDB ID plus substitution is not a unique experiment; constructs, assay methods, buffers, data-source definitions and unrecorded conditions can differ. Quantiles are variation among database rows, not instrument error or predictive-model uncertainty.')
    output.write_text(json.dumps(result,indent=2)+'\n');return result
