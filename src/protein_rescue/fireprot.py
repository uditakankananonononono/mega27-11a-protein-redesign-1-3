"""Inspect contradictions and map FireProtDB PDB structures without crossing DDG sign conventions."""
from pathlib import Path
import json
import re
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from .features import AA


def analyze(path:Path, structures:Path, output:Path):
    df=pd.read_csv(path)
    per=[]; valid=set(); chains_checked=0
    for pdb,group in df.groupby('WWPDB'):
        path_pdb=structures/f'{pdb}.pdb'
        if not path_pdb.exists():continue
        try:
            model=PDBParser(QUIET=True).get_structure(pdb,str(path_pdb))[0]
        except Exception:continue
        chain_maps=[]
        for c in model:
            m={r.id[1]:seq1(r.resname,custom_map={'UNK':'X'}) for r in c if r.id[0]==' ' and 'CA' in r}
            chain_maps.append(m);chains_checked+=1
        direct=0
        for row in group.itertuples():
            mt=re.fullmatch(r'([A-Z])(\d+)([A-Z])',row.SUBSTITUTION)
            if mt and any(m.get(int(mt[2]))==mt[1] for m in chain_maps):
                direct+=1;valid.add(str(row.EXPERIMENT_ID))
        per.append(dict(pdb=pdb,rows=len(group),direct_numbering_matches=direct,chains=len(chain_maps)))
    df['experiment']=df.EXPERIMENT_ID.astype(str)
    contradiction=[]
    same_conditions_discordant=0
    for (pdb,mutation),g in df.groupby(['WWPDB','SUBSTITUTION']):
        if len(g)<2:continue
        neg=g[g.DDG < 0];pos=g[g.DDG > 0]
        if len(neg) and len(pos):
            if any((h.DDG.min()<0 and h.DDG.max()>0) for _,h in g.groupby(['PH','EXP_TEMPERATURE'],dropna=False)):
                same_conditions_discordant += 1
            contradiction.append(dict(pdb=pdb,mutation=mutation,experiments=len(g),min_ddg=float(g.DDG.min()),max_ddg=float(g.DDG.max()),source_experiment_ids=g.EXPERIMENT_ID.astype(str).tolist()))
    result=dict(source=str(path),rows=len(df),distinct_pdbs=int(df.WWPDB.nunique()),structures_available=len(per),structures_with_direct_matches=sum(p['direct_numbering_matches']>0 for p in per),direct_numbering_experiments=len(valid),chains_checked=chains_checked,
        opposite_sign_same_pdb_mutation_pairs=len(contradiction),same_recorded_ph_and_temperature_discordant_pairs=same_conditions_discordant,opposite_sign_examples=contradiction[:20],per_pdb=per,
        caveats=['FireProtDB DDG convention differs from S669 and cannot be pooled without explicit sign calibration.','Repeated rows may differ in measurement conditions or source proteins; opposite signs do not by themselves prove contradictory experimental findings.','PDB numbering must match independently; absent matching residue numbers require sequence alignment, not assumed conversion.'])
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n')
    return result
