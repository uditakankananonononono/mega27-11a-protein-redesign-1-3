"""Audit whether PDB-disjoint split still admits identical reference sequences."""
from pathlib import Path
import json
import pandas as pd
from collections import defaultdict


def audit(csv:Path,benchmark:Path,output:Path):
    df=pd.read_csv(csv)
    data=json.loads(benchmark.read_text())
    used=set(data['split']['train_pdbs'])|set(data['split']['test_pdbs'])
    df['pdb']=df.name.str.extract(r'^rcsb_([A-Za-z0-9]{4})_')[0]
    ref=df[df.pdb.isin(used)][['pdb','wt_seq']].drop_duplicates()
    grouped=defaultdict(set)
    for row in ref.itertuples():grouped[row.wt_seq].add(row.pdb)
    train=set(data['split']['train_pdbs']);test=set(data['split']['test_pdbs'])
    duplicates=[dict(sequence_length=len(seq),pdbs=sorted(ids),train=sorted(ids&train),test=sorted(ids&test)) for seq,ids in grouped.items() if len(ids)>1]
    crossing=[d for d in duplicates if d['train'] and d['test']]
    result=dict(unique_pdbs=len(used),unique_reference_sequences=len(grouped),duplicated_sequences=len(duplicates),
                duplicated_across_train_test=len(crossing),crossing=crossing,all_duplicates=duplicates,
                caveat='Exact reference-sequence identity only; high homology, similar folds, homologous domains or overlap with published predictor training sets remain unassessed.')
    output.write_text(json.dumps(result,indent=2)+'\n');return result
