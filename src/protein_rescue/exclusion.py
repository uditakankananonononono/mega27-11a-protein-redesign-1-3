"""Explain which S669 rows cannot be mapped to an observed residue graph."""
from pathlib import Path
import json
import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from .features import AA, PAT, align_index


def classify(row,structures:Path,cache:dict):
    match=PAT.match(str(row.name))
    if not match:return 'unparsed_identifier'
    pdb,chain,wt,position,mut=match.groups();label_position=int(position);ref=str(row.wt_seq);alt=str(row.mut_seq)
    if not (wt in AA and mut in AA and len(ref)==len(alt)):return 'invalid_single_substitution'
    differing=[i for i,(a,b) in enumerate(zip(ref,alt)) if a!=b]
    if len(differing)!=1:return 'invalid_single_substitution'
    index=differing[0]
    if ref[index]!=wt or alt[index]!=mut:return 'identifier_aa_disagrees_with_sequences'
    file=structures/f'{pdb}.pdb'
    if not file.is_file():return 'pdb_file_missing'
    if pdb not in cache:
        try:cache[pdb]=PDBParser(QUIET=True).get_structure(pdb,str(file))[0]
        except Exception:cache[pdb]=None
    model=cache[pdb]
    if model is None:return 'pdb_parse_failed'
    if chain not in model:return 'chain_missing'
    residues=[r for r in model[chain] if r.id[0]==' ' and 'CA' in r and seq1(r.resname,custom_map={'UNK':'X'}) in AA]
    if not residues:return 'no_resolved_ca'
    sequence=''.join(seq1(r.resname) for r in residues)
    mapped=align_index(ref,sequence,index)
    if mapped is None:return 'mutated_site_unresolved'
    if sequence[mapped]!=wt:return 'mutated_site_identity_mismatch'
    if residues[mapped].id[1]!=label_position:return 'pdb_number_disagrees_with_identifier'
    return 'candidate_graph'


def audit(csv:Path,structures:Path,output:Path):
    df=pd.read_csv(csv);cache={};counts={};examples={}
    for row in df.itertuples():
        reason=classify(row,structures,cache)
        counts[reason]=counts.get(reason,0)+1
        if len(examples.get(reason,[]))<5: examples.setdefault(reason,[]).append(str(row.name))
    result=dict(rows=len(df),categories=counts,examples=examples,
      caveat='candidate_graph means site and sequence checks passed, but the final graph also requires >=4 local C-alpha residues. Original graph loader remains source of truth for accepted count.')
    output.write_text(json.dumps(result,indent=2)+'\n');return result
