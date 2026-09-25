"""Build a compact provenance-preserving sample of PDB-linked FireProtDB DDG records."""
import argparse
from pathlib import Path
import pandas as pd

COLS=['EXPERIMENT_ID','SUBSTITUTION','DDG','WWPDB','UNIPROTKB','SOURCE_DATASET','PH','EXP_TEMPERATURE','SEQUENCE_LENGTH']

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,default=Path('data/raw/fireprotdb.csv'))
    p.add_argument('--output',type=Path,default=Path('data/fireprot_pdb_ddg.csv'))
    a=p.parse_args()
    kept=[]; scanned=0
    for chunk in pd.read_csv(a.source,usecols=COLS,dtype=str,chunksize=10000,on_bad_lines='skip'):
        scanned+=len(chunk)
        chunk=chunk.dropna(subset=['WWPDB','SUBSTITUTION','DDG'])
        chunk=chunk[chunk.WWPDB.str.fullmatch('[A-Za-z0-9]{4}',na=False) & chunk.SUBSTITUTION.str.fullmatch('[A-Z][0-9]+[A-Z]',na=False)]
        chunk['DDG']=pd.to_numeric(chunk['DDG'],errors='coerce')
        chunk=chunk.dropna(subset=['DDG'])
        kept.append(chunk)
    out=pd.concat(kept,ignore_index=True)
    out.to_csv(a.output,index=False)
    print(f'scanned={scanned}, filtered={len(out)}, distinct_PDBs={out.WWPDB.nunique()}, output={a.output}')

if __name__=='__main__': main()
