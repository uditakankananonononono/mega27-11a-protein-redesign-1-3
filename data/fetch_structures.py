"""Fetch immutable, accession-addressed public PDB inputs referenced by S669."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.request import urlopen
import hashlib
import json
import pandas as pd


def fetch_one(pdb, target):
    url=f'https://files.rcsb.org/download/{pdb}.pdb'
    try:
        data=urlopen(url,timeout=20).read()
        if b'ATOM  ' not in data:
            raise ValueError('No ATOM record')
        (target / f'{pdb}.pdb').write_bytes(data)
        return dict(accession=pdb,url=url,bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),status='fetched')
    except Exception as exc:
        return dict(accession=pdb,url=url,status='failed',error=repr(exc))


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--csv',type=Path,default=Path('data/raw/S669.csv'))
    p.add_argument('--dest',type=Path,default=Path('data/raw'))
    p.add_argument('--manifest',type=Path,default=Path('data/structure_fetch_manifest.json'))
    args=p.parse_args()
    ids=sorted(set(pd.read_csv(args.csv).name.str.extract(r'^rcsb_([A-Za-z0-9]{4})_')[0].dropna()) | {'1TSR'})
    args.dest.mkdir(parents=True,exist_ok=True)
    results=[]
    with ThreadPoolExecutor(max_workers=8) as pool:
        jobs={pool.submit(fetch_one,id,args.dest):id for id in ids if not (args.dest / f'{id}.pdb').exists()}
        for job in as_completed(jobs):
            results.append(job.result())
    for id in ids:
        path=args.dest / f'{id}.pdb'
        if path.is_file() and id not in {r['accession'] for r in results}:
            data=path.read_bytes()
            results.append(dict(accession=id,url=f'https://files.rcsb.org/download/{id}.pdb',bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),status='present'))
    args.manifest.write_text(json.dumps(sorted(results,key=lambda r:r['accession']),indent=2)+'\n')
    print('Structures:',len(ids),'local:',sum(r['status']!='failed' for r in results),'failed:',sum(r['status']=='failed' for r in results))


if __name__=='__main__': main()
