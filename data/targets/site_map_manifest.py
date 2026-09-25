"""Regenerate disease-site canonical-to-structure map without network access."""
from pathlib import Path
import json
from protein_rescue.targets import map_site


def main():
    sites={'TP53':[123,143,168,175,220,235,239,240,245,248,249,268,273],
           'SOD1':[5,91,94], 'PTEN':[107,130,335]}
    records=[map_site(g,p) for g,positions in sites.items() for p in positions]
    out=Path('results/disease_site_numbering.json')
    out.write_text(json.dumps(records,indent=2)+'\n')
    print(f'wrote {len(records)} canonical sites: {sum(r["status"]=="matched" for r in records)} matched')


if __name__=='__main__':main()
