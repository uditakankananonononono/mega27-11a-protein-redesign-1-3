# Protein rescue, item 11a

Work in progress. No rescue mutation is experimentally validated by this repository. The S669 table contains experimental stability changes and third-party predictions; its aggregate performance must never be represented as a held-out test for a predictor trained on S669. Evaluation groups by PDB identifier.

Source data: [S669](https://github.com/Gonglab-THU/GeoStab/blob/main/data/ddG/S669.csv), original benchmark [Zenodo](https://zenodo.org/records/7568094), structures [RCSB](https://www.rcsb.org/structure/1TSR). Data is tracked with source attribution; measured findings must be read from `results/` after reproducible runs.

## Run

```bash
pip install -e .
python data/fetch_structures.py
protein-rescue benchmark --epochs 15
protein-rescue audit
python data/prepare_fireprot.py --source data/raw/fireprotdb.csv
protein-rescue fireprot-audit
protein-rescue inspect rcsb_1A0F_A_S11A_6_56
pytest -q
```

The raw FireProtDB bulk export used for the present audit timed out and is a partial snapshot, so the observed subset is not the full database. For a complete API export, page requests by offset and check both expected count and end-of-data. The raw FireProtDB bulk export is downloaded manually from its [official export](https://loschmidt.chemi.muni.cz/fireprotdb/download/) or by GET of its documented `/api/search?format=csv&sort=` route; it is not tracked because its size exceeds 300 MB. Downloaded PDB files are ignored in Git and their source addresses and checksums are recorded in `data/structure_fetch_manifest.json` for the initial S669 panel. FireProtDB PDB fetches are reproducible from identifiers in `data/fireprot_pdb_ddg.csv` using the RCSB endpoint. Inspect the source license and the exact model/data purpose before redistribution.

## Status and limits

The exploratory manuscript is in `paper/manuscript.pdf` (rendered six pages, 14 numbered equations). This is a draft below the requested 50-page paper gate. The GNN loses to published stored S669 predictions. FireProtDB DDG records were audited but not trained together with S669 due to uncertain sign and condition conventions. The CLI supplies a benchmark and data audit, not clinically valid mutation recommendations. See `results/` JSON for measured results and failures. No therapeutic rescue mutation is claimed.
