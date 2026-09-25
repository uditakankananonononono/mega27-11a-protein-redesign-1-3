# Protein rescue, item 11a

Work in progress. The dataset, service, formula and page-count gates are audited, but a validated rescue finding remains open. No rescue mutation is experimentally validated by this repository. The S669 table contains experimental stability changes and third-party predictions; its aggregate performance must never be represented as a held-out test for a predictor trained on S669. Evaluation groups by PDB identifier.

Source data: [S669](https://github.com/Gonglab-THU/GeoStab/blob/main/data/ddG/S669.csv), original benchmark [Zenodo](https://zenodo.org/records/7568094), structures [RCSB](https://www.rcsb.org/structure/1TSR). Data is tracked with source attribution; measured findings must be read from `results/` after reproducible runs.

## Run

```bash
pip install -e .
python data/fetch_structures.py
protein-rescue benchmark --epochs 15 --graph-cache results/checkpoints  # only if independently regenerated from the exact CSV, PDB files and checked mapping code
# without the cache: protein-rescue benchmark --epochs 15 (may exceed a short shell timeout)
protein-rescue audit
python data/prepare_fireprot.py --source data/raw/fireprotdb.csv
protein-rescue fireprot-audit
protein-rescue inspect rcsb_1A0F_A_S11A_6_56
protein-rescue map-site SOD1 5  # UniProt A5 -> mature PDB 2C9V A4
protein-rescue cycle 1.69 2.75 1.72  # experimental R249S/H168R interaction; not prediction
protein-rescue rank-check --context results/disease_target_abstention.json
pytest -q
```

The raw FireProtDB bulk export used for the present audit timed out and is a partial snapshot, so the observed subset is not the full database. For a complete API export, page requests by offset and check both expected count and end-of-data. The raw FireProtDB bulk export is downloaded manually from its [official export](https://loschmidt.chemi.muni.cz/fireprotdb/download/) or by GET of its documented `/api/search?format=csv&sort=` route; it is not tracked because its size exceeds 300 MB. Downloaded PDB files are ignored in Git and their source addresses and checksums are recorded in `data/structure_fetch_manifest.json` for the initial S669 panel. FireProtDB PDB fetches are reproducible from identifiers in `data/fireprot_pdb_ddg.csv` using the RCSB endpoint. Inspect the source license and the exact model/data purpose before redistribution.

## Status and limits

The corrected exploratory manuscript is in `paper/manuscript.pdf`: 20 rendered pages, 18 displayed equations, nine figures and full-data audit tables. It now uses embedded original Times New Roman for body and headings; scientific equations retain their specialized math fonts. The only bundled artifact is the PDF, not the licensed TTF files. The font installer and EULA remain in the program Drive folder. The 669-row S669 mapping yields 94 PDBs; a PDB-group holdout has 603 training and 66 test variants across 70 and 24 disjoint PDBs. On that test split, GNN MAE is 1.422 versus stored GeoDDG-Seq 1.310 and DDMut 1.290 kcal/mol, so this GNN did not outperform those published predictions. Comparator historical training overlap remains unknown. The 163 accession-backed PDBs were fetched and used: 94 S669, 74 FireProtDB direct mutation-number matches (11 overlap), and six mutant-structure comparisons. The 40 distinct scientific-service ledger excludes GitHub/Drive and collapses Ensembl endpoints; service count is not a study count or a measure of evidentiary strength. The repo is published privately to https://github.com/uditakankananonononono/mega27-11a-protein-redesign-1-3, with backup history in the program Drive folder. The old 420-row result and PDF are explicitly quarantined in `results/legacy_420_subset` and `paper/legacy_420_subset` because the first site mapper confused PDB label numbering with sequence offsets. They are not the current benchmark.

FireProtDB DDG records were audited but not pooled with S669 because assay definitions, conditions and signs need harmonization. The CLI supplies a benchmark, data audit, residue mapper and experimental double-mutant cycle calculator, not clinically valid mutation recommendations. No therapeutic rescue mutation is claimed.
