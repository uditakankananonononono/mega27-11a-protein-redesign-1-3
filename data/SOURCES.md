GeoStab S669.csv: https://github.com/Gonglab-THU/GeoStab/blob/main/data/ddG/S669.csv ; experimental ddG labels and published comparator predictions; retrieval September 25 2026.
PDB structures: https://files.rcsb.org/download/{PDB_ID}.pdb ; use only successful unique IDs corresponding to rows of S669, plus disease target 1TSR. A PDB downloaded but not processed is not counted as used.
S669 benchmark publication/data: https://zenodo.org/records/7568094 ; https://doi.org/10.1093/bib/bbab555 .

FireProtDB official CSV export: https://loschmidt.chemi.muni.cz/fireprotdb/download/ ; API https://loschmidt.chemi.muni.cz/fireprotdb/api/search?format=csv&sort= ; documentation for negative-stabilizing DDG https://loschmidt.chemi.muni.cz/fireprotdb/help/ . Export retrieval September 25 2026; source bulk file is not committed. Filtered experimental rows are committed in fireprot_pdb_ddg.csv, and conditions are not silently combined.
The source-of-truth for the 120 processed accessions and their PDB file checksums is results/accession_count.json; this excludes fetched-but-unmapped inputs.
