from pathlib import Path
from protein_rescue.targets import map_site
ROOT=Path(__file__).resolve().parents[1]


def test_target_numbering_real_structures():
    if not (ROOT/'data/raw/2C9V.pdb').exists():return
    result=map_site('SOD1',5,ROOT/'data/raw',ROOT/'data/targets')
    assert (result['canonical_aa'],result['pdb_residue_number'],result['pdb_aa'],result['status'])==('A',4,'A','matched')
    assert map_site('SOD1',94,ROOT/'data/raw',ROOT/'data/targets')['pdb_residue_number']==93
    assert map_site('TP53',220,ROOT/'data/raw',ROOT/'data/targets')['pdb_residue_number']==220
    assert map_site('PTEN',130,ROOT/'data/raw',ROOT/'data/targets')['pdb_residue_number']==130
