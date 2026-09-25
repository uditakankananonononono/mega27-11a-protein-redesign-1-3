from pathlib import Path
import pytest
from protein_rescue.structure_compare import compare
ROOT=Path(__file__).resolve().parents[1]


def test_observed_y220c_structure_alignment():
    a=ROOT/'data/raw/1TSR.pdb'; b=ROOT/'data/raw/8QWN.pdb'
    if not (a.exists() and b.exists()):return
    x=compare(a,b,'A',220)
    assert (x['wild_aa'],x['mutant_aa'])==('Y','C')
    assert x['fit_ca']>=170
    assert 0<x['fit_rmsd_angstrom']<5
    assert x['site_ca_displacement_angstrom']>=0
