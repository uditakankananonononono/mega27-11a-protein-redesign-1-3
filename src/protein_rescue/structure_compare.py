"""Compare observed mutant and wild-type backbones by common numbered C-alpha residues.

A coordinate comparison is not a folding free-energy or functional rescue prediction.
"""
from pathlib import Path
import json
import numpy as np
from Bio.PDB import PDBParser, Superimposer
from Bio.SeqUtils import seq1


def compare(wild:Path,mutant:Path,chain:str,mutation_position:int):
    parser=PDBParser(QUIET=True)
    w=parser.get_structure(wild.stem,str(wild))[0][chain]
    m=parser.get_structure(mutant.stem,str(mutant))[0][chain]
    common=sorted(set(r.id for r in w if r.id[0]==' ' and 'CA' in r) & set(r.id for r in m if r.id[0]==' ' and 'CA' in r))
    if len(common)<20:raise ValueError('Fewer than 20 identically numbered C-alpha residues')
    matched=[r for r in common if seq1(w[r].resname,custom_map={'UNK':'X'})==seq1(m[r].resname,custom_map={'UNK':'X'})]
    if len(matched)<20:raise ValueError('Fewer than 20 identity-matched common C-alpha residues')
    sup=Superimposer();sup.set_atoms([w[r]['CA'] for r in matched],[m[r]['CA'] for r in matched])
    rotation,translation=sup.rotran
    differences={str(r[1]):float(np.linalg.norm(w[r]['CA'].coord - (m[r]['CA'].coord@rotation+translation))) for r in common}
    nearby=[r for r in common if r[1]==mutation_position]
    return dict(wild_pdb=wild.stem,mutant_pdb=mutant.stem,chain=chain,mutation_position=mutation_position,
      wild_aa=seq1(w[nearby[0]].resname) if nearby else None,mutant_aa=seq1(m[nearby[0]].resname) if nearby else None,
      common_ca=len(common),fit_ca=len(matched),fit_rmsd_angstrom=float(sup.rms),site_ca_displacement_angstrom=differences.get(str(mutation_position)),
      median_ca_displacement_angstrom=float(np.median(list(differences.values()))),
      per_residue_ca_displacement_angstrom=differences,
      caveat='Static crystal backbone coordinates with matched residue numbering; sample conditions, ligands and crystal packing may differ. This is not a mutant free energy, binding, epistasis or functional rescue measurement.')


def panel(output:Path,structures:Path):
    pairs=[('1TSR','8QWN','A',220),('1TSR','6SHZ','A',220),('1TSR','4AGO','A',220),('2C9V','1UXM','A',4),('2C9V','3GZO','A',93),('2C9V','3GZQ','A',4)]
    result=[compare(structures/f'{w}.pdb',structures/f'{m}.pdb',c,pos) for w,m,c,pos in pairs]
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result
