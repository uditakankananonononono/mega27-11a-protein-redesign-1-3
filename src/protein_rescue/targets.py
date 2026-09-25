"""Verified canonical-to-PDB residue mapping for disease-target proteins."""
from pathlib import Path
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from .features import AA

TARGETS={
    'TP53':dict(uniprot='P04637',pdb='1TSR',chain='A'),
    'SOD1':dict(uniprot='P00441',pdb='2C9V',chain='A'),
    'PTEN':dict(uniprot='P60484',pdb='1D5R',chain='A'),
}


def map_site(gene:str,canonical_position:int,structures:Path=Path('data/raw'),sequences:Path=Path('data/targets')):
    if gene not in TARGETS: raise ValueError(f'Unknown target: {gene}; available: {", ".join(TARGETS)}')
    t=TARGETS[gene]
    fasta=sequences/f"{t['uniprot']}.fasta"
    sequence=''.join(fasta.read_text().splitlines()[1:])
    if not 1<=canonical_position<=len(sequence):raise ValueError(f'Canonical position {canonical_position} outside 1..{len(sequence)}')
    model=PDBParser(QUIET=True).get_structure(t['pdb'],str(structures/f"{t['pdb']}.pdb"))[0]
    residues=[r for r in model[t['chain']] if r.id[0]==' ' and 'CA' in r and seq1(r.resname,custom_map={'UNK':'X'}) in AA]
    observed=''.join(seq1(r.resname) for r in residues)
    aligner=PairwiseAligner();aligner.mode='global';aligner.match_score=2;aligner.mismatch_score=-1;aligner.open_gap_score=-5;aligner.extend_gap_score=-.5
    canonical_indices,structure_indices=aligner.align(sequence,observed)[0].indices
    matches=structure_indices[canonical_indices==canonical_position-1]
    site=int(matches[0]) if len(matches) and matches[0]>=0 else None
    out=dict(gene=gene,uniprot=t['uniprot'],canonical_position=canonical_position,canonical_aa=sequence[canonical_position-1],pdb=t['pdb'],chain=t['chain'],pdb_residue_number=None,pdb_insertion_code=None,pdb_aa=None,status='not_resolved')
    if site is not None:
        residue=residues[site]
        out.update(pdb_residue_number=residue.id[1],pdb_insertion_code=residue.id[2].strip() or None,pdb_aa=observed[site],status='matched' if observed[site]==out['canonical_aa'] else 'mismatch')
    return out
