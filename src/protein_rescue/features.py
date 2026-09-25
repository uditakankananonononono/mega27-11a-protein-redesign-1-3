"""Experimental mutation rows to aligned residue contact graphs; no network in this module."""
from __future__ import annotations
import re
from pathlib import Path
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

AA = 'ACDEFGHIKLMNPQRSTVWY'
PAT = re.compile(r'^rcsb_([A-Za-z0-9]{4})_([^_]+)_([A-Z])(\d+)([A-Z])_')


def align_index(reference: str, observed: str, ref_index: int) -> int | None:
    """Map zero-based reference index to observed index using global sequence alignment."""
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    a, b = aligner.align(reference, observed)[0].indices
    hits = b[a == ref_index]
    return int(hits[0]) if len(hits) and hits[0] >= 0 else None


def graph_for_row(row: pd.Series, structures: Path, radius: float = 12., nodes: int = 32):
    m = PAT.match(str(row['name']))
    if not m:
        return None
    pdb, chain, wt, index, mutant = m.groups()
    index = int(index) - 1
    ref = str(row['wt_seq'])
    alt = str(row['mut_seq'])
    if not (wt in AA and mutant in AA and 0 <= index < len(ref) == len(alt) and ref[index] == wt and alt[index] == mutant
            and sum(a != b for a, b in zip(ref, alt)) == 1):
        return None
    path = structures / f'{pdb}.pdb'
    if not path.is_file():
        return None
    model = PDBParser(QUIET=True).get_structure(pdb, str(path))[0]
    if chain not in model:
        return None
    residues = [r for r in model[chain] if r.id[0] == ' ' and 'CA' in r and seq1(r.resname, custom_map={'UNK': 'X'}) in AA]
    if not residues:
        return None
    sequence = ''.join(seq1(r.resname) for r in residues)
    mapped = align_index(ref, sequence, index)
    if mapped is None or sequence[mapped] != wt:
        return None
    xyz = np.array([r['CA'].coord for r in residues], dtype=np.float32)
    distances = np.linalg.norm(xyz - xyz[mapped], axis=1)
    chosen = np.flatnonzero(distances <= radius)
    chosen = chosen[np.argsort(distances[chosen])[:nodes]]
    if len(chosen) < 4:
        return None
    local_xyz = xyz[chosen]
    pairdist = np.linalg.norm(local_xyz[:, None] - local_xyz[None, :], axis=-1)
    adjacency = (pairdist <= 8.).astype(np.float32)
    np.fill_diagonal(adjacency, 1.)
    inv = 1 / np.maximum(adjacency.sum(1), 1)
    adjacency = (np.sqrt(inv)[:, None] * adjacency * np.sqrt(inv)[None, :]).astype(np.float32)
    aa = np.eye(20, dtype=np.float32)[[AA.index(sequence[i]) for i in chosen]]
    site = (chosen == mapped).astype(np.float32)[:, None]
    radial = (distances[chosen] / radius).astype(np.float32)[:, None]
    mutant_context = np.zeros((len(chosen), 20), dtype=np.float32)
    mutant_context[chosen == mapped, AA.index(mutant)] = 1.
    x = np.concatenate([aa, mutant_context, site, radial], axis=1)
    return dict(x=x, a=adjacency, y=float(row['ddG']), name=str(row['name']), pdb=pdb)


def load_graphs(csv_path: Path, structure_dir: Path):
    data = pd.read_csv(csv_path)
    graphs = []
    reasons = dict(total=len(data), valid=0, unavailable=0)
    for _, row in data.iterrows():
        g = graph_for_row(row, structure_dir)
        if g is None:
            reasons['unavailable'] += 1
        else:
            graphs.append(g)
            reasons['valid'] += 1
    return graphs, reasons


def tabular_features(g):
    site = int(np.argmax(g['x'][:, 40]))
    aa = g['x'][:, :20]
    return np.r_[aa[site], g['x'][site, 20:40], aa.mean(0), g['x'][:, 41].mean(), len(aa)]
