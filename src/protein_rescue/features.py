"""Experimental mutation rows to aligned residue contact graphs; no network in this module."""
from __future__ import annotations
import re
from functools import lru_cache
from pathlib import Path
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

AA = 'ACDEFGHIKLMNPQRSTVWY'
PAT = re.compile(r'^rcsb_([A-Za-z0-9]{4})_([^_]+)_([A-Z])(\d+)([A-Z])_')


@lru_cache(maxsize=4096)
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


@lru_cache(maxsize=512)
def _cached_structure(path: str):
    return PDBParser(QUIET=True).get_structure(Path(path).stem, path)[0]


def graph_for_row(row: pd.Series, structures: Path, radius: float = 12., nodes: int = 32):
    m = PAT.match(str(row['name']))
    if not m:
        return None
    pdb, chain, wt, index, mutant = m.groups()
    label_position = int(index)
    ref = str(row['wt_seq'])
    alt = str(row['mut_seq'])
    if wt not in AA or mutant not in AA or len(ref) != len(alt):
        return None
    differing = [i for i, (a,b) in enumerate(zip(ref, alt)) if a != b]
    if len(differing) != 1:
        return None
    index = differing[0]
    if ref[index] != wt or alt[index] != mutant:
        return None
    # Name positions may use construct/PDB numbering rather than sequence offset.
    # Require that resolved PDB residue numbering agrees with the mutation label.
    path = structures / f'{pdb}.pdb'
    if not path.is_file():
        return None
    model = _cached_structure(str(path))
    if chain not in model:
        return None
    residues = [r for r in model[chain] if r.id[0] == ' ' and 'CA' in r and seq1(r.resname, custom_map={'UNK': 'X'}) in AA]
    if not residues:
        return None
    sequence = ''.join(seq1(r.resname) for r in residues)
    mapped = align_index(ref, sequence, index)
    if mapped is None or sequence[mapped] != wt or residues[mapped].id[1] != label_position:
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
