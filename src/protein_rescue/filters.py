"""Predictor abstention rules based on measured out-of-domain risk, not confidence theater."""
from pathlib import Path
import json
import numpy as np


def safe_to_rank(variant:dict,minimum_proteins:int=10,minimum_rows:int=20,require_target_supported:bool=True):
    problems=[]
    if variant.get('pdb_accession') is None:problems.append('Missing verified PDB accession')
    if variant.get('sequence_aligned') is not True:problems.append('Canonical and PDB sequence not verified')
    if variant.get('numbering_match') is not True:problems.append('Mutation numbering not matched to PDB')
    if variant.get('mutant_assay_ddg_sign_known') is not True:problems.append('Assay DDG convention not verified')
    if int(variant.get('relevant_training_proteins') or 0) < minimum_proteins:problems.append('Too few relevant training proteins')
    if int(variant.get('relevant_training_rows') or 0) < minimum_rows:problems.append('Too few relevant training mutation rows')
    if require_target_supported and variant.get('target_double_mutant_assay_supported') is not True:
        problems.append('No target-specific double-mutant assay validation')
    return dict(eligible=not problems,reasons=problems,
                note='This checklist is an abstention guard, not a probability of rescue or model calibration.')
