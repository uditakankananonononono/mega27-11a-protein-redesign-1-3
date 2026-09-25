"""Transparent double-mutant-cycle calculator for experimental DDG values.

This is a calculation on measured inputs, not a predictor of unmeasured rescue.
"""
from __future__ import annotations
import math
from pathlib import Path
import json


def double_mutant_cycle(single_a: float, single_b: float, combined: float,
                        se_a: float | None = None, se_b: float | None = None, se_combined: float | None = None):
    vals=(single_a,single_b,combined)
    if any(not math.isfinite(v) for v in vals):raise ValueError('All DDG inputs must be finite')
    values=(se_a,se_b,se_combined)
    if any(v is not None and (not math.isfinite(v) or v<0) for v in values):raise ValueError('Uncertainty values must be finite, nonnegative')
    if any(v is not None for v in values) and any(v is None for v in values):raise ValueError('Supply all three independent SEs or none')
    predicted=single_a+single_b
    result=dict(additive_ddg=predicted,observed_ddg=combined,interaction_ddg=combined-predicted,
      interpretation='A nonzero interaction is a nonadditive difference in the measurement convention supplied; it does not establish functional rescue or causal contact.')
    if all(v is not None for v in values):
        result['interaction_se_independence_assumption']=math.sqrt(sum(v*v for v in values))
    return result


def literature_controls(source:Path,output:Path):
    x=json.loads(source.read_text())
    d={m['variant']:m for m in x['mutants']}
    pairings=[('G245S','N239Y','G245S/N239Y'),('R249S','H168R','R249S/H168R')]
    result=dict(source=x['provenance'],models='Experimental inputs from source, no GNN inference',controls={})
    for a,b,ab in pairings:
        result['controls'][ab]=double_mutant_cycle(d[a]['ddg'],d[b]['ddg'],d[ab]['ddg'],d[a]['se'],d[b]['se'],d[ab]['se'])
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n')
    return result
