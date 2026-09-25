import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def test_accession_ledger_unique_with_reproducible_split():
    ledger=json.loads((ROOT/'results/accession_count.json').read_text())
    assert ledger['union']==len(ledger['accessions'])==163
    assert len({r['pdb'] for r in ledger['accessions']})==163
    assert ledger['s669_mutation_mapped_structures']+ledger['fireprot_direct_mutation_mapped_structures']-ledger['overlap']==157
    assert ledger['mutant_structure_comparison_accessions']==6
    assert all(len(r['sha256'])==64 for r in ledger['accessions'])
    benchmark=json.loads((ROOT/'results/s669_group_holdout.json').read_text())
    train=set(benchmark['split']['train_pdbs'])
    test=set(benchmark['split']['test_pdbs'])
    assert train.isdisjoint(test)
    assert len(train|test)==ledger['s669_mutation_mapped_structures']


def test_paper_uses_measured_benchmark_and_no_false_win():
    paper=(ROOT/'paper/manuscript.tex').read_text()
    assert 'did not outperform' in paper
    assert paper.count('\\begin{equation}')>=10
