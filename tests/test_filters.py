from protein_rescue.filters import safe_to_rank


def test_abstain_when_target_not_validated():
    v={'pdb_accession':'1TSR','sequence_aligned':True,'numbering_match':True,'mutant_assay_ddg_sign_known':True,'relevant_training_proteins':52,'relevant_training_rows':420,'target_double_mutant_assay_supported':False}
    x=safe_to_rank(v)
    assert not x['eligible'] and 'No target-specific double-mutant assay validation' in x['reasons']
    v['target_double_mutant_assay_supported']=True
    assert safe_to_rank(v)['eligible']


def test_partial_context_refuses_rank():
    result=safe_to_rank({'pdb_accession':'2C9V'})
    assert not result['eligible'] and len(result['reasons'])>=5
