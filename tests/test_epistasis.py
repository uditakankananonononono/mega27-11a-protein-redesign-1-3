import math
import pytest
from protein_rescue.epistasis import double_mutant_cycle


def test_published_p53_suppressor_interactions():
    x=double_mutant_cycle(1.22,-1.37,0.0,.22,.23,.22)
    assert x['interaction_ddg']==pytest.approx(.15)
    assert x['interaction_se_independence_assumption']==pytest.approx(math.sqrt(.22**2+.23**2+.22**2))
    y=double_mutant_cycle(1.69,2.75,1.72)
    assert y['interaction_ddg']==pytest.approx(-2.72)


def test_partial_se_rejected():
    with pytest.raises(ValueError):double_mutant_cycle(1,2,3,se_a=.1)
    with pytest.raises(ValueError):double_mutant_cycle(float('nan'),2,3)
