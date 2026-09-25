import json
import pandas as pd
import pytest
from protein_rescue.audit import audit


def test_group_audit(tmp_path):
    pred=[dict(name=f'{p}-{i}',pdb=p,experimental_ddg=float(i),gnn=1.+i) for p,n in [('A',4),('B',2)] for i in range(n)]
    results=tmp_path/'pred.json'; results.write_text(json.dumps(dict(predictions=pred)))
    csv=tmp_path/'source.csv'; pd.DataFrame([dict(name=r['name'],**{'GeoDDG-Seq_dir':r['experimental_ddg'],'GeoDDG-3D_dir':r['experimental_ddg']}) for r in pred]).to_csv(csv,index=False)
    output=tmp_path/'audit.json'
    result=audit(results,csv,output,samples=100)
    assert result['rows']==6 and result['protein_groups']==2
    assert result['row_weighted_delta_mae']==1
    assert result['protein_equal_delta_mae']==1
    assert output.is_file()
