import pandas as pd
import pytest
from athena_humanquant.core import donor_difference

def fixture():
    return pd.DataFrame([{'donor':str(i),'field':str(j),'condition':'control' if i<2 else 'treated','value':i+j} for i in range(4) for j in range(2)])

def test_known_donor_difference_and_bootstrap_reproducible():
    r,d=donor_difference(fixture(),'treated','control',draws=100)
    assert r['difference']==2 and r['n_treatment_donors']==2 and len(d)==4
    assert r==donor_difference(fixture(),'treated','control',draws=100)[0]

def test_replicating_cells_does_not_change_donor_weight_or_precision():
    x=fixture();r,_=donor_difference(x,'treated','control',draws=100)
    duplicated=pd.concat([x,x.loc[x.donor.eq('0')]],ignore_index=True)
    assert donor_difference(duplicated,'treated','control',draws=100)[0]==r

def test_conflicting_condition_rejected():
    x=fixture();x.loc[0,'condition']='treated'
    with pytest.raises(ValueError):donor_difference(x,'treated','control',draws=100)

def test_one_donor_is_not_cell_level_replication():
    x=fixture();x.loc[x.condition.eq('treated'),'donor']='same'
    with pytest.raises(ValueError):donor_difference(x,'treated','control',draws=100)


def test_field_weights_are_equal_even_with_different_cell_counts():
    x=fixture()
    # Replicate all cells from one field, leaving its mean unchanged.
    replicated=pd.concat([x,x.loc[(x.donor=='0') & (x.field=='0')]],ignore_index=True)
    result,_=donor_difference(x,'treated','control',draws=100)
    assert donor_difference(replicated,'treated','control',draws=100)[0]==result


def test_additional_condition_is_not_silently_ignored():
    x=fixture();extra=x.iloc[[0]].assign(donor='new',condition='other')
    with pytest.raises(ValueError,match='exactly'):
        donor_difference(pd.concat([x,extra]),'treated','control',draws=100)


@pytest.mark.parametrize('draws',[100.5,True,0,99])
def test_draw_count_is_validated(draws):
    with pytest.raises(ValueError,match='draws'):
        donor_difference(fixture(),'treated','control',draws=draws)


def test_schema_and_blank_ids_rejected():
    with pytest.raises(ValueError,match='Missing required columns'):
        donor_difference(fixture().drop(columns='field'),'treated','control')
    x=fixture();x.loc[0,'donor']=' '
    with pytest.raises(ValueError,match='nonempty'):
        donor_difference(x,'treated','control')


def test_leave_one_donor_out_reveals_fragility_at_three_per_group():
    x=pd.DataFrame({'donor':list('abcdef'),'field':['1']*6,
        'condition':['control']*3+['treated']*3,'value':[1,2,30,10,11,12]})
    result,_=donor_difference(x,'treated','control',draws=100)
    assert result['difference']==0
    assert result['leave_one_donor_out_max']==9.5
    assert result['leave_one_donor_out_min']==-5
