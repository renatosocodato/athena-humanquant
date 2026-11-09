import subprocess
import json
import sys
from pathlib import Path
import pytest

def test_example_and_no_overwrite(tmp_path):
    out=tmp_path/'example'
    command=[sys.executable,'-m','athena_humanquant','example','--output',str(out)]
    subprocess.run(command,check=True,capture_output=True)
    assert (out/'result.json').is_file()
    result=json.loads((out/'result.json').read_text())
    assert (result['bootstrap_seed'],result['treatment'],result['control']) == (2026,'treated','control')
    second=subprocess.run(command,capture_output=True)
    assert second.returncode!=0


def test_invalid_input_does_not_leave_output_directory(tmp_path):
    source=tmp_path/'invalid.csv';source.write_text('invalid\n1\n')
    out=tmp_path/'out'
    result=subprocess.run([sys.executable,'-m','athena_humanquant','run','--input',str(source),'--output',str(out)],capture_output=True,text=True)
    assert result.returncode==2
    assert not out.exists()
    assert 'Traceback' not in result.stderr


def test_package_specific_help_has_no_unrelated_metadata_argument():
    result=subprocess.run([sys.executable,'-m','athena_humanquant','run','--help'],capture_output=True,text=True,check=True)
    assert '--metadata' not in result.stdout


def test_csv_preserves_distinct_leading_zero_donor_identifiers(tmp_path):
    import pandas as pd
    source=tmp_path/'input.csv'
    pd.DataFrame({'donor':['001','01','002','02'],'field':['1']*4,
        'condition':['control','control','treated','treated'],'value':[1,2,4,5]}).to_csv(source,index=False)
    output=tmp_path/'output'
    subprocess.run([sys.executable,'-m','athena_humanquant','run','--input',str(source),
        '--output',str(output),'--draws','100'],check=True,capture_output=True)
    donors=pd.read_csv(output/'donors.csv',dtype={'donor':'string'})
    assert set(donors.donor)=={'001','01','002','02'}
