"""Execute preparation only; fail if any execution flag becomes enabled."""
import ast
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
ROOT=Path(__file__).resolve().parent
path=ROOT/'triga_single_FE101_B1_1965_5x5_CLEAN.ipynb'
nb=json.loads(path.read_text())
def validate(output_root):
    ns={'__name__':'__fe101_5x5_preparation__','FE101_CASE_DIR':ROOT,'FE101_OUTPUT_ROOT':Path(output_root)}
    flags=['RUN_DEPLETION','RUN_POSTPROCESS','RUN_BENCHMARKS']
    for i,c in enumerate(nb['cells']):
        if c['cell_type']!='code': continue
        assert not any(ns.get(k,False) for k in flags)
        code=''.join(c['source'])
        ast.parse(code)
        exec(compile(code,f'{path.name}:cell_{i}','exec'),ns)
        assert not any(ns.get(k,False) for k in flags)
    assert 'operator' not in ns and 'integrator' not in ns
    assert not ns['RUN_DIR'].exists()
    for p,h in ns['protected_hashes'].items():
        assert hashlib.sha256((ROOT.parent/p).read_bytes()).hexdigest()==h
    print('PASS: all execution switches false; no operator/integrator/run directory; protected inputs unchanged')
    return ns

if __name__ == "__main__":
    with TemporaryDirectory(prefix="fe101_5x5_depletion_prepare_") as output_root:
        validate(output_root)
