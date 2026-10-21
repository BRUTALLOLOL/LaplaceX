"""Check external controls against analytic peaks and the installed release."""
import json
from pathlib import Path
import sys
import numpy as np
import laplacex as lx

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'.build/research-native'))
import _research_native as native

t=np.linspace(.01,20,600)
checks=[]
for rate in [.3,.5,1.]:
    y=np.exp(-rate*t)
    p=native.proposals(t,y,'linear')
    strongest=max(zip(p['rates'],p['votes']),key=lambda pair:pair[1])[0]
    assert abs(strongest/rate-1)<.005
    checks.append({'analytic_single_rate':rate,'linear_peak_rate':strongest})
for rates,amplitudes in [([.3,3],[2,.5]),([.5,.7],[2,.5]),([.03,2,60],[.7,1.3,.5])]:
    times=t if len(rates)==2 else np.geomspace(.002,400,600)
    y=np.exp(-np.outer(times,rates))@amplitudes
    reference=lx.decompose(times,y,auto_deflate=False)
    actual=native.fit(times,y,'log')
    np.testing.assert_allclose(np.sort(actual['rates']),np.sort(reference.rates),rtol=1e-12,atol=1e-14)
    checks.append({'release_equivalence_rates':rates,'passed':True})
for route in ['log','linear']:
    assert native.proposals(t,np.ones(len(t)),route)['rates']==[]
checks.append({'constant_cancellation_both_coordinates':True})
out=ROOT/'benchmarks/research/results/native_checks.json'
out.write_text(json.dumps(checks,indent=2)+'\n',encoding='utf8')
print('Passed analytic peak, constant-cancellation, and frozen-release equivalence checks.')
