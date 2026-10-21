"""Same-process shuffled timing controls; no manuscript generation."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import sys,json,random,time,hashlib,platform
from pathlib import Path
import numpy as np
import laplacex as lx
import scipy
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'.build/detection-native'),str(ROOT/'benchmarks'),str(ROOT/'benchmarks/research')]
import _detection_native as native
from compare import decomposition_cases
from baselines import variable_projection
calls={key:(lambda t,y,key=key:native.fit(t,y,key)) for key in ('before','after','grid')}
calls['varpro']=variable_projection
def public_fit(t,y):
    c=lx.decompose(t,y)
    return {'rates':c.rates,'amplitudes':c.amplitudes}
calls['public_default']=public_fit
rng=random.Random(2026102208);rows=[]
for name,t,clean,y,truth,amps,sigma,seed in decomposition_cases():
 for fn in calls.values():fn(t,y)
 records={key:[] for key in calls}
 for _ in range(9):
  keys=list(calls);rng.shuffle(keys)
  for key in keys:
   start=time.perf_counter();fit=calls[key](t,y);elapsed=time.perf_counter()-start
   rates=np.array(fit['rates'] if isinstance(fit,dict) else fit.rates);amps=np.array(fit['amplitudes'] if isinstance(fit,dict) else fit.amplitudes);ix=np.argsort(rates);rates=rates[ix];amps=amps[ix]
   pred=np.exp(-np.outer(t,rates))@amps;count=len(rates)==len(truth)
   error=float(np.max(abs(rates/truth-1))) if count else None;rmse=float(np.sqrt(np.mean((pred-clean)**2)));maxerr=float(max(abs(pred-clean))/max(abs(clean)))
   passed=count and (error<=.1 and rmse<=sigma if sigma else error<=1e-5 and maxerr<=1e-8)
   records[key].append({'seconds':elapsed,'rates':rates.tolist(),'count':count,'rate_error':error,'clean_rmse':rmse,'passed':bool(passed)})
 row={'id':name,'noise_sigma':sigma,'methods':{key:{'median_seconds':float(np.median([r['seconds'] for r in runs])),'runs':runs} for key,runs in records.items()}};rows.append(row)
 print(name,[(key,round(row['methods'][key]['median_seconds']*1000,3),sum(v['passed'] for v in records[key])) for key in calls],flush=True)
Path('benchmarks/detection/results').mkdir(exist_ok=True)
Path('benchmarks/detection/results/timing.json').write_text(json.dumps({'metadata':{'source_sha256':hashlib.sha256((ROOT/'cpp/src/critical.cpp').read_bytes()).hexdigest(),'platform':platform.platform(),'seed':2026102208,'repeats':9,'blas_threads':1,'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'baselines_sha256':hashlib.sha256((ROOT/'benchmarks/research/baselines.py').read_bytes()).hexdigest()},'results':rows},indent=2)+'\n')
