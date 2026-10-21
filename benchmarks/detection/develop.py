"""Development replay of the exposed 960-signal sweep; no baseline records overwritten."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import json,sys,time,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'.build/detection-native'))
import _detection_native as native
old=json.loads((ROOT/'benchmarks/research/results/sweep.json').read_text())['results']
rng=np.random.default_rng(2026102101);rows=[]
for grid in ('uniform','log'):
 t=(np.linspace if grid=='uniform' else np.geomspace)(.01,20,600)
 for ratio in (1.1,1.25,1.4,2,3,10):
  truth=np.array([.5,.5*ratio]);clean=np.exp(-np.outer(t,truth))@np.array([1.,.5])
  for sigma in (.0001,.001,.01,.05):
   for repeat in range(20):
    y=clean+rng.normal(0,sigma,len(t));start=time.perf_counter();fit=native.fit(t,y);elapsed=time.perf_counter()-start
    rates=np.array(fit['rates']);amps=np.array(fit['amplitudes']);ix=np.argsort(rates);rates=rates[ix];amps=amps[ix]
    pred=np.exp(-np.outer(t,rates))@amps
    count=len(rates)==2;error=float(max(abs(rates/truth-1))) if count else None
    rmse=float(np.sqrt(np.mean((pred-clean)**2)))
    previous=old[len(rows)]
    rows.append({'grid':grid,'ratio':ratio,'sigma':sigma,'repeat':repeat,'rates':rates.tolist(),'amplitudes':amps.tolist(),'correct_count':count,'rate_error':error,'clean_rmse':rmse,'joint':count and error<=.1 and rmse<=sigma,'seconds':elapsed,'before_count':previous['methods']['native_log']['correct_count'],'varpro_count':previous['methods']['varpro_bic']['correct_count']})
   print(grid,ratio,sigma,sum(r['correct_count'] for r in rows),len(rows),flush=True)
Path('.build/detection-development.json').write_text(json.dumps(rows,indent=2)+'\n')
records=[]
for current,previous in zip(rows,old):
    after={'status':'ok','seconds':current['seconds'],'rates':current['rates'],
           'amplitudes':current['amplitudes'],'count':len(current['rates']),
           'correct_count':current['correct_count'],
           'max_rate_relative_error':current['rate_error'],
           'clean_rmse':current['clean_rmse']}
    records.append({'grid':current['grid'],'rate_ratio':current['ratio'],
        'noise_sigma':current['sigma'],'repeat':current['repeat'],
        'methods':{'before':previous['methods']['native_log'],
                   'after':after,'varpro':previous['methods']['varpro_bic']}})
output=ROOT/'benchmarks/detection/results/development.json'
output.parent.mkdir(exist_ok=True)
output.write_text(json.dumps({'metadata':{
    'seed':2026102101,'status':'exposed development observations',
    'source_sha256':hashlib.sha256((ROOT/'cpp/src/critical.cpp').read_bytes()).hexdigest(),
    'baselines_sha256':hashlib.sha256((ROOT/'benchmarks/research/baselines.py').read_bytes()).hexdigest(),
    'before_and_varpro_source':'benchmarks/research/results/sweep.json',
    'timing_policy':'Before and baseline timings are archived; do not compare development runtimes.'},
    'results':records},indent=2)+'\n',encoding='utf-8')
print('counts',sum(r['correct_count'] for r in rows),'joint',sum(r['joint'] for r in rows),'seconds',sum(r['seconds'] for r in rows))
print('rate p95/max',np.percentile([r['rate_error'] for r in rows if r['correct_count']],[95,100]))
print('lost',[(r['grid'],r['ratio'],r['sigma'],r['repeat'],r['rates']) for r in rows if r['before_count'] and not r['correct_count']][:20])
