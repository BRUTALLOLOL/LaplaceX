"""Time one fixed observation per independent sweep cell; preserve all outcomes."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import hashlib,json,random
from pathlib import Path
import numpy as np
from validate import HERE,ROOT,metadata,native_fit,variable_projection,evaluate
rng=np.random.default_rng(2026102217);order=random.Random(2026102220);rows=[]
validation=json.loads((HERE/'results/holdout.json').read_text())['results'];index=0
for grid in ('uniform','log'):
 t=(np.linspace if grid=='uniform' else np.geomspace)(.01,20,600)
 for ratio in (1.1,1.25,1.4,2,3,10):
  truth=np.array([.5,.5*ratio]);clean=np.exp(-np.outer(t,truth))@np.array([1.,.5])
  for sigma in (.0001,.001,.01,.05):
   for repeat in range(20):
    y=clean+rng.normal(0,sigma,len(t));previous=validation[index];index+=1
    if repeat:continue
    calls={route:(lambda route=route:native_fit(t,y,route)) for route in ('before','after','grid')};calls['varpro']=lambda:variable_projection(t,y)
    for call in calls.values():call()
    result={key:[] for key in calls}
    for _ in range(3):
     keys=list(calls);order.shuffle(keys)
     for key in keys:
      fit=evaluate(calls[key],t,clean,y,truth)
      assert fit.get('count')==previous['methods'][key].get('count'),(grid,ratio,sigma,key)
      result[key].append(fit)
    rows.append({'grid':grid,'rate_ratio':ratio,'noise_sigma':sigma,'repeat':0,'methods':{key:{'median_seconds':float(np.median([r['seconds'] for r in runs])),'runs':runs} for key,runs in result.items()}})
    print(grid,ratio,sigma,len(rows),flush=True)
    (HERE/'results/balanced_timing.json').write_text(json.dumps({'metadata':{**metadata(2026102217),'timing_protocol_sha256':hashlib.sha256((HERE/'TIMING_PROTOCOL.md').read_bytes()).hexdigest(),'selection':'repeat 0 from every sweep cell; no correctness filtering','rounds':3,'shuffle_seed':2026102220},'results':rows},indent=2)+'\n',encoding='utf-8')
