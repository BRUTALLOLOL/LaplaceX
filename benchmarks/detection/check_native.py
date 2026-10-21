"""Independent native projection/Jacobian checks and signed/null regression sweep."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'.build/detection-native'))
import _detection_native as native
checks=[]
for rates in ([.5,1.8],[.5,.55],[.03,2.,60.]):
 t=np.geomspace(.01,20,160);y=.8*np.exp(-.4*t)+.3*np.exp(-2*t)+.01*np.sin(t);x=np.log(rates)
 d=native.projected(t,y,x);j=np.array(d['jacobian_columns']).T;h=1e-5
 numerical=np.column_stack([(np.array(native.projected(t,y,x+np.eye(len(x))[k]*h)['residual'])-np.array(native.projected(t,y,x-np.eye(len(x))[k]*h)['residual']))/(2*h) for k in range(len(x))])
 np.testing.assert_allclose(j,numerical,rtol=2e-5,atol=2e-9)
 e=np.exp(-np.outer(t,rates));a=np.linalg.lstsq(e,y,rcond=None)[0]
 np.testing.assert_allclose(d['amplitudes'],a,rtol=1e-10,atol=1e-12)
 checks.append({'rates':rates,'jacobian_max_absolute_error':float(np.max(abs(j-numerical))),'amplitudes_match_numpy_svd':True})
rng=np.random.default_rng(2026102209);rows=[]
for grid in ('uniform','log'):
 t=(np.linspace(.01,20,600) if grid=='uniform' else np.geomspace(.001,1000,600))
 low,high,sep=(.2,6,2) if grid=='uniform' else (.02,60,3)
 for count in (1,2,3,4):
  for sigma in (0,.001,.01):
   for repeat in range(10):
    while True:
     truth=np.sort(np.exp(rng.uniform(np.log(low),np.log(high),count)))
     if count==1 or min(truth[1:]/truth[:-1])>=sep:break
    a=np.exp(rng.uniform(np.log(.2),np.log(2),count))*rng.choice([-1,1],count)
    clean=np.exp(-np.outer(t,truth))@a;y=clean+rng.normal(0,sigma,len(t));result={}
    for route in ('before','after'):
     fit=native.fit(t,y,route);rates=np.array(fit['rates']);amps=np.array(fit['amplitudes']);ix=np.argsort(rates);rates=rates[ix];amps=amps[ix];pred=np.exp(-np.outer(t,rates))@amps
     correct=len(rates)==count;error=float(max(abs(rates/truth-1))) if correct else None;rmse=float(np.sqrt(np.mean((pred-clean)**2)));maximum=float(max(abs(pred-clean))/max(abs(clean)))
     result[route]={'rates':rates.tolist(),'amplitudes':amps.tolist(),'correct_count':correct,'rate_error':error,'clean_rmse':rmse,'joint':bool(correct and (error<=.1 and rmse<=sigma if sigma else error<=1e-5 and maximum<=1e-8))}
    rows.append({'grid':grid,'count':count,'noise_sigma':sigma,'repeat':repeat,'rates':truth.tolist(),'amplitudes':a.tolist(),'methods':result})
null=[]
for grid in ('uniform','log'):
 t=(np.linspace(.01,20,600) if grid=='uniform' else np.geomspace(.001,1000,600))
 for repeat in range(20):
  noise=rng.normal(0,.001,len(t))
  for offset in (0,1):
   methods={}
   for route in ('before','after'):
    fit=native.fit(t,noise+offset,route);rates=fit['rates'];methods[route]={'rates':rates,'passed':len(rates)==offset and (offset==0 or rates[0]==0)}
   null.append({'grid':grid,'repeat':repeat,'offset':offset,'methods':methods})
(HERE/'results').mkdir(exist_ok=True)
(HERE/'results/checks.json').write_text(json.dumps({'source_sha256':hashlib.sha256((ROOT/'cpp/src/critical.cpp').read_bytes()).hexdigest(),'protocol_sha256':hashlib.sha256((HERE/'STRESS_PROTOCOL.md').read_bytes()).hexdigest(),'seed':2026102209,'jacobian_checks':checks,'signed':rows,'null':null},indent=2)+'\n',encoding='utf-8')
for route in ('before','after'):print(route,'signed counts',sum(r['methods'][route]['correct_count'] for r in rows),'joint',sum(r['methods'][route]['joint'] for r in rows),'null',sum(r['methods'][route]['passed'] for r in null),len(null))
