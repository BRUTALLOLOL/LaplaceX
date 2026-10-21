"""Independent validation and exposed-data replay for the detection revision.

Build controls first, install the candidate wheel, then run --section holdout.
The seed is fixed in PROTOCOL.md. Original research records are never modified.
"""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,hashlib,json,platform,sys
from pathlib import Path
from time import perf_counter
import numpy as np
import scipy
import laplacex as lx
from laplacex import _laplacex
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'.build/detection-native'),str(ROOT/'benchmarks/research')]
import _detection_native as native
from baselines import Fit,variable_projection
from run_study import evaluate,legacy_cases

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def metadata(seed):
 return {'seed':seed,'protocol_sha256':digest(HERE/'PROTOCOL.md'),'script_sha256':digest(__file__),'source_sha256':digest(ROOT/'cpp/src/critical.cpp'),'baselines_sha256':digest(ROOT/'benchmarks/research/baselines.py'),'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),'installed_extension_sha256':digest(_laplacex.__file__),'installed_version_label':lx.__version__,'status':'unreleased candidate; not published','blas_threads':1,'native_build':json.loads((ROOT/'.build/detection-native/build.json').read_text())}

def native_fit(t,y,route):
 d=native.fit(t,y,route);r=np.array(d['rates']);a=np.array(d['amplitudes']);order=np.argsort(r);res=np.exp(-np.outer(t,r))@a-y
 return Fit(r[order],a[order],float(res@res))

def public_fit(t,y):
 d=lx.decompose(t,y);r=d.rates;a=d.amplitudes;order=np.argsort(r);res=d.reconstruct(t)-y
 return Fit(r[order],a[order],float(res@res))

def methods(t,y,baseline=True):
 calls={route:(lambda route=route:native_fit(t,y,route)) for route in ('before','after','grid')}
 calls['public_default']=lambda:public_fit(t,y)
 if baseline:calls['varpro']=lambda:variable_projection(t,y)
 return calls

def save(name,seed,rows):
 (HERE/'results').mkdir(exist_ok=True)
 (HERE/f'results/{name}.json').write_text(json.dumps({'metadata':metadata(seed),'results':rows},indent=2,allow_nan=False)+'\n',encoding='utf-8')

def sweep(seed,name):
 rng=np.random.default_rng(seed);rows=[]
 for grid in ('uniform','log'):
  t=(np.linspace if grid=='uniform' else np.geomspace)(.01,20,600)
  for ratio in (1.1,1.25,1.4,2,3,10):
   rates=np.array([.5,.5*ratio]);amps=np.array([1.,.5]);clean=np.exp(-np.outer(t,rates))@amps
   for sigma in (.0001,.001,.01,.05):
    for repeat in range(20):
     y=clean+rng.normal(0,sigma,len(t));results={key:evaluate(call,t,clean,y,rates) for key,call in methods(t,y).items()}
     rows.append({'grid':grid,'rate_ratio':ratio,'noise_sigma':sigma,'repeat':repeat,'rates':rates.tolist(),'amplitudes':amps.tolist(),'methods':results})
    save(name,seed,rows);print(name,grid,ratio,sigma,len(rows),flush=True)
 save(name,seed,rows)

def legacy():
 rows=[]
 for info,t,clean,y,rates,amps in legacy_cases():
  results={key:evaluate(call,t,clean,y,rates) for key,call in methods(t,y,False).items()}
  rows.append({**info,'rates':rates.tolist(),'amplitudes':amps.tolist(),'methods':results})
  if len(rows)%40==0:save('legacy',None,rows);print('legacy',len(rows),flush=True)
 save('legacy',None,rows)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--section',choices=['holdout','development','legacy'],required=True);args=p.parse_args()
 if args.section=='legacy':legacy()
 else:sweep(2026102217 if args.section=='holdout' else 2026102101,args.section)
