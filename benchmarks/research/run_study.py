"""Run the frozen research protocol; results are checkpointed by section.

python benchmarks/research/run_study.py --section all
The installed release stays unchanged. Build the external native harness first.
"""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[key]='1'
import argparse
import hashlib
import json
import platform
from pathlib import Path
import random
import sys
from time import perf_counter
import warnings

import numpy as np
import scipy
import laplacex as lx
from laplacex import _laplacex as release_native

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUT=HERE/'results'
sys.path.insert(0,str(ROOT/'.build/research-native'))
import _research_native as native
from baselines import Fit, cramer_rao, matrix_pencil, variable_projection


def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def metadata():
    return {'protocol_sha256':digest(HERE/'PROTOCOL.md'),'run_script_sha256':digest(__file__),
        'baselines_sha256':digest(HERE/'baselines.py'),'python':sys.version,'numpy':np.__version__,
        'scipy':scipy.__version__,'platform':platform.platform(),'processor':platform.processor(),
        'laplacex_version':lx.__version__,'installed_extension_sha256':digest(release_native.__file__),
        'native_harness':json.loads((ROOT/'.build/research-native/build.json').read_text()),
        'max_components':5,'retry':False,'blas_threads':1}


def save(section, rows):
    OUT.mkdir(exist_ok=True)
    (OUT/f'{section}.json').write_text(json.dumps({'metadata':metadata(),'results':rows},
        indent=2,allow_nan=False)+'\n',encoding='utf8')


def native_fit(t,y,route='log',points=256):
    d=native.fit(t,y,route,points)
    rates=np.asarray(d['rates']); amps=np.asarray(d['amplitudes']); order=np.argsort(rates)
    residual=np.exp(-np.outer(t,rates))@amps-y
    return Fit(rates[order],amps[order],float(residual@residual))


def public_fit(t,y):
    c=lx.decompose(t,y)
    order=np.argsort(c.rates); residual=c.reconstruct(t)-y
    return Fit(c.rates[order],c.amplitudes[order],float(residual@residual))


def methods(t,y,count=None, sensitivity=False):
    result={'native_log':lambda:native_fit(t,y),'native_linear':lambda:native_fit(t,y,'linear'),
            'native_grid':lambda:native_fit(t,y,'grid'),
            'varpro_bic':lambda:variable_projection(t,y)}
    if sensitivity:
        result['public_default']=lambda:public_fit(t,y)
        result['native_linear_1024']=lambda:native_fit(t,y,'linear',1024)
    if count is not None: result['varpro_known_k']=lambda:variable_projection(t,y,count)
    if np.allclose(np.diff(t),t[1]-t[0],rtol=1e-7,atol=1e-12):
        result['pencil_varpro_bic']=lambda:matrix_pencil(t,y)
    return result


def evaluate(call,t,clean,y,true_rates=None):
    start=perf_counter()
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always'); fit=call()
        elapsed=perf_counter()-start
        pred=fit.values(t)
        if not np.all(np.isfinite(pred)) or not np.all(np.isfinite(fit.rates)):
            raise ValueError('nonfinite fit')
        row={'status':'ok','seconds':elapsed,'rates':fit.rates.tolist(),
             'amplitudes':fit.amplitudes.tolist(),'count':len(fit.rates),
             'observed_rmse':float(np.sqrt(np.mean((pred-y)**2))),
             'warnings':list(dict.fromkeys(str(w.message) for w in caught))}
        if clean is not None:
            row['clean_rmse']=float(np.sqrt(np.mean((pred-clean)**2)))
            row['normalized_max_error']=float(np.max(abs(pred-clean))/max(np.max(abs(clean)),1e-100))
        if true_rates is not None:
            row['correct_count']=len(fit.rates)==len(true_rates)
            row['max_rate_relative_error']=float(np.max(abs(fit.rates/true_rates-1))) if row['correct_count'] else None
            row['log_rate_errors']=np.log(fit.rates/true_rates).tolist() if row['correct_count'] and np.all(fit.rates>0) else None
        return row
    except Exception as error:
        return {'status':'failed','error':str(error),'seconds':perf_counter()-start}


def legacy_cases():
    for filename in ('validation','validation_fresh','validation_holdout'):
        archived=json.loads((ROOT/f'benchmarks/comparison/{filename}.json').read_text())
        rng=np.random.default_rng(archived['seed']); j=0
        for grid in ('log','uniform'):
            for index in range(80):
                count=1+index%4
                lo,hi,separation=(.02,60,3) if grid=='log' else (.2,6,2)
                while True:
                    rates=np.sort(np.exp(rng.uniform(np.log(lo),np.log(hi),count)))
                    if count==1 or min(rates[1:]/rates[:-1])>=separation: break
                amps=np.exp(rng.uniform(np.log(.2),np.log(2),count))
                t=np.geomspace(.001,1000,600) if grid=='log' else np.linspace(.01,20,600)
                clean=np.exp(-np.outer(t,rates))@amps
                sigma=0 if index<40 else .001
                y=clean+rng.normal(0,sigma,len(t))
                previous=archived['results'][j]; j+=1
                np.testing.assert_array_equal(rates,previous['rates'])
                np.testing.assert_array_equal(amps,previous['amplitudes'])
                yield {'id':f"{archived['seed']}-{grid}-{index}",'seed':archived['seed'],
                       'grid':grid,'index':index,'noise_sigma':sigma,
                       'original_correct_count':previous['correct_count'],
                       'original_rate_error':previous['max_rate_relative_error']},t,clean,y,rates,amps


def run_legacy():
    rows=[]
    for info,t,clean,y,rates,amps in legacy_cases():
        result={name:evaluate(fn,t,clean,y,rates) for name,fn in methods(t,y,len(rates),True).items()}
        bound=cramer_rao(t,amps,rates,info['noise_sigma']) if info['noise_sigma'] else None
        rows.append({**info,'rates':rates.tolist(),'amplitudes':amps.tolist(),'cramer_rao':bound,'methods':result})
        if len(rows)%20==0:
            save('legacy',rows);print('legacy',len(rows),flush=True)
    save('legacy',rows)


def run_sweep():
    rng=np.random.default_rng(2026102101); rows=[]
    for grid in ('uniform','log'):
        t=(np.linspace if grid=='uniform' else np.geomspace)(.01,20,600)
        for ratio in (1.1,1.25,1.4,2,3,10):
            rates=np.array([.5,.5*ratio]); amps=np.array([1.,.5])
            clean=np.exp(-np.outer(t,rates))@amps
            for sigma in (.0001,.001,.01,.05):
                bound=cramer_rao(t,amps,rates,sigma)
                for repeat in range(20):
                    y=clean+rng.normal(0,sigma,len(t))
                    result={name:evaluate(fn,t,clean,y,rates) for name,fn in methods(t,y,2).items()}
                    rows.append({'grid':grid,'rate_ratio':ratio,'noise_sigma':sigma,'repeat':repeat,
                                 'cramer_rao':bound,'methods':result})
                save('sweep',rows);print('sweep',grid,ratio,sigma,len(rows),flush=True)
    save('sweep',rows)


def run_timing():
    # Import the archived fixture definitions; do not alter their settings.
    sys.path.insert(0,str(ROOT/'benchmarks'))
    from compare import decomposition_cases
    rows=[]; rng=random.Random(2026102102)
    for name,t,clean,y,rates,amps,sigma,seed in decomposition_cases():
        calls=methods(t,y,sensitivity=True)
        for call in calls.values(): evaluate(call,t,clean,y,rates)
        records={key:[] for key in calls}
        for _ in range(7):
            order=list(calls);rng.shuffle(order)
            for key in order: records[key].append(evaluate(calls[key],t,clean,y,rates))
        rows.append({'id':name,'noise_sigma':sigma,'rates':rates.tolist(),
            'proposals':{route:native.proposals(t,y,route) for route in ('log','linear')},
            'methods':{key:{'runs':runs,'median_seconds':float(np.median([r['seconds'] for r in runs]))}
                       for key,runs in records.items()}})
        save('timing',rows);print('timing',name,flush=True)
    save('timing',rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--section',choices=['all','legacy','sweep','timing'],default='all')
    args=parser.parse_args()
    for section in (['legacy','sweep','timing'] if args.section=='all' else [args.section]):
        globals()['run_'+section]()
