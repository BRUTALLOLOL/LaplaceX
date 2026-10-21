"""Prepare and fit two public measured fluorescence tails with fixed rules."""
from __future__ import annotations
import hashlib
import io
import json
from pathlib import Path
import urllib.request
import zipfile

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
DATA=HERE/'data'
URL='https://zenodo.org/records/20792957/files/Emission%20Decay.zip?download=1'
MD5='65718629f717491b3c043f1f87e22b16'


def prepare():
    DATA.mkdir(exist_ok=True)
    archive=ROOT/'.build/experimental/emission.zip'
    archive.parent.mkdir(parents=True,exist_ok=True)
    if not archive.exists():
        urllib.request.urlretrieve(URL,archive)
    raw=archive.read_bytes()
    assert hashlib.md5(raw).hexdigest()==MD5
    sources=[]
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        for suffix,label,units in [('S1-DF.csv','delayed','microseconds'),('S1-PF.csv','prompt','nanoseconds')]:
            name=next(n for n in z.namelist() if n.endswith(suffix))
            content=z.read(name)
            (DATA/f'{label}_original.csv').write_bytes(content)
            sources.append({'label':label,'archive_member':name,'units':units,
                            'sha256':hashlib.sha256(content).hexdigest()})
    provenance={'doi':'10.5281/zenodo.20792957','url':URL,'license':'CC-BY-4.0',
        'authors':['Stephan Kupfer','Alexander Tombrink','Benedikt Callies','Daniel Kowalczyk','Benedikt Wiedemann'],
        'title':'Supplementary Data: 4CzIPN Photosensitizes the Cobaloxime-catalyzed Light-Driven Hydrogen Evolution Reaction',
        'version':'v3','published':'2026-06-22','archive_md5':MD5,
        'archive_sha256':hashlib.sha256(raw).hexdigest(),'files':sources,
        'protocol_sha256':hashlib.sha256((HERE/'EXPERIMENT_PROTOCOL.md').read_bytes()).hexdigest()}
    (DATA/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf8')
    return provenance


def tail(label,extra_gate=0):
    original=np.genfromtxt(DATA/f'{label}_original.csv',delimiter=',',names=True,encoding='utf8')
    names=original.dtype.names
    t=original[names[0]]; irf=original['IRF']; decay=original['Decay']
    background=float(np.median(decay[-len(decay)//10:]))
    peak=int(np.argmax(decay))
    gate=max(peak+5,int(np.flatnonzero(irf>.01*irf.max())[-1])+1)+extra_gate
    groups=(len(t)-gate)//4
    times=t[gate:gate+groups*4].reshape(-1,4).mean(axis=1)
    values=decay[gate:gate+groups*4].reshape(-1,4).mean(axis=1)-background
    return times,values,{'gate_index':gate,'peak_index':peak,'background':background,
                        'bin_average':4,'samples':groups,'time_min':float(times[0]),'time_max':float(times[-1])}


def main():
    from run_study import evaluate, methods, metadata
    provenance=prepare(); rows=[]
    for label in ['delayed','prompt']:
        for extra in [0,5]:
            t,y,preprocessing=tail(label,extra)
            outcomes={name:evaluate(fn,t,None,y) for name,fn in methods(t,y,sensitivity=True).items()}
            rows.append({'label':label,'extra_gate_bins':extra,'preprocessing':preprocessing,
                'times':t.tolist(),'observations':y.tolist(),'methods':outcomes})
    (HERE/'results').mkdir(exist_ok=True)
    (HERE/'results/experimental.json').write_text(json.dumps({'metadata':metadata(),
        'provenance':provenance,'results':rows},indent=2,allow_nan=False)+'\n',encoding='utf8')
    for row in rows:
        print(row['label'],row['extra_gate_bins'],{k:(v.get('count'),v.get('observed_rmse')) for k,v in row['methods'].items()})


if __name__=='__main__': main()
