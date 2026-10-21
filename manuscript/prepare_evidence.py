"""Regenerate article assets from the completed frozen research records.

Run: .venv/Scripts/python manuscript/prepare_evidence.py
No fits, timings, library builds, or published numerical code are changed.
"""
from pathlib import Path
import hashlib
import json
import runpy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
runpy.run_path(str(ROOT/'benchmarks/research/summarize.py'),run_name='__main__')

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,
                     'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,2,figsize=(8,2.8),gridspec_kw={'wspace':.3})
t=np.linspace(0,100,2500); x=np.linspace(-7,6,2500)
for rate,color in zip([.03,2,60],['#277687','#ae6633','#596093']):
    axes[0].plot(t,np.exp(-rate*t),color=color,label=f'rate {rate:g}')
    z=x+np.log(rate); axes[1].plot(x,np.exp(z-np.exp(z)),color=color)
axes[0].set(xlabel='Time',ylabel='Component amplitude',title='Ordinary time')
axes[1].set(xlabel='Logarithmic time',ylabel='Negative log-time derivative',title='Translated derivative profiles')
axes[0].legend(frameon=False,fontsize=8)
for ax in axes: ax.grid(alpha=.15)
for suffix in ['pdf','png']: fig.savefig(HERE/f'figures/log_time_profiles.{suffix}',dpi=220,bbox_inches='tight')
plt.close(fig)

evidence=json.loads((ROOT/'benchmarks/research/results/summary.json').read_text(encoding='utf8'))
evidence['release_commit']='ff164d309a89b72eaafb2c58eb83430399ba9b24'
evidence['manuscript_source_sha256']=hashlib.sha256((HERE/'main.tex').read_bytes()).hexdigest()
evidence['policy']='Assets derived from frozen records; no new fits or timings. Research code and manuscript accompany version 1.1.0.'
(HERE/'evidence.json').write_text(json.dumps(evidence,indent=2,allow_nan=False)+'\n',encoding='utf8')
print('Saved four manuscript figures, three tables, and evidence.json.')
