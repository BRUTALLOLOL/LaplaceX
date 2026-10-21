"""Render the detection study and current COMPARISON section from raw outcomes."""
from __future__ import annotations
import hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
SOURCE=hashlib.sha256((ROOT/'cpp/src/critical.cpp').read_bytes()).hexdigest()

def load(name):
 d=json.loads((HERE/f'results/{name}.json').read_text(encoding='utf-8'))
 assert d['metadata']['source_sha256']==SOURCE,(name,'source mismatch')
 return d

def passes(row,sigma):
 if row.get('status')!='ok' or not row.get('correct_count'):return False
 if sigma:return row['clean_rmse']<=sigma and row['max_rate_relative_error']<=.1
 return row['normalized_max_error']<=1e-8 and row['max_rate_relative_error']<=1e-5

def stats(rows,method):
 selected=[r for r in rows if method in r['methods']]
 noisy=[r['methods'][method]['max_rate_relative_error'] for r in selected if r['noise_sigma'] and r['methods'][method].get('correct_count')]
 return {'eligible':len(selected),'counts':sum(bool(r['methods'][method].get('correct_count')) for r in selected),'joint':sum(passes(r['methods'][method],r['noise_sigma']) for r in selected),'failures':sum(r['methods'][method]['status']!='ok' for r in selected),'conditional_noisy_rate_p95':float(np.percentile(noisy,95)) if noisy else None,'conditional_noisy_rate_max':float(max(noisy)) if noisy else None}

def main():
 candidate=json.loads((HERE/'CANDIDATE.json').read_text());assert candidate['source_sha256']==SOURCE
 datasets={name:load(name) for name in ('development','holdout','legacy','timing')}
 assert [len(datasets[k]['results']) for k in ('development','holdout','legacy','timing')]==[960,960,480,12]
 assert datasets['holdout']['metadata']['seed']==2026102217
 balanced=load('balanced_timing');assert len(balanced['results'])==48
 assert balanced['metadata']['timing_protocol_sha256']==hashlib.sha256((HERE/'TIMING_PROTOCOL.md').read_bytes()).hexdigest()
 checks=json.loads((HERE/'results/checks.json').read_text());assert checks['source_sha256']==SOURCE
 assert checks['protocol_sha256']==hashlib.sha256((HERE/'STRESS_PROTOCOL.md').read_bytes()).hexdigest()
 assert datasets['holdout']['metadata']['protocol_sha256']==hashlib.sha256((HERE/'PROTOCOL.md').read_bytes()).hexdigest()
 assert datasets['holdout']['metadata']['script_sha256']==hashlib.sha256((HERE/'validate.py').read_bytes()).hexdigest()
 summary={name:{m:stats(d['results'],m) for m in d['results'][0]['methods']} for name,d in datasets.items() if name!='timing'}
 names={'before':'1.1.0 wavelets','after':'1.2.0 wavelets','grid':'1.2.0 native rate grid','varpro':'Variable projection + BIC','public_default':'1.2.0 public default'}
 lines=['# LaplaceX 1.2.0: component detection','',
 'Measured on 21 October 2026 using the numerical algorithms released in 1.2.0.',
 'Baseline core: commit `3708dd2` (version 1.1.0).',
 f'Numerical source SHA-256: `{SOURCE}`.','',
 '## Count recovery and parameter accuracy','',
 '| Dataset | Method | Correct count | Count + accuracy | Fit exceptions |',
 '|---|---|---:|---:|---:|']
 for dataset,label in [('development','Exposed development'),('holdout','Independent seed'),('legacy','Original mixtures')]:
  for method,st in summary[dataset].items():
   lines.append(f"| {label} | {names[method]} | {st['counts']}/{st['eligible']} | {st['joint']}/{st['eligible']} | {st['failures']} |")
 lines+=['',
 'The development sweep reuses seed 2026102101; its before/baseline measurements',
 'come from the unchanged original research records. It is not independent evidence.',
 'The independent sweep uses seed **2026102217**, fixed before its outcomes,',
 'and reruns every comparator on exactly the same observations: two components,',
 'amplitudes 1 and 0.5, rate ratios 1.1–10, Gaussian standard deviations',
 '0.0001–0.05, uniform/geometric sampling, 600 observations, 20 realizations per cell.',
 'All 960 outputs are retained. A rejected first candidate and its partly examined',
 'seed-2026102207 observations were treated as development, never as final validation.','',
 'A noisy joint success requires the true count, maximum rate error ≤10%, and',
 'clean-signal RMSE ≤noise standard deviation. Clean cases require the true count,',
 'maximum rate error ≤1e-5, and normalized maximum reconstruction error ≤1e-8.',
 'The native wavelet routes disable the Python rate-grid retry; public default',
 'results are reported separately. The original 480-mixture variable-projection',
 'count result remains 479/480 in the frozen research report; it was not rerun here.','',
 '### Count recovery does not guarantee accurate rates','',
 '| Independent sweep, correct-count fits only | Rate error, 95th percentile | Worst rate error |',
 '|---|---:|---:|']
 for method in ('before','after','varpro'):
  st=summary['holdout'][method];lines.append(f"| {names[method]} | {100*st['conditional_noisy_rate_p95']:.2f}% | {100*st['conditional_noisy_rate_max']:.2f}% |")
 hold=datasets['holdout']['results'];paired={}
 for endpoint in ('count','joint'):
  def ok(r,m):return bool(r['methods'][m].get('correct_count')) if endpoint=='count' else passes(r['methods'][m],r['noise_sigma'])
  paired[endpoint]={'after_only':sum(ok(r,'after') and not ok(r,'varpro') for r in hold),'varpro_only':sum(ok(r,'varpro') and not ok(r,'after') for r in hold)}
 lines+=['',f"Paired independent count outcomes: revised-only successes {paired['count']['after_only']}, variable-projection-only successes {paired['count']['varpro_only']}.",
 f"Paired joint outcomes: revised-only successes {paired['joint']['after_only']}, variable-projection-only successes {paired['joint']['varpro_only']}.",
 'These are descriptive results on fixed synthetic conditions, not proof of',
 'statistically significant or general accuracy superiority. Close noisy modes',
 'can remain unidentifiable even when a two-mode fit is selected.','',
 '## Speed with accuracy checked on every timed fit','',
 'Nine shuffled round-robin repetitions after one warmup per method and case;',
 'one BLAS thread. Times below are medians in milliseconds. The native routes',
 'share the wrapper and compiler; the revised wavelet/grid routes also share',
 'the projected optimizer, admissibility retry, and model-selection rules.',
 'Variable projection uses the original SciPy/Python implementation with two',
 'deterministic starts per count, K=1–5 and BIC. Its comparison includes language',
 'and implementation costs. Validation-run elapsed times are not speed evidence.','',
 '| Case | 1.1.0 native | Revised wavelets | Revised grid | Variable projection | Grid / revised | Varpro / revised |',
 '|---|---:|---:|---:|---:|---:|---:|']
 timing=datasets['timing']['results'];ratios=[];regressions=[]
 for row in timing:
  d=row['methods'];assert all(len(v['runs'])==9 and all(r['passed'] for r in v['runs']) for v in d.values()),row['id']
  b,a,g,v=(d[m]['median_seconds']*1000 for m in ('before','after','grid','varpro'))
  ratios.append({'id':row['id'],'before_ms':b,'after_ms':a,'grid_over_after':g/a,'varpro_over_after':v/a,'before_over_after':b/a})
  if a>b:regressions.append((row['id'],a/b))
  lines.append(f"| {row['id']} | {b:.3f} | {a:.3f} | {g:.3f} | {v:.3f} | {g/a:.2f}× | {v/a:.2f}× |")
 lines+=['',
 'Every one of the nine fits per fixture and method passed the count and accuracy',
 'targets. The public Python default is also timed and retained in the raw file.',
 'The exact single-exponential shortcut is shared by both native routes; its',
 'speed is not evidence for a wavelet advantage. The remaining native grid ratios',
 'are the cleaner evidence for proposal-search savings.','',
 'Per-case changes relative to 1.1.0 include slowdowns: '+', '.join(f'{name} ({factor:.2f}× the old time)' for name,factor in regressions)+'.',
 'On these twelve fixtures the revised wavelets remain faster than both comparators,',
 'but do not preserve latency on every existing workload.','',
 'Hardware: Intel Core i7-7700HQ, Windows, Balanced power plan; clock frequency was',
 'not fixed. See the original hardware record and raw metadata for Python, NumPy,',
 'SciPy, compiler command, source hashes, and the installed extension checksum.','',
 '## Changes and checks','',
 '- Fit linear amplitudes by scaled, reorthogonalized QR at every nonlinear rate',
 '  step. The projected rate Jacobian includes the residual-dependent correction;',
 '  normal equations are not formed.',
 '- Search both sides of an unresolved one-mode fit even when residual wavelet',
 '  ridges fall below the noise threshold. Screen cheaply, then retain the existing',
 '  full-data BIC/AICc decision.',
 '- Cap screening variance by both measured noise and fitted residual variance.',
 '  An incomplete model must not hide a weak signed component by calling model',
 '  mismatch noise.',
 '- Retry the original joint optimizer only for projected results violating the',
 '  existing amplitude/coalescence limits. No new component count is prescribed.','']
 for method in ('before','after'):
  signed=checks['signed'];null=checks['null'];c=sum(r['methods'][method]['correct_count'] for r in signed);j=sum(r['methods'][method]['joint'] for r in signed);n=sum(r['methods'][method]['passed'] for r in null)
  lines.append(f"- {names[method]}: signed mixtures {c}/240 correct counts, {j}/240 joint successes; noise/constant checks {n}/80.")
 lines+=['',
 'The signed supplement is development/regression coverage, not another untouched',
 'validation claim. The two null/constant failures occur in both versions; this',
 'revision does not eliminate all false modes. Native QR amplitudes were checked',
 'against NumPy SVD least squares and projected Jacobians against central finite',
 'differences, including close rates. **163 tests pass** against the rebuilt wheel.','',
 '![Independent counts and native timing controls](detection/results/comparison.png)','',
 '## Reproduction and raw records','',
 '```powershell',
 'python -m pip install numpy scipy matplotlib pybind11',
 '# Windows build helper also requires: python -m pip install ziglang',
 '# The frozen case generator imports the original research harness.',
 'python benchmarks/research/build_native.py --compiler .venv/Lib/site-packages/ziglang/zig.exe',
 '# On Linux/macOS use --compiler c++ instead.',
 'python benchmarks/detection/build_native.py',
 '# Install LaplaceX 1.2.0 before validating the public API.',
 'python benchmarks/detection/develop.py',
 'python benchmarks/detection/check_native.py',
 'python benchmarks/detection/validate.py --section legacy',
 'python benchmarks/detection/validate.py --section holdout',
 '# Run timing separately, after validation and other CPU-intensive jobs finish.',
 'python benchmarks/detection/time_controls.py',
 'python benchmarks/detection/time_balanced.py',
 'python benchmarks/detection/summarize.py',
 '```','',
 'Protocols: [main](detection/PROTOCOL.md), [signed/null checks](detection/STRESS_PROTOCOL.md).',
 'Supplemental [balanced timing protocol](detection/TIMING_PROTOCOL.md) and [raw records](detection/results/balanced_timing.json).',
 'Raw data: [development](detection/results/development.json),',
 '[independent sweep](detection/results/holdout.json), [legacy replay](detection/results/legacy.json),',
 '[timing](detection/results/timing.json), [native/signed/null checks](detection/results/checks.json).',
 'Use a Git checkout containing commit 3708dd2 for the before control. Source',
 'checksums refer to checkout bytes, whose line endings are specified in .gitattributes.',
 'The original measurements used a pre-release build carrying the 1.1.0 label.',
 'Frozen raw metadata retains that provenance. The numerical source checksum',
 'matches release 1.2.0; only packaging/version metadata changed afterward.','']
 balanced_lines=['### Balanced timing on the independent sweep', '',
 'Take repeat 0 from every grid/ratio/noise cell, fixed before timing: 48 cases.',
 'Use three shuffled rounds after warmup and retain wrong counts and all errors.',
 'This is a subset of the accuracy sweep, not another independent accuracy sample.',
 'The totals below sum the per-case median elapsed times; they are not a single',
 'wall-clock run. No case was selected for passing or for being fast.', '',
 '| Method | Sum of case medians, ms | Median case, ms | Correct counts | Count + accuracy |',
 '|---|---:|---:|---:|---:|']
 balanced_summary={}
 for method in ('before','after','grid','varpro'):
  cases=balanced['results'];medians=np.array([r['methods'][method]['median_seconds'] for r in cases])
  st={'sum_medians_seconds':float(sum(medians)),'median_case_seconds':float(np.median(medians)),
      'counts':sum(bool(r['methods'][method]['runs'][0].get('correct_count')) for r in cases),
      'joint':sum(passes(r['methods'][method]['runs'][0],r['noise_sigma']) for r in cases)}
  balanced_summary[method]=st
  balanced_lines.append(f"| {names[method]} | {1000*st['sum_medians_seconds']:.2f} | {1000*st['median_case_seconds']:.3f} | {st['counts']}/48 | {st['joint']}/48 |")
 ratio=balanced_summary['varpro']['sum_medians_seconds']/balanced_summary['after']['sum_medians_seconds']
 grid_ratio=balanced_summary['after']['sum_medians_seconds']/balanced_summary['grid']['sum_medians_seconds']
 slower=sum(r['methods']['after']['median_seconds']>r['methods']['before']['median_seconds'] for r in balanced['results'])
 balanced_lines+=['',f'The revised wavelet route uses **{ratio:.2f}× less cumulative time** than the Python variable-projection baseline.',
 f'However, the matched native grid uses **{grid_ratio:.2f}× less cumulative time than wavelets** on this harder subset.',
 f'The revised wavelet route is slower than 1.1.0 in {slower}/48 cells.',
 'The results support retaining a speed advantage over this Python comparator,',
 'not a universal wavelet-search speed advantage. Matching native-grid detection',
 'accuracy does not imply beating its runtime. All raw timed outputs are retained.','']
 where=next(i for i,v in enumerate(lines) if v.startswith('Hardware:'))
 lines[where:where]=balanced_lines
 section='\n'.join(lines)
 (HERE/'COMPARISON_SECTION.md').write_text(section,encoding='utf-8')
 (HERE/'REPORT.md').write_text(section.replace('(detection/','('),encoding='utf-8')
 summary.update({'source_sha256':SOURCE,'paired_holdout':paired,'timing_ratios':ratios,'balanced_timing':balanced_summary})
 (HERE/'results/summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
 # Publication-independent benchmark figure; no article paths or files touched.
 fig,axes=plt.subplots(1,2,figsize=(11,4.6),gridspec_kw={'width_ratios':[1,1.6]})
 for i,(key,label) in enumerate([('before','1.1.0'),('after','1.2.0'),('varpro','Variable projection')]):
  st=summary['holdout'][key];axes[0].bar(i,st['counts'],color=['#8c979e','#267589','#ad6533'][i]);axes[0].text(i,st['counts']+7,str(st['counts']),ha='center',fontsize=9)
 axes[0].set(ylim=(0,1000),ylabel='Correct component counts / 960',title='Independent sweep');axes[0].set_xticks(range(3),['1.1.0','Revised','Variable\nprojection'])
 y=np.arange(12)
 for key,label,color,delta in [('before','1.1.0 native','#8c979e',-.15),('after','Revised wavelets','#267589',0),('grid','Revised native grid','#ad6533',.15)]:
  med=np.array([r['methods'][key]['median_seconds']*1000 for r in timing]);q=np.array([np.percentile([v['seconds']*1000 for v in r['methods'][key]['runs']],[25,75]) for r in timing]);axes[1].errorbar(med,y+delta,xerr=np.array([med-q[:,0],q[:,1]-med]),fmt='o',ms=3,label=label,color=color,capsize=2)
 axes[1].set_yticks(y,[r['id'] for r in timing],fontsize=8);axes[1].invert_yaxis();axes[1].set_xscale('log');axes[1].set(xlabel='Milliseconds, median and interquartile range',title='Matched native timing controls');axes[1].legend(fontsize=8,loc='best')
 for ax in axes:ax.spines[['top','right']].set_visible(False)
 fig.tight_layout();fig.savefig(HERE/'results/comparison.png',dpi=180);plt.close(fig)
 print('Saved detection report, current comparison section, summary and benchmark figure.')
 print(json.dumps(summary['holdout'],indent=2))

if __name__=='__main__':main()
