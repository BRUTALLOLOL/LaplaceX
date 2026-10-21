"""Derive article assets and the research report from completed frozen records."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
ART=ROOT/'manuscript'
FIG=ART/'figures'; TABLE=ART/'tables'
LABELS={'native_log':'Log wavelets','native_linear':'Linear-coordinate wavelets',
        'native_grid':'Native rate grid','varpro_bic':'Variable projection + BIC',
        'pencil_varpro_bic':'Pencil + variable projection',
        'public_default':'Public default (retry enabled)',
        'native_linear_1024':'Linear coordinates, 1024 points'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,
                     'axes.spines.top':False,'axes.spines.right':False})


def load(name): return json.loads((HERE/f'results/{name}.json').read_text())


def passes(row,sigma):
    if row['status']!='ok' or not row.get('correct_count'): return False
    if sigma: return row['clean_rmse']<=sigma and row['max_rate_relative_error']<=.1
    return row['normalized_max_error']<=1e-8 and row['max_rate_relative_error']<=1e-5


def statistics(rows,method):
    selected=[r for r in rows if method in r['methods']]
    ok=[r['methods'][method] for r in selected if r['methods'][method].get('correct_count')]
    noisy=[r['methods'][method]['max_rate_relative_error'] for r in selected
           if r['noise_sigma'] and r['methods'][method].get('correct_count')]
    return {'eligible':len(selected),'correct_counts':len(ok),
            'passes_count_and_accuracy':sum(passes(r['methods'][method],r['noise_sigma']) for r in selected),
            'failures':sum(r['methods'][method]['status']!='ok' for r in selected),
            'noisy_correct_counts':len(noisy),'noisy_rate_error_median':float(np.median(noisy)),
            'noisy_rate_error_p95':float(np.quantile(noisy,.95)),
            'noisy_rate_error_max':float(max(noisy))}


def table(name,caption,label,columns,header,rows):
    (TABLE/f'{name}.tex').write_text('\n'.join([
        r'\begin{table}[htbp]',r'\centering\small',r'\caption{'+caption+r'}\label{'+label+'}',
        r'\begin{tabular}{@{}'+columns+r'@{}}',r'\toprule',header+r' \\',r'\midrule',
        *[row+r' \\' for row in rows],r'\bottomrule',r'\end{tabular}',r'\end{table}','']),encoding='utf8')


def figure(fig,name):
    for suffix in ['pdf','png']: fig.savefig(FIG/f'{name}.{suffix}',bbox_inches='tight',dpi=210)
    plt.close(fig)


def main():
    legacy=load('legacy'); sweep=load('sweep'); timing=load('timing'); experiment=load('experimental')
    assert len(legacy['results'])==480 and len(sweep['results'])==960 and len(timing['results'])==12
    protocol=hashlib.sha256((HERE/'PROTOCOL.md').read_bytes()).hexdigest()
    for dataset in [legacy,sweep,timing]:
        assert dataset['metadata']['protocol_sha256']==protocol
        assert dataset['metadata']['run_script_sha256']==hashlib.sha256((HERE/'run_study.py').read_bytes()).hexdigest()
        assert dataset['metadata']['baselines_sha256']==hashlib.sha256((HERE/'baselines.py').read_bytes()).hexdigest()
    keys=['native_log','native_linear','native_grid','varpro_bic','pencil_varpro_bic']
    summary={'legacy':{key:statistics(legacy['results'],key) for key in keys+['public_default','native_linear_1024']},
             'sweep':{key:statistics(sweep['results'],key) for key in keys},
             'metadata':timing['metadata'],'timing':[],'cramer_rao_cells':[]}
    rows=[]
    short={'native_log':'Log wavelets','native_linear':'Linear wavelets','native_grid':'Native grid',
           'varpro_bic':'Variable projection','pencil_varpro_bic':'Pencil + refinement'}
    for key in keys:
        old=summary['legacy'][key]; new=summary['sweep'][key]
        rows.append(f"{short[key]} & {old['correct_counts']}/{old['eligible']} & "
                    f"{100*old['noisy_rate_error_p95']:.2f} & {new['correct_counts']}/{new['eligible']} & "
                    f"{new['passes_count_and_accuracy']}/{new['eligible']}")
    table('research_counts','Count recovery on the reused 480-mixture set and the new 960-observation sweep. '
          'Legacy rate errors are the conditional 95th percentile of the largest per-mixture relative error, '
          'in percent, for noisy correctly counted cases. The final column also requires the stated '
          'reconstruction and rate-error targets. Pencil denominators include only uniform grids.',
          'tab:research-counts','lrrrr','Method & Legacy counts & Legacy 95th (\\%) & Sweep counts & Sweep pass',rows)

    labels=['Single','Two uniform','Wide three A','Wide two','Wide three B','Close pair',
            'Weak fast','Weak middle','Weak slow','Noise 7','Noise 19','Noise 41']
    timing_rows=[]
    for case,label in zip(timing['results'],labels,strict=True):
        ms={key:r['median_seconds']*1000 for key,r in case['methods'].items()}
        validity={key:all(passes(r,case['noise_sigma']) for r in record['runs'])
                  for key,record in case['methods'].items()}
        ratios={key:ms[key]/ms['native_log'] for key in ms}
        summary['timing'].append({'id':case['id'],'label':label,'ms':ms,'passing':validity,'ratios':ratios})
        def cell(key): return f"{ms[key]:.2f}"+('' if validity[key] else r'$^{\dagger}$')
        timing_rows.append(f"{label} & {cell('native_log')} & {cell('native_linear')} & "
                          f"{cell('native_grid')} & {cell('varpro_bic')}")
    table('research_timing','Median milliseconds from seven shuffled repetitions after warmup. '
          'Log, linear, and grid routes use the same compiled optimizer and wrapper; variable projection '
          'uses SciPy callbacks. All fixtures are development-exposed. $\\dagger$ marks a route that '
          'fails the count/accuracy target in at least one repetition; its time is still reported.',
          'tab:research-timing','lrrrr','Fixture & Log & Linear & Grid & VarPro',timing_rows)
    fig,axes=plt.subplots(1,2,figsize=(8.0,3.9),sharey=True,gridspec_kw={'wspace':.13})
    for ax,key,title,color in zip(axes,['native_grid','native_linear'],
                                  ['Same-core rate grid / log','Same-core linear / log'],['#277687','#ae6633'],strict=True):
        values=[r['ratios'][key] for r in summary['timing']]
        ax.barh(np.arange(12),values,color=color,height=.65)
        ax.axvline(1,color='black',linestyle='--',linewidth=.8)
        ax.set_xscale('log');ax.set_xlim(.5,max(values)*2)
        ax.set_title(title,fontsize=10);ax.set_xlabel('Time ratio; >1 favors log proposals')
        ax.set_yticks(np.arange(12),labels);ax.grid(axis='x',alpha=.15);ax.set_axisbelow(True)
        for i,value in enumerate(values):
            label=f'{value:.2f}'+('' if summary['timing'][i]['passing'][key] else ' *')
            ax.text(value*(1.05 if value>=1 else .94),i,label,va='center',
                    ha='left' if value>=1 else 'right',fontsize=8)
    axes[0].invert_yaxis();figure(fig,'matched_controls')

    ratios=[1.1,1.25,1.4,2,3,10];sigmas=[.0001,.001,.01,.05]
    fig,axes=plt.subplots(3,2,figsize=(7.8,7.8),gridspec_kw={'hspace':.55,'wspace':.35})
    for row,grid in enumerate(['uniform','log']):
        for col,key in enumerate(['native_log','varpro_bic']):
            a=np.zeros((4,6))
            for i,sigma in enumerate(sigmas):
                for j,ratio in enumerate(ratios):
                    cell=[r for r in sweep['results'] if r['grid']==grid and r['noise_sigma']==sigma and r['rate_ratio']==ratio]
                    assert len(cell)==20
                    a[i,j]=sum(r['methods'][key].get('correct_count',False) for r in cell)
            ax=axes[row,col];ax.imshow(a,vmin=0,vmax=20,cmap='cividis',aspect='auto')
            for i in range(4):
                for j in range(6): ax.text(j,i,str(int(a[i,j])),ha='center',va='center',color='white' if a[i,j]<10 else 'black',fontsize=8)
            ax.set_xticks(range(6),[str(r) for r in ratios]);ax.set_yticks(range(4),['0.0001','0.001','0.01','0.05'])
            ax.set_xlabel('Adjacent rate ratio');ax.set_ylabel('Gaussian noise SD')
            ax.set_title(f'{grid.capitalize()} grid: {short[key]} (out of 20)',fontsize=9)
        for sigma in sigmas:
            curve=[]
            for ratio in ratios:
                cell=[r for r in sweep['results'] if r['grid']==grid and r['noise_sigma']==sigma and r['rate_ratio']==ratio]
                errors=[r['methods']['varpro_known_k'].get('log_rate_errors') for r in cell]
                assert all(e is not None for e in errors)
                errors=np.array(errors); bound=np.array(cell[0]['cramer_rao']['log_rate_sd'])
                rmse=np.sqrt(np.mean(errors**2,axis=0));bias=np.mean(errors,axis=0)
                summary['cramer_rao_cells'].append({'grid':grid,'rate_ratio':ratio,'sigma':sigma,
                    'log_rate_cr_sd':bound.tolist(),'known_k_log_rate_rmse':rmse.tolist(),
                    'known_k_log_rate_bias':bias.tolist(),'rmse_over_cr_sd':(rmse/bound).tolist()})
                curve.append(rmse[1]/bound[1])
            axes[2,row].plot(range(6),curve,marker='o',ms=3,label=f'{sigma:g}')
        ax=axes[2,row];ax.axhline(1,color='gray',linestyle='--',linewidth=.8)
        ax.set_yscale('log');ax.set_xticks(range(6),[str(r) for r in ratios])
        ax.set_xlabel('Tested rate ratio');ax.set_ylabel('Known-K log-rate RMSE / CR SD')
        ax.set_title(grid.capitalize()+' grid: faster rate',fontsize=9);ax.legend(title='Noise SD',fontsize=7,title_fontsize=7,frameon=False,ncol=2)
        ax.grid(alpha=.15)
    figure(fig,'noise_sweep')

    fig,axes=plt.subplots(2,2,figsize=(8,5),gridspec_kw={'hspace':.35,'wspace':.3})
    exp_rows=[]
    for col,label in enumerate(['delayed','prompt']):
        case=next(r for r in experiment['results'] if r['label']==label and r['extra_gate_bins']==0)
        t=np.array(case['times']); y=np.array(case['observations']);scale=max(abs(y))
        axes[0,col].scatter(t,y/scale,s=5,color='#888888',alpha=.45,label='Measured tail')
        for key,color in [('native_log','#277687'),('varpro_bic','#ae6633')]:
            r=case['methods'][key]; pred=np.exp(-np.outer(t,r['rates']))@r['amplitudes']
            axes[0,col].plot(t,pred/scale,color=color,label=f"{short[key]}: {r['count']} modes")
            axes[1,col].plot(t,(pred-y)/scale,color=color,linewidth=.7,label=short[key])
        units='microseconds' if label=='delayed' else 'nanoseconds'
        axes[0,col].set_title(label.capitalize()+' fluorescence');axes[0,col].set_ylabel('Normalized intensity')
        axes[0,col].set_yscale('symlog',linthresh=.0002)
        axes[0,col].legend(fontsize=7,frameon=False)
        axes[1,col].set_ylabel('Normalized residual');axes[1,col].set_xlabel(f'Time ({units})')
        for ax in axes[:,col]:
            ax.set_xscale('log');ax.grid(alpha=.15)
            ticks=[.05,.1,.5,1,5,10] if label=='delayed' else [20,50,100,200,400]
            ax.set_xticks(ticks,[str(v) for v in ticks])
            ax.xaxis.set_minor_formatter(NullFormatter())
        for gate in [0,5]:
            other=next(r for r in experiment['results'] if r['label']==label and r['extra_gate_bins']==gate)
            a=other['methods']['native_log'];b=other['methods']['varpro_bic']
            exp_rows.append(f"{label.capitalize()}, +{gate} bins & {a['count']} & {a['observed_rmse']:.3f} & {b['count']} & {b['observed_rmse']:.3f}")
    figure(fig,'experimental_tails')
    table('experimental','Descriptive unweighted fits to measured fluorescence tails. RMSE is in original '
          'intensity units after four-bin averaging. The additional gate is a predeclared sensitivity check. '
          'Component truth is unavailable; lower residuals do not establish physically correct components.',
          'tab:experimental','lrrrr','Tail / gate & Log K & Log RMSE & VarPro K & VarPro RMSE',exp_rows)

    worst=max(legacy['results'],key=lambda r:r['methods']['native_log'].get('max_rate_relative_error') or 0)
    summary['worst_legacy_case']={key:worst[key] for key in ['id','rates','amplitudes','cramer_rao','methods']}
    summary['sources_sha256']={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
         for p in [HERE/'PROTOCOL.md',HERE/'EXPERIMENT_PROTOCOL.md',HERE/'baselines.py',HERE/'run_study.py',
                   HERE/'build_native.py',HERE/'experimental.py',HERE/'check_native.py',
                   *sorted((HERE/'results').glob('*.json'))] if p.name!='summary.json'}
    (HERE/'results/summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf8')
    report=['# Expanded research evaluation','',
        'The numerical core is frozen at release 1.0.0. These local research records do not imply a new package release.',
        'The [protocol](PROTOCOL.md) was fixed before inspecting outcomes. No core algorithm was tuned to the new sweep.','',
        '## Count recovery','',
        '| Method | Legacy correct / eligible | New sweep correct / eligible | Sweep count + accuracy pass |',
        '|---|---:|---:|---:|']
    for key in keys:
        a=summary['legacy'][key];b=summary['sweep'][key]
        report.append(f"| {LABELS[key]} | {a['correct_counts']}/{a['eligible']} | {b['correct_counts']}/{b['eligible']} | {b['passes_count_and_accuracy']}/{b['eligible']} |")
    report += ['', 'Legacy cases are reused validation/regression data. New sweep cells have 20 noise realizations each.',
        'Rate errors are conditional on correct counts; the last column counts all eligible observations, including incorrect counts as failures.',
        'Pencil methods are evaluated only on originally uniform samples. Known-K oracles are not scored as automatic count recovery.','',
        '## Matched native timing','', '| Fixture | Log ms | Linear ms | Grid ms | VarPro ms | Grid/log | Linear/log |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for r in summary['timing']:
        ms=r['ms'];report.append(f"| {r['id']} | {ms['native_log']:.3f} | {ms['native_linear']:.3f} | {ms['native_grid']:.3f} | {ms['varpro_bic']:.3f} | {r['ratios']['native_grid']:.2f} | {r['ratios']['native_linear']:.2f} |")
    report += ['', 'Seven shuffled repetitions after warmup; one BLAS thread. Every repetition and its accuracy is retained in timing.json.',
        'A ratio above one favors the log route. Fast failed fits are not accuracy-preserving speedups.',
        'The single-mode shortcut is shared. Existing pilot/holdout fixture names do not indicate untouched data.',
        'The coordinate control changes both finite-grid allocation and the geometry of smoothing while targeting the same -t*f\'(t) peaks.',
        'It is not a claim that all linear-time wavelet methods must fail or that complexity is logarithmic.','',
        '## Noise, separation, and local information','',
        'Rates are 0.5 and 0.5*r, amplitudes 1 and 0.5. Ratios 1.1, 1.25, 1.4, 2, 3, 10; noise SD 0.0001, 0.001, 0.01, 0.05.',
        'Both sampling grids use 600 observations on [0.01,20]. Seed 2026102101. All 960 outcomes are retained.',
        'CR bounds assume known count, known independent Gaussian noise, and local unbiased estimation. Amplitudes are nuisance parameters.',
        'Known-K empirical log-rate RMSE and bias use all 20 realizations per cell, without conditioning on automatic count success.',
        'A single error, a biased estimator, and a selected model are not governed directly by the displayed local standard deviations.','',
        'The legacy worst-rate-error case is '+worst['id']+'. Native and known-K variable projection reproduce the same poorly determined fit.',
        '## Experimental application','',
        'Two unquenched TCSPC tails from [Kupfer et al., Zenodo v3](https://doi.org/10.5281/zenodo.20792957), CC BY 4.0.',
        'Original files and hashes are in [data](data/README.md). [Preprocessing rules](EXPERIMENT_PROTOCOL.md) were fixed before fitting.',
        'Component counts and rates are not known. Gate sensitivity and unequal residuals are reported, not treated as recovery accuracy.',
        'Photon noise is heteroscedastic and the IRF is gated rather than convolved. This is a descriptive application outside the Gaussian validation model.',
        'The rounded experimental time columns fail the strict uniform-grid eligibility check; pencil fits are not included for those curves.','',
        '## Configuration sensitivity','']
    for key in ['public_default','native_linear_1024']:
        r=summary['legacy'][key];report.append(f"- {LABELS[key]}: {r['correct_counts']}/{r['eligible']} legacy counts. This is a secondary configuration, not a replacement for the main row.")
    report += ['', 'The 1024-point linear detector performs worse here; increasing resolution is not a recovery guarantee.',
        '## Reproducibility','',
        'CPU: Intel Core i7-7700HQ, nominal 2.80 GHz; Balanced power plan. No fixed-clock or thermal-control claim.',
        'Compiler and source hashes, Python/NumPy/SciPy versions, exceptions, warnings and raw fits are in results/.',
        'The external harness compiles a copy of cpp/src/critical.cpp and a checked coordinate transform; it does not replace the installed extension.',
        'Independent analytic-Jacobian, Fisher-information, baseline, and native-equivalence checks are included.','',
        '```powershell',r'.venv\Scripts\python benchmarks\research\build_native.py --compiler .venv\Lib\site-packages\ziglang\zig.exe',
        r'.venv\Scripts\python benchmarks\research\check_native.py',r'.venv\Scripts\python benchmarks\research\run_study.py --section all',
        r'.venv\Scripts\python benchmarks\research\experimental.py',r'.venv\Scripts\python benchmarks\research\summarize.py','```','']
    (HERE/'REPORT.md').write_text('\n'.join(report),encoding='utf8')
    print('Saved expanded report, summary, three figures, and three tables.')


if __name__=='__main__': main()
