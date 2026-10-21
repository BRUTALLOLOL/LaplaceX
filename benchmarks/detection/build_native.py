"""Build before/after native controls; generated sources stay in .build/."""
import hashlib, json, os, subprocess, sys, sysconfig
from pathlib import Path
import pybind11
ROOT=Path(__file__).resolve().parents[2]
BUILD=ROOT/'.build/detection-native'; BUILD.mkdir(exist_ok=True)
before=subprocess.check_output(['git','show','3708dd2:cpp/src/critical.cpp']).decode('utf-8')
before=before.replace('namespace laplacex {','namespace detection_before {\nusing laplacex::DecomposeConfig;\nusing laplacex::ExponentialComponent;\nusing laplacex::validate_samples;',1)
(BUILD/'before.cpp').write_text(before,encoding='utf-8')
(BUILD/'after.cpp').write_bytes((ROOT/'cpp/src/critical.cpp').read_bytes())
wrapper=r'''
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "after.cpp"
#include "before.cpp"
namespace py=pybind11;
using Array=py::array_t<double,py::array::c_style|py::array::forcecast>;
std::vector<double> vec(Array a) {if(a.ndim()!=1)throw std::invalid_argument("1D required");return {a.data(),a.data()+a.size()};}
py::dict fit(Array ta,Array ya,const std::string& route,double sigma,int budget) {
 auto t=vec(ta),y=vec(ya);laplacex::DecomposeConfig c;c.noise_sigma=sigma;
 std::vector<laplacex::ExponentialComponent> result;
 {py::gil_scoped_release release;
 if(route=="before")result=detection_before::decompose(t,y,c,budget);
 else if(route=="grid_before")result=detection_before::decompose_deflate(t,y,c,budget,1e-6);
 else if(route=="grid")result=laplacex::decompose_deflate(t,y,c,budget);
 else result=laplacex::decompose(t,y,c,budget);}
 std::vector<double> rates,amps;for(auto v:result){rates.push_back(v.rate);amps.push_back(v.amplitude);}
 py::dict d;d["rates"]=rates;d["amplitudes"]=amps;return d;
}
py::dict refine(Array ta,Array ya,std::vector<double> rates) {
 auto t=vec(ta),y=vec(ya); auto f=laplacex::refine(t,y,rates,true);
 py::dict d;d["rates"]=f.rates;d["amplitudes"]=f.amplitudes;d["rss"]=f.rss;return d;
}
py::dict diagnose(Array ta,Array ya) {
 auto t=vec(ta),y=vec(ya);double scale=0;for(double v:y)scale=std::max(scale,std::abs(v));for(auto& v:y)v/=scale;
 laplacex::DecomposeConfig c;double sigma=laplacex::estimate_noise(t,y);laplacex::LogWavelets det(t,c);auto ridges=det.candidates(y,sigma);
 std::sort(ridges.begin(),ridges.end(),[](auto a,auto b){return a.votes==b.votes?a.strength>b.strength:a.votes>b.votes;});
 std::vector<double> initial;
 for(auto p:ridges){if(p.votes<8)continue;bool distinct=true;for(double r:initial)if(std::abs(std::log(r/p.rate))<.3)distinct=false;if(distinct)initial.push_back(p.rate);if(initial.size()==5)break;}
 auto a=laplacex::refine_joint(t,y,initial,true),b=laplacex::refine(t,y,initial,true);
 py::dict d;d["initial"]=initial;d["sigma"]=sigma;d["joint_rates"]=a.rates;d["projected_rates"]=b.rates;d["joint_rss"]=a.rss;d["projected_rss"]=b.rss;d["joint_amplitudes"]=a.amplitudes;d["projected_amplitudes"]=b.amplitudes;
 py::list trials;
 for(auto p:det.candidates(b.residual,sigma)){
  auto rates=b.rates;rates.push_back(p.rate);auto j=laplacex::refine_joint(t,y,rates,true),v=laplacex::refine(t,y,rates,true);
  py::dict c;c["proposal"]=p.rate;c["joint_rss"]=j.rss;c["projected_rss"]=v.rss;c["joint_rates"]=j.rates;c["projected_rates"]=v.rates;c["projected_amps"]=v.amplitudes;trials.append(c);
 }d["trials"]=trials;return d;
}
py::dict projected(Array ta,Array ya,std::vector<double> log_rates) {
 auto t=vec(ta),y=vec(ya);laplacex::Fit f;
 for(double x:log_rates)f.rates.push_back(std::exp(x));laplacex::Projection qr;
 if(!laplacex::project_fit(f,t,y,qr))throw std::invalid_argument("rank deficient");
 auto jac=laplacex::projected_jacobian(f,t,qr);
 for(auto& v:f.residual)v=-v;
 py::dict d;d["residual"]=f.residual;d["amplitudes"]=f.amplitudes;d["jacobian_columns"]=jac;return d;
}
PYBIND11_MODULE(_detection_native,m){
 m.def("fit",&fit,py::arg("t"),py::arg("y"),py::arg("route")="after",py::arg("sigma")=-1.,py::arg("budget")=5);
 m.def("refine",&refine);m.def("diagnose",&diagnose);m.def("projected",&projected);
}
'''
(BUILD/'wrapper.cpp').write_text(wrapper,encoding='utf-8')
compiler=str(ROOT/'.venv/Lib/site-packages/ziglang/zig.exe') if os.name=='nt' else 'c++'
cmd=[compiler]+(['c++','-target','x86_64-windows-gnu','-Wno-nullability-completeness'] if os.name=='nt' else [])
cmd+=['-std=c++17','-O3','-shared',str(BUILD/'wrapper.cpp'),'-I'+str(ROOT/'cpp/include'),'-I'+pybind11.get_include(),'-I'+sysconfig.get_path('include'),'-o',str(BUILD/('_detection_native'+sysconfig.get_config_var('EXT_SUFFIX')))]
if os.name=='nt':cmd+=['-L'+str(sysconfig.get_config_var('LIBDIR')),'-lpython'+str(sys.version_info.major)+str(sys.version_info.minor)]
else:
 cmd+=['-fPIC']
 if sys.platform=='darwin':cmd+=['-undefined','dynamic_lookup']
env=dict(os.environ);env['ZIG_GLOBAL_CACHE_DIR']=str(ROOT/'.build/zig-cache')
subprocess.run(cmd,check=True,env=env)
(BUILD/'build.json').write_text(json.dumps({'command':cmd,'before_commit':'3708dd2','after_source_sha256':hashlib.sha256((BUILD/'after.cpp').read_bytes()).hexdigest()},indent=2)+'\n',encoding='utf-8')
print('Built detection controls.')
