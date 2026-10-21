"""Compile external matched controls without changing the installed library.

Windows: python benchmarks/research/build_native.py --compiler path/to/zig.exe
Other platforms: --compiler c++ (requires Python development libraries).
Generated copies and binaries live in the ignored .build/research-native folder.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import sysconfig

import pybind11

ROOT = Path(__file__).resolve().parents[2]
BUILD = ROOT / '.build/research-native'


def replace_once(source, old, new):
    assert source.count(old) == 1, old
    return source.replace(old, new)


def linear_source(source):
    source = replace_once(source, 'namespace laplacex {',
                          'namespace linear_control {\nusing laplacex::DecomposeConfig;\n'
                          'using laplacex::ExponentialComponent;\nusing laplacex::validate_samples;')
    source = replace_once(source, 'Vector x, alpha;', 'Vector x, alpha;\n    Vector times;\n    double slope, offset;')
    source = replace_once(source, 'x.resize(n); alpha.resize(n); left.resize(n);',
        'x.resize(n); alpha.resize(n); left.resize(n); times.resize(n);\n'
        '        slope = (std::log(t.back())-std::log(t.front()))/(t.back()-t.front());\n'
        '        offset = slope*t.front()-std::log(t.front());')
    source = replace_once(source, 'double ti = std::exp(x[i]);',
                          'double ti = t.front()+i*(t.back()-t.front())/(n-1);\n        times[i] = ti;')
    source = replace_once(source, 'double maximum = 0;\n            for (double value : derivative)',
        'for (int i = 0; i < n; ++i) {\n'
        '                double coordinate = slope*times[i];\n'
        '                mexican[i] = kernel.scale*kernel.scale*derivative[i]+coordinate*mexican[i];\n'
        '                derivative[i] *= coordinate;\n'
        '            }\n'
        '            double maximum = 0;\n            for (double value : derivative)')
    source = replace_once(source, 'cfg.snr_threshold*sigma*noise_norm(kernel, i)',
                          'cfg.snr_threshold*sigma*noise_norm(kernel, i)*slope*times[i]')
    source = replace_once(source, 'double log_rate = -(x[i]+fraction*dx)-.5*kernel.scale*kernel.scale;',
        'double peak_time = times[i]+fraction*(times[i+1]-times[i]);\n'
        '                double log_rate = -std::log(peak_time);')
    return source


WRAPPER = r'''
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "original.cpp"
#include "linear.cpp"
namespace py = pybind11;
using Array = py::array_t<double, py::array::c_style | py::array::forcecast>;
std::vector<double> vec(Array a) {
  if(a.ndim()!=1) throw std::invalid_argument("one-dimensional samples required");
  return {a.data(), a.data()+a.size()};
}
py::dict fit(Array ta, Array ya, const std::string& route, int points, double sigma, int budget) {
  auto t=vec(ta), y=vec(ya); laplacex::DecomposeConfig c;
  c.n_log_points=points; c.noise_sigma=sigma;
  std::vector<laplacex::ExponentialComponent> result;
  { py::gil_scoped_release release;
    if(route=="linear") result=linear_control::decompose(t,y,c,budget);
    else if(route=="grid") result=laplacex::decompose_deflate(t,y,c,budget);
    else if(route=="log") result=laplacex::decompose(t,y,c,budget);
    else throw std::invalid_argument("unknown route");
  }
  std::vector<double> rates, amplitudes;
  for(const auto& comp:result) { rates.push_back(comp.rate); amplitudes.push_back(comp.amplitude); }
  py::dict out; out["rates"]=rates; out["amplitudes"]=amplitudes; return out;
}
py::dict proposals(Array ta, Array ya, const std::string& route, int points) {
  auto t=vec(ta), y=vec(ya); laplacex::DecomposeConfig c; c.n_log_points=points;
  std::vector<double> rates; std::vector<int> votes;
  if(route=="linear") {
    linear_control::LogWavelets d(t,c);
    for(auto p:d.candidates(y,0)) { rates.push_back(p.rate); votes.push_back(p.votes); }
  } else {
    laplacex::LogWavelets d(t,c);
    for(auto p:d.candidates(y,0)) { rates.push_back(p.rate); votes.push_back(p.votes); }
  }
  py::dict out; out["rates"]=rates; out["votes"]=votes; return out;
}
PYBIND11_MODULE(_research_native,m) {
  m.def("fit",&fit,py::arg("t"),py::arg("y"),py::arg("route")="log",
        py::arg("points")=256,py::arg("sigma")=-1.,py::arg("budget")=5);
  m.def("proposals",&proposals,py::arg("t"),py::arg("y"),py::arg("route")="log",py::arg("points")=256);
}
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', required=True)
    args = parser.parse_args()
    BUILD.mkdir(parents=True, exist_ok=True)
    original = (ROOT/'cpp/src/critical.cpp').read_text(encoding='utf8')
    for name, text in [('original.cpp',original),('linear.cpp',linear_source(original)),('wrapper.cpp',WRAPPER)]:
        (BUILD/name).write_text(text,encoding='utf8')
    output = BUILD/('_research_native'+sysconfig.get_config_var('EXT_SUFFIX'))
    cmd=[args.compiler]
    if Path(args.compiler).stem.lower()=='zig': cmd += ['c++','-target','x86_64-windows-gnu','-Wno-nullability-completeness']
    cmd += ['-std=c++17','-O3','-shared',str(BUILD/'wrapper.cpp'),
            '-I'+str(ROOT/'cpp/include'), '-I'+pybind11.get_include(),
            '-I'+sysconfig.get_path('include'),'-o',str(output)]
    if os.name=='nt': cmd += ['-L'+str(sysconfig.get_config_var('LIBDIR')),'-lpython'+str(sys.version_info.major)+str(sys.version_info.minor)]
    else:
        cmd += ['-fPIC']
        if sys.platform=='darwin': cmd += ['-undefined','dynamic_lookup']
    env=dict(os.environ); env['ZIG_GLOBAL_CACHE_DIR']=str(ROOT/'.build/zig-cache')
    subprocess.run(cmd,check=True,cwd=ROOT,env=env)
    meta={'command':cmd, 'critical_cpp_sha256':hashlib.sha256((ROOT/'cpp/src/critical.cpp').read_bytes()).hexdigest(),
          'linear_copy_sha256':hashlib.sha256((BUILD/'linear.cpp').read_bytes()).hexdigest(),
          'binary_sha256':hashlib.sha256(output.read_bytes()).hexdigest()}
    (BUILD/'build.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf8')
    print(output)


if __name__=='__main__': main()
