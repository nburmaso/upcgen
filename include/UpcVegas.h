//////////////////////////////////////////////////////////////////////////
// Copyright (C) 2021-2025, Nazar Burmasov, Evgeny Kryshen
//
// E-mail of the corresponding author: nazar.burmasov@cern.ch
//
// This file is a part of Upcgen
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//////////////////////////////////////////////////////////////////////////

// Simplistic interface to the newest VEGAS by G.P.Lepage
// (G. P. Lepage, J. Comput. Phys. 27(1978) 192, J. Comput. Phys. 439 (2021) 110386)
// via pybind11 (https://github.com/pybind/pybind11)

#pragma once

#include <cassert>
#include <chrono>
#include <functional>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>

#include "plog/Appenders/ColorConsoleAppender.h"
#include "plog/Formatters/TxtFormatter.h"
#include "plog/Init.h"
#include "plog/Initializers/RollingFileInitializer.h"
#include "plog/Log.h"

namespace py = pybind11;
using namespace py::literals;
using namespace std::chrono;

class UpcVegas
{
 public:
  UpcVegas()
  {
    py::initialize_interpreter();
  }

  ~UpcVegas()
  {
    // py::finalize_interpreter();
  }

  // initialize integrator and perform warm-up
  void init()
  {
    py::object gvar_ranseed = py::module_::import("gvar").attr("ranseed");
    if (seed == 0) {
      seed = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    }
    gvar_ranseed(seed);
    py::module_ vegas = py::module_::import("vegas");
    py::object Integrator = py::module_::import("vegas").attr("Integrator");
    py::list list = py::cast(xlims);
    integrator = std::make_shared<py::object>(Integrator(list));
    // py::print(vegas);
    // py::print(integrator.get());
    func_py = std::make_shared<py::cpp_function>(func_cpp);
    PLOG_INFO << "VEGAS warm-up";
    (*integrator)(*func_py,
                  "nitn"_a=20, "neval"_a=10000, "alpha"_a=0.25, "adapt"_a=true);
//                  "nitn"_a=20, "neval"_a=1000000, "alpha"_a=0.25, "adapt"_a=true);
  }

  // actual integration
  double integrate()
  {
    PLOG_INFO << "Integrating with VEGAS";
    double res, chi2, dof;
    py::object result = (*integrator)(*func_py,
                                      "nitn"_a=100, "neval"_a=10000, "alpha"_a=0.05, "adapt"_a=true);
//                                      "nitn"_a=100, "neval"_a=100000, "alpha"_a=0.05, "adapt"_a=true);
    py::print(py::getattr(result,"summary")());
    res = py::cast<double>(py::getattr(result, "mean"));
    chi2 = py::cast<double>(py::getattr(result, "chi2"));
    dof = py::cast<double>(py::getattr(result, "dof"));
    PLOG_INFO << std::setprecision(6) << "Integrating with VEGAS: res=" << res << ", chi2/dof=" << chi2 / dof;
    return res;
  }

  // generate sample with shape of [sampleNEvts, ndim]
  auto generateSample()
  {
    if (sample != nullptr)
      py::getattr(*sample, "clear")();
    sample = std::make_shared<py::list>(integrator->attr("random")());
    sampleEvId = 0;
    sampleNEvts = static_cast<int>(py::len(*sample));
    // printf("generated sample of [%zd]\n", py::len(*sample));
  }

  auto generate()
  {
    if (sampleEvId == sampleNEvts) {
      generateSample();
    }
    py::tuple t = (*sample)[sampleEvId];
    auto x = py::cast<std::vector<double>>(t[0]);
    sampleEvId++;
    return x;
  }

  uint64_t seed{1};
  std::vector<std::vector<double>> xlims{};
  std::shared_ptr<py::object> integrator{nullptr};
  std::function<double(std::vector<double>&)> func_cpp{};
  std::shared_ptr<py::cpp_function> func_py{nullptr};
  std::shared_ptr<py::list> sample{nullptr};
  int sampleNEvts{0};
  int sampleEvId{0};
};
