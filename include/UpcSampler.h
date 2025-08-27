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

#pragma once

#include <algorithm>
#include <cassert>
#include <chrono>
#include <map>
#include <random>

#include <vector>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_rng.h>

using namespace std::chrono;

class UpcSampler1D
{
 public:
  UpcSampler1D(const std::vector<double>& dist,
               const std::vector<double>& binEdges,
               uint64_t seed = std::numeric_limits<uint64_t>::max())
  {
    auto nBinsX = static_cast<int>(dist.size());
    auto nBinEdgesX = static_cast<int>(binEdges.size());
    hist = gsl_histogram_alloc(dist.size());
    gsl_histogram_set_ranges(hist, binEdges.data(), nBinEdgesX);
    rng = gsl_rng_alloc(gsl_rng_default);
    if (seed == std::numeric_limits<uint64_t>::max()) {
      seed = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    }
    gsl_rng_set(rng, seed);
    for (int i = 0; i < nBinsX; ++i) {
      hist->bin[i] = dist[i];
    }
    hpdf = gsl_histogram_pdf_alloc(nBinsX);
    int status = gsl_histogram_pdf_init(hpdf, hist);
    assert(("Could not initialize 1D sampler", status != GSL_EDOM));
  }

  ~UpcSampler1D()
  {
    gsl_rng_free(rng);
    gsl_histogram_free(hist);
    gsl_histogram_pdf_free(hpdf);
  }

  double operator()() const
  {
    return gsl_histogram_pdf_sample(hpdf, gsl_rng_uniform(rng));
  }

  gsl_rng* rng{nullptr};            // random generator
  gsl_histogram* hist{nullptr};     // input histogram
  gsl_histogram_pdf* hpdf{nullptr}; // pdf to be sampled from
};

class UpcSampler2D
{
 public:
  UpcSampler2D(const std::vector<std::vector<double>>& dist,
               const std::vector<double>& binEdgesX,
               const std::vector<double>& binEdgesY,
               uint64_t seed = std::numeric_limits<uint64_t>::max())
  {
    nBinsX = static_cast<int>(dist.size());
    nBinsY = static_cast<int>(dist[0].size());
    auto nBinEdgesX = static_cast<int>(binEdgesX.size());
    auto nBinEdgesY = static_cast<int>(binEdgesY.size());
    edgesX = binEdgesX;
    edgesY = binEdgesY;
    hist = gsl_histogram2d_alloc(nBinsX, nBinsY);
    gsl_histogram2d_set_ranges(hist,
                               binEdgesX.data(), nBinEdgesX,
                               binEdgesY.data(), nBinEdgesY);
    rng = gsl_rng_alloc(gsl_rng_default);
    if (seed == std::numeric_limits<uint64_t>::max()) {
      seed = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    }
    gsl_rng_set(rng, seed);
    for (int i = 0; i < nBinsX; ++i) {
      for (int j = 0; j < nBinsY; ++j) {
        hist->bin[i * nBinsY + j] = dist[i][j]; // inside hist bins are stored as single 1D array
      }
    }
    hpdf = gsl_histogram2d_pdf_alloc(nBinsX, nBinsY);
    int status = gsl_histogram2d_pdf_init(hpdf, hist);
    assert(("Could not initialize 2D sampler", status != GSL_EDOM));
  }

  ~UpcSampler2D()
  {
    gsl_rng_free(rng);
    gsl_histogram2d_free(hist);
    gsl_histogram2d_pdf_free(hpdf);
  }

  void operator()(double& x, double& y) const
  {
    gsl_histogram2d_pdf_sample(hpdf, gsl_rng_uniform(rng), gsl_rng_uniform(rng), &x, &y);
  }

  // from 0 to nBinsX
  int getBinX(double x)
  {
    return int(nBinsX * (x - edgesX.front())) / (edgesX.back() - edgesX.front());
  }

  // from 0 to nBinsY
  int getBinY(double y)
  {
    return int(nBinsY * (y - edgesY.front())) / (edgesY.back() - edgesY.front());
  }

  int nBinsX;
  int nBinsY;
  std::vector<double> edgesX;
  std::vector<double> edgesY;
  gsl_rng* rng;                       // random generator
  gsl_histogram2d* hist{nullptr};     // input histogram
  gsl_histogram2d_pdf* hpdf{nullptr}; // pdf to be sampled from
};
