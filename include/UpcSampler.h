//////////////////////////////////////////////////////////////////////////
// Copyright (C) 2021-2024, Nazar Burmasov, Evgeny Kryshen
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

// simple interface class for sampling from root histograms using stl
// can work with a 2d distributions in rectangular region

#ifndef UPCGENERATOR_INCLUDE_UPCSAMPLER_H_
#define UPCGENERATOR_INCLUDE_UPCSAMPLER_H_

#include <map>
#include <random>
#include <vector>

#include "TH1D.h"
#include "TH2D.h"

class UpcSampler
{
 public:
  UpcSampler() = default;

  ~UpcSampler() = default;

  std::mt19937_64 randGen;

  // 1d distribution for y variable (usually independent of x)
  std::discrete_distribution<int>* distr1D;
  std::vector<double>* yGrid;

  // vector of 1d distributions along x for each y
  std::vector<std::discrete_distribution<int>>* distr2D;
  std::vector<double>* xGrid;

  void setYGrid(std::vector<double>& yGrid)
  {
    this->yGrid = &yGrid;
  }

  void setXGrid(std::vector<double>& xGrid)
  {
    this->xGrid = &xGrid;
  }

  // set seed for internal generator
  void setSeed(long int seed)
  {
    std::random_device rd;
    randGen = std::mt19937_64(seed == 0 ? rd() : seed);
  }

  // set 1d distribution from 1d histogram
  void setDistr1D(const TH1D* hist1D, int ny)
  {
    std::vector<double> distr(ny, 0.);
    for (int i = 1; i <= ny; i++) {
      distr[i - 1] = hist1D->GetBinContent(i);
    }
    distr1D = new std::discrete_distribution<int>(distr.begin(), distr.end());
  }

  // set vector of 1d distributions from 2d histogram
  void setDistr2D(const TH2D* hist2D, int nx, int ny)
  {
    distr2D = new std::vector<std::discrete_distribution<int>>(ny);
    std::vector<double> distr(nx, 0.);
    for (int j = 1; j <= ny; j++) {
      for (int i = 1; i <= nx; i++) {
        distr[i - 1] = 0;
        distr[i - 1] = hist2D->GetBinContent(i, j);
      }
      (*distr2D)[j - 1] = std::discrete_distribution<int>(distr.begin(), distr.end());
    }
  }

  // sample y value
  void sample1D(double& y)
  {
    y = (*yGrid)[(*distr1D)(randGen)];
  }

  // sample y bin
  void sample1D(int& iy)
  {
    iy = (*distr1D)(randGen);
  }

  // sample y and x values
  // first, y is sampled independently,
  // then x is sampled for the obtained y
  void sample2D(double& y, double& x)
  {
    int iy = (*distr1D)(randGen);
    y = (*yGrid)[iy];
    int ix = (*distr2D)[iy](randGen);
    x = (*xGrid)[ix];
  }

  // sample y and x bins -- same method as for the values
  void sample2D(int& iy, int& ix)
  {
    iy = (*distr1D)(randGen);
    ix = (*distr2D)[iy](randGen);
  }

  // sample x for a given y
  void sampleXAtY(int iy, double& x)
  {
    int ix = (*distr2D)[iy](randGen);
    x = (*xGrid)[ix];
  }
};

#endif // UPCGENERATOR_INCLUDE_UPCSAMPLER_H_
