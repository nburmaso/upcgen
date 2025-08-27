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

/// a class implementing photonuclear cross sections

#pragma once

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "UpcElemProcess.h"

class TGraph;
class TSpline3;
class UpcCrossSection;

class UpcPhotoNuclearVM : public UpcElemProcess
{
 public:
  explicit UpcPhotoNuclearVM(int partPDG, int shadowingOpt, int dghtPDG);

  ~UpcPhotoNuclearVM() override;

  double calcCrossSectionY(double m) override;

  double calcCrossSectionZM(double z, double m) override { return 0.; }

  double calcCrossSectionMPolS(double m) override { return 0.; }

  double calcCrossSectionZMPolS(double z, double m) override { return 0.; }

  double calcCrossSectionMPolPS(double m) override { return 0.; }

  double calcCrossSectionZMPolPS(double z, double m) override { return 0.; }

 private:

  int fShadowing{0};           // gluon shadowing calculation: 0 = IA, 1 = EPS09, 2 = LTA weak, 3 = LTA strong, 4 = LTA (jpsi and psip)
  double fMu2{1.};             // resolution scale (GeV^2)
  TF1* fFormFactorSq{nullptr}; // squared form factor of a nuclei
  TF1* fDsDt0{nullptr};        // cross section for photonuclear VM production off proton

  // Rg interpolation helpers
  TGraph* fRgJpsiLta{nullptr};
  TGraph* fRgJpsiLtaWeak{nullptr};
  TGraph* fRgJpsiLtaStrong{nullptr};
  TSpline3* fSp3JpsiLta{nullptr};
  TSpline3* fSp3JpsiLtaWeak{nullptr};
  TSpline3* fSp3JpsiLtaStrong{nullptr};

  // gsl interpolation helpers
  const int nx = 90;
  const int nq = 7;
  double xgrid[90];
  double q2grid[7];
  double rgrid[7 * 90];
  gsl_interp_accel* fXacc{nullptr};
  gsl_interp_accel* fYacc{nullptr};
  gsl_spline2d* fSpline{nullptr};

  double getRgEps09(int order, int pset, double x, double Q2);
  double getRgLta(int type, double x);
  double getRgLtaVG(double x);
};
