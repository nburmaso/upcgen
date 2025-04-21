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

/// a class implementing photonuclear cross sections

#pragma once

#include "UpcElemProcess.h"

class UpcCrossSection;

class UpcPhotoNuclearVM : public UpcElemProcess
{
 public:
  explicit UpcPhotoNuclearVM(int partPDG, int shadowingOpt);

  ~UpcPhotoNuclearVM() override;

  double calcCrossSectionY(double m) override;

  double calcCrossSectionZM(double z, double m) override { return 0.; }

  double calcCrossSectionMPolS(double m) override { return 0.; }

  double calcCrossSectionZMPolS(double z, double m) override { return 0.; }

  double calcCrossSectionMPolPS(double m) override { return 0.; }

  double calcCrossSectionZMPolPS(double z, double m) override { return 0.; }

 private:

  int fShadowing{0}; // gluon shadowing calculation: 0 = IA, 1 = LTA, 2 = EPS09
  double fMu2{1.}; // resolution scale (GeV^2)
  TF1* fFormFactorSq{nullptr}; // squared form factor of a nuclei
  TF1* fDsDt0{nullptr}; // cross section for photonuclear VM production off proton

  double getRgEps09(int order, int pset, double x, double Q2);
};
