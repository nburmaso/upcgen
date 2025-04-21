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
  explicit UpcPhotoNuclearVM(int partPDG, int shadowingOpt, double sqrts);

  ~UpcPhotoNuclearVM() override;

  double calcCrossSectionY(double m) override;

  double calcCrossSectionZM(double z, double m) override { return 0.; }

  double calcCrossSectionMPolS(double m) override { return 0.; }

  double calcCrossSectionZMPolS(double z, double m) override { return 0.; }

  double calcCrossSectionMPolPS(double m) override { return 0.; }

  double calcCrossSectionZMPolPS(double z, double m) override { return 0.; }

 private:

  int fShadowing{0};
  TF1* fFdsdt0{nullptr};
  TF1* fFormFac{nullptr};
  double fSqrts{5020.};
};
