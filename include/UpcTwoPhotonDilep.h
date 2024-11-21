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

/// a class implementing cross sections for dilepton photoproduction
/// gamma+gamma -> l+l-

#ifndef UPCGENERATOR_INCLUDE_UPCTWOPHOTONDILEP_H_
#define UPCGENERATOR_INCLUDE_UPCTWOPHOTONDILEP_H_

#include "UpcElemProcess.h"

class UpcTwoPhotonDilep : public UpcElemProcess
{
 public:
  explicit UpcTwoPhotonDilep(int partPDG);

  ~UpcTwoPhotonDilep() override = default;

  // anomalous magnetic moment
  double aLep{0};

  double calcCrossSectionM(double m) override;

  double calcCrossSectionZM(double z, double m) override;

  double calcCrossSectionMPolS(double m) override;

  double calcCrossSectionZMPolS(double z, double m) override;

  double calcCrossSectionMPolPS(double m) override;

  double calcCrossSectionZMPolPS(double z, double m) override;
};

#endif // UPCGENERATOR_INCLUDE_UPCTWOPHOTONDILEP_H_
