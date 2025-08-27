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

#ifndef UPCGENERATOR_INCLUDE_UPCTWOPHOTONALP_H_
#define UPCGENERATOR_INCLUDE_UPCTWOPHOTONALP_H_

#include "UpcElemProcess.h"

class UpcTwoPhotonALP : public UpcElemProcess
{
 public:
  UpcTwoPhotonALP(double mass, double width, int spin = 0) : width{width}, spin{spin}
  {
    partPDG = 51; // code for spin-0 axion-like particle according to PDG MC numbering
    mPart = mass;
    isCharged = false;
  }

  ~UpcTwoPhotonALP() override = default;

  double width{0.001};

  int spin{0};

  double calcCrossSectionM(double m) override;

  double calcCrossSectionZM(double z, double m) override { return 0; };

  // polarized cross sections are not available for this process
  double calcCrossSectionMPolS(double m) override { return 0; };

  double calcCrossSectionZMPolS(double z, double m) override { return 0; };

  double calcCrossSectionMPolPS(double m) override { return 0; };

  double calcCrossSectionZMPolPS(double z, double m) override { return 0; };
};

#endif // UPCGENERATOR_INCLUDE_UPCTWOPHOTONALP_H_
