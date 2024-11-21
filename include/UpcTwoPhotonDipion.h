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

/// a class implementing cross sections for dipion photoproduction
/// gamma+gamma -> pi0 pi0
/// disclaimer: cross section calculated by Mariola Klusek-Gawenda et al
/// see https://arxiv.org/pdf/1302.4204.pdf for more information

#ifndef UPCGENERATOR_INCLUDE_UPCTWOPHOTONDIPION_H_
#define UPCGENERATOR_INCLUDE_UPCTWOPHOTONDIPION_H_

#include "UpcElemProcess.h"

#include "TH1D.h"
#include "TH2D.h"

class UpcTwoPhotonDipion : public UpcElemProcess
{
 public:
  UpcTwoPhotonDipion(bool doMassCut, double lowMCut, double hiMCut);

  ~UpcTwoPhotonDipion() override
  {
    delete hCrossSectionM;
    delete hCrossSectionZM;
  }

  // for this process, cross sections are stored in files
  TH1D* hCrossSectionM{nullptr};
  TH2D* hCrossSectionZM{nullptr};

  double calcCrossSectionM(double m) override;

  double calcCrossSectionZM(double z, double m) override;

  // polarized cross sections are not available for this process
  double calcCrossSectionMPolS(double m) override { return 0; };

  double calcCrossSectionZMPolS(double z, double m) override { return 0; };

  double calcCrossSectionMPolPS(double m) override { return 0; };

  double calcCrossSectionZMPolPS(double z, double m) override { return 0; };
};

#endif // UPCGENERATOR_INCLUDE_UPCTWOPHOTONDIPION_H_
