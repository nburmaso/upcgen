//////////////////////////////////////////////////////////////////////////
// Copyright (C) 2021, Nazar Burmasov, Evgeny Kryshen
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

/// disclaimer: cross section calculated by Mariola Klusek-Gawenda et al
/// see https://arxiv.org/pdf/1302.4204.pdf for more information

#include "UpcTwoPhotonDipion.h"

#include "TFile.h"

UpcTwoPhotonDipion::UpcTwoPhotonDipion()
{
  mPart = 0.1349770; // pi0 mass from PDG
  partPDG = 111;
  isCharged = false;

  // load cross sections from files
  auto* f_m = new TFile("../cross_sections/pi0pi0/cross_section_m.root");
  hCrossSectionM = (TH1D*)f_m->Get("hCrossSectionM");
  hCrossSectionM->SetDirectory(nullptr);
  f_m->Close();

  auto* f_zm = new TFile("../cross_sections/pi0pi0/cross_section_zm.root");
  hCrossSectionZM = (TH2D*)f_zm->Get("hCrossSectionZM");
  hCrossSectionZM->SetDirectory(nullptr);
  f_zm->Close();
}

double UpcTwoPhotonDipion::calcCrossSectionM(double m)
{
  double cs = hCrossSectionM->Interpolate(m);
  return cs;
  // [nb]
}

double UpcTwoPhotonDipion::calcCrossSectionZM(double z, double m)
{
  double cs = hCrossSectionZM->Interpolate(z, m);
  return cs;
}
