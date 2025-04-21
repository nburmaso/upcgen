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

#include "UpcPhotoNuclearVM.h"
#include "UpcPhysConstants.h"

UpcPhotoNuclearVM::UpcPhotoNuclearVM(int partPDG)
{
  this->partPDG = partPDG;
  // populating lepton mass map
  // data from PDG
  if (partPDG == 443) { // jpsi
    mPart = 3.0969;
    isCharged = false;
  }
}

double UpcPhotoNuclearVM::calcCrossSectionY(double y, int shadowing_option)
{
  // gamma-p cross section for J/psi, psi(2S) and Upsilon (dsigma/dt|t=0)
  double csGammaP = 3.14;
  // scale
  double mu2 = 1; // 3, 4, 22.4;
  // shadowing factor
  double m2 = mPart*mPart;
  double sqrts = 1;
  double Wgp2 = sqrts*mPart*exp(-y);
  double x = m2/Wgp2;
  double R = 1; // (x, mu2);
  

  return csGammaP*R*R;
}
