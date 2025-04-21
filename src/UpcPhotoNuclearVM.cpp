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

#include <cmath>

#include "TF1.h"

#include "UpcCrossSection.h"
#include "UpcPhysConstants.h"
#include "UpcPhotoNuclearVM.h"

double formFactor(double* xx, double*)
{
  double x = xx[0];
  double ff = UpcCrossSection::calcFormFac(x);
  return ff * ff;
}

UpcPhotoNuclearVM::UpcPhotoNuclearVM(int partPDG, int shadowingOpt, double sqrts)
{
  this->partPDG = partPDG;
  fShadowing = shadowingOpt;
  fSqrts = sqrts;
  isCharged = false;
  fFormFac = new TF1("ff", formFactor, 1e-12, 100, 1);
  if (partPDG == 443) { // jpsi
    mPart = 3.0969;
    fFdsdt0 = new TF1("dsdtJpsi","x>0.9383+3.0969 ? (342*(1-(0.9383+3.0969)^2/x^2)^1.5*(x^2/10000)^0.40) : 0 ",0,300);
  } else if (partPDG == 100443) { // psi(2s)
    mPart = 3.6861;
    fFdsdt0 = new TF1("dsdtJpsi","x>0.9383+3.0969 ? (342*(1-(0.9383+3.0969)^2/x^2)^1.5*(x^2/10000)^0.40) : 0 ",0,300);
  } else if (partPDG == 553) { // upsilon(1s)
    mPart = 9.3987;
    fFdsdt0 = new TF1("dsdtUpsilon","x>0.9383+9.4603 ? (0.902*(1-(0.9383+9.4603)^2/x^2)^1.5*(x^2/10000)^0.447) : 0 ",0,1000);
  } else {
    throw std::invalid_argument("Unsupported particle");
  }
}

UpcPhotoNuclearVM::~UpcPhotoNuclearVM()
{
  delete fFdsdt0;
  delete fFormFac;
}

double UpcPhotoNuclearVM::calcCrossSectionY(double y)
{
  double mhc2 = mPart;
  double w = mhc2 / 2. * std::exp(y);
  double W = std::sqrt(w * fSqrts);
  double csGammaP = fFdsdt0->Eval(W);

  double m2 = mhc2 * mhc2;
  double Wgp2 = W * W;
  double x = m2 / Wgp2;

  double tmin = x * x / Wgp2;
//  printf("x=%.9f, wgp2=%.9f, tmin=%.9f\n", x, Wgp2, tmin);
  double PhiA = fFormFac->Integral(tmin, 100.);

  // shadowing factor
  double mu2 = 1; // 3, 4, 22.4;
  double R = 1; // (x, mu2); // IA
  if (fShadowing == 1) {
    // LTA
  }
  if (fShadowing == 2) {
    // EPS09
  }
  return csGammaP * R * R * PhiA * 1e-6;
}
