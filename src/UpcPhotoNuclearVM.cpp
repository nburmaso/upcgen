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

#ifndef CROSS_SEC_DIR
#define CROSS_SEC_DIR "undefined"
#endif

double formFactorSq(double* xx, double*)
{
  double x = xx[0];
  double ff = UpcCrossSection::calcFormFac(x);
  return ff * ff;
}

double dsdt(double* xx, double* par)
{
  double Wgp = xx[0];
  double Wgp2 = Wgp * Wgp;
  double mmin = par[0];
  double c0 = par[1];
  double pw = par[2];
  double res = 0.;
  if (Wgp > mmin) {
    res = c0 * std::pow(1. - (mmin * mmin) / Wgp2, 1.5) * std::pow(Wgp2 * 1e-4, pw);
  }
  return res;
}

UpcPhotoNuclearVM::UpcPhotoNuclearVM(int partPDG, int shadowingOpt)
{
  this->partPDG = partPDG;
  fShadowing = shadowingOpt;
  isCharged = false;
  fFormFactorSq = new TF1("fFormFactorSq", formFactorSq, 1e-12, 100, 1);
  fDsDt0 = new TF1("fDsDt0", dsdt, 0., 1000., 3);
  double c0;
  double pw = 0.4;
  if (partPDG == 443) { // jpsi
    mPart = 3.0969;
    c0 = 342.;
    fMu2 = 3.;
  } else if (partPDG == 100443) { // psi(2s)
    mPart = 3.6861;
    c0 = 56.8;
    fMu2 = 4.;
  } else if (partPDG == 553) { // upsilon(1s)
    mPart = 9.3987;
    c0 = 0.902;
    pw = 0.447;
    fMu2 = 22.4;
  } else {
    throw std::invalid_argument("Unsupported particle");
  }
  double mmin = phys_consts::mProt + mPart;
  fDsDt0->FixParameter(0, mmin);
  fDsDt0->FixParameter(1, c0);
  fDsDt0->FixParameter(2, pw);
}

UpcPhotoNuclearVM::~UpcPhotoNuclearVM()
{
  delete fFormFactorSq;
  delete fDsDt0;
}

// Modified version of Cern Library
// interpolation routine E100
void luovi(const double* f, const double* arg, int mmm, double z, double &sum)
{
  auto* cof = new double[mmm+1];
  int jndex, index;

  int mm = mmm < 20 ? mmm : 20;
  int m = mm - 1;
  for (int i=1; i<=mm; i++) cof[i] = f[i];
  for (int i=1; i<=m; i++) {
    for (int j=i; j<=m; j++){
      jndex = mm - j;
      index = jndex + i;
      cof[index] = (cof[index]-cof[index-1])/(arg[index]-arg[jndex]);
    }
  }
  sum = cof[mm];
  for (int i=1; i<=m; i++){
    index = mm - i;
    sum = (z-arg[index])*sum + cof[index];
  }
  delete[] cof;
}

/// original reader by Thomas Ullrich (thomas.ullrich@bnl.gov)
/// https://research.hip.fi/qcdtheory/nuclear-pdfs/eps09/
/// order: 1=LO, 2=NLO   ; integer
/// pset : 1...31        ; integer
///        1     = central (best) fit
///        2,3   = error sets S{+1}, S{-1}
///        4,5   = error sets S{+2}, S{-2}
///        ...   ...
///        30,31 = error sets {S+15}, {S-15}
/// A    : atomic number  ; integer
/// x    : Bjorken-x      ; double precision
/// Q2   : scale in GeV^2 ; double precision
double UpcPhotoNuclearVM::getRgEps09(int order, int pset, double x, double Q2)
{
  const double Q2min = 1.69;
  const double Q2max = 1000000;
  const double Qsteps = 50;

  int A  = UpcCrossSection::A;

  double x_i = 0.000001, arg[4+1], fu[4+1], res, fg[3+1];

  int xlinsteps = 25;
  int xlogsteps = 25;

  static double allvalues[31+1][8+1][50+1][50+1];
  static int psetlast = -1, Alast = -1, orderlast = -1;

  // Freeze x if it's < 10E-6 or > 1
  // Freeze Q^2 if it's < 1.69 or > 10E+6

  if (x < x_i) x = x_i;
  if (x > 1) x = 1;
  if (Q2 < Q2min) Q2 = Q2min;
  if (Q2 > Q2max) Q2 = Q2max;

  // If the set specifications have been changed, read the tables again
  if (A != Alast || order != orderlast) {
    std::string csDir = CROSS_SEC_DIR;
    std::string filenimi = csDir;
    if (order == 1)
      filenimi += "/vm/eps09/lo/EPS09LOR_" + std::to_string(A);
    else
      filenimi += "/vm/eps09/nlo/EPS09NLOR_" + std::to_string(A);

    std::ifstream ifs(filenimi);

    if (!ifs) {
      PLOG_ERROR << "Missing file: " << filenimi;
      std::exit(1);
    }

    std::string dummy;
    for (int setnumber = 1; setnumber <= 31; setnumber++) {
      for (int k = 0; k <= 50; ++k) {
        ifs >> dummy;
        for (int t = 0; t <= 50; ++t)
          for (int p = 1; p <= 8; ++p)
            ifs >> allvalues[setnumber][p][k][t];
      }
    }

    ifs.close();

    psetlast  = pset;
    Alast     = A;
    orderlast = order;
  }

  // Find out the position in the loglog Q^2-grid
  double realQ  = Qsteps * (std::log(std::log(Q2)/std::log(Q2min))) / (std::log(std::log(Q2max)/std::log(Q2min)));
  int Qpoint = static_cast<int>(realQ);
  if (Qpoint <= 0) Qpoint = 1;
  if (Qpoint >= static_cast<int>(Qsteps + 0.5) - 1) Qpoint = static_cast<int>(Qsteps + 0.5) - 1;

  double LSTEP = (1. / xlogsteps) * std::log(0.1 / x_i);

  // Interpolate the grids
  double result;
  int t = 8; // skipping everything except Rg

  // Find the position in the x-grid
  double n_x;
  if (x <= 0.1) n_x = ((1. / LSTEP) * log(x / x_i));
  else n_x = ((x - 0.1) * xlinsteps / (1. - 0.1) + xlogsteps);
  int xpoint = static_cast<int>(n_x);

  if (xpoint <= 0) xpoint = 1;
  if (xpoint >= (xlinsteps + xlogsteps) - 4) xpoint = (xlinsteps + xlogsteps) - 4;

  for (int k = 1; k <= 4; ++k) {
    if (xpoint-2+k < xlogsteps)
      arg[k] = x_i * exp(LSTEP * (xpoint - 2 + k));
    else
      arg[k] = 0.1 + (xpoint - 2 + k - xlogsteps) * (1 - 0.1) / xlinsteps;
  }

  for (int j= 1; j <= 3; ++j) {
    fu[1] = allvalues[pset][t][Qpoint - 2 + j][xpoint - 1];
    fu[2] = allvalues[pset][t][Qpoint - 2 + j][xpoint];
    fu[3] = allvalues[pset][t][Qpoint - 2 + j][xpoint + 1];
    fu[4] = allvalues[pset][t][Qpoint - 2 + j][xpoint + 2];
    luovi(fu, arg, 4, x, res);
    fg[j] = res;
  }

  arg[1] = Qpoint - 1;
  arg[2] = Qpoint;
  arg[3] = Qpoint+1;

  luovi(fg, arg, 3, realQ, res);

  double Rg  = res > 0 ? res : 0;

  return Rg;
}

double UpcPhotoNuclearVM::calcCrossSectionY(double y)
{
  // exclusive photonuclear cross section off proton
  double w = mPart / 2. * std::exp(y); // photon energy
  double beamE = 0.5 * UpcCrossSection::sqrts; // beam energy in lab cms
  double Wgp = std::sqrt(4. * w * beamE); // photon energy in nucleon cms
  double csGammaP = fDsDt0->Eval(Wgp);

  // integrated squared form factor
  double m2 = mPart * mPart;
  Double_t tmin = (0.0625 * m2 * m2 / (w * w * UpcCrossSection::g1 * UpcCrossSection::g1));
  double PhiA = fFormFactorSq->Integral(tmin, 100.);

  // shadowing factor
  double Wgp2 = Wgp * Wgp;
  double x = m2 / Wgp2;
  double Rg = 1; // (x, mu2); // IA
  if (fShadowing == 1) {
    // LTA
  }
  if (fShadowing == 2) {
    Rg = getRgEps09(1, 1, x, fMu2);
  }
  return csGammaP * Rg * Rg * PhiA * 1e-6;
}
