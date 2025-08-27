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

#include <cmath>

#include "TF1.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TSpline.h"

#include "UpcCrossSection.h"
#include "UpcPhysConstants.h"
#include "UpcPhotoNuclearVM.h"

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

UpcPhotoNuclearVM::UpcPhotoNuclearVM(int partPDG, int shadowingOpt, int dghtPDG)
{
  this->partPDG = partPDG;
  this->dghtPDG = dghtPDG;
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
  if (dghtPDG == 11) {
    mDght = phys_consts::mEl;
  } else if (dghtPDG == 13) { // muon
    mDght = phys_consts::mMu;
  } else if (dghtPDG == 2212) {
    mDght = phys_consts::mProt;
  } else {
    throw std::invalid_argument("Unsupported decay mode");
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
  delete fRgJpsiLtaWeak;
  delete fRgJpsiLtaStrong;
  delete fSp3JpsiLtaWeak;
  delete fSp3JpsiLtaStrong;
}

// Modified version of Cern Library
// interpolation routine E100
void luovi(const double* f, const double* arg, int mmm, double z, double &sum)
{
  auto* cof = new double[mmm+1];
  int jndex, index;

  int mm = mmm < 20 ? mmm : 20;
  int m = mm - 1;
  for (int i = 1; i <= mm; ++i)
    cof[i] = f[i];
  for (int i = 1; i <= m; ++i) {
    for (int j = i; j <= m; ++j){
      jndex = mm - j;
      index = jndex + i;
      cof[index] = (cof[index] - cof[index-1]) / (arg[index] - arg[jndex]);
    }
  }
  sum = cof[mm];
  for (int i = 1; i <= m; ++i){
    index = mm - i;
    sum = (z - arg[index]) * sum + cof[index];
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

  int A = UpcCrossSection::A;

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
    std::string filenimi;
    if (order == 1)
      filenimi = Form("%s/vm/eps09/lo/EPS09LOR_%d", CROSS_SEC_DIR, A);
    else
      filenimi = Form("%s/vm/eps09/nlo/EPS09NLOR_%d", CROSS_SEC_DIR, A);

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
  arg[3] = Qpoint + 1;

  luovi(fg, arg, 3, realQ, res);

  double Rg  = res > 0 ? res : 0;

  return Rg;
}

double UpcPhotoNuclearVM::getRgLta(int type, double x)
{
  const int A = UpcCrossSection::A;

  static bool isInit = false;

  if (!isInit) {
    std::string csDir = CROSS_SEC_DIR;
    // using files produced by L. Frankfurt, V. Guzey and M. Strikman (Phys.Rept. 512 (2012) 255-393, arXiv: 1106.2091)
    // initialize GSL spline interpolator
    const gsl_interp2d_type* interpType = gsl_interp2d_bicubic;
    std::string fname = Form("%s/vm/lta/QCDEvolution_pb%dproton_2009_model%d.dat", CROSS_SEC_DIR, A, type == 0 ? 2 : 1);
    std::ifstream ifs(fname);
    if (!ifs) {
      PLOG_ERROR << "Missing file: " << fname;
      std::exit(1);
    }
    // read-in variables: x, Q^2, quark and gluon ratios, F_2 ratio
    double xx, q2, dv, uv, ubar, dbar, sbar, cbar, glue, f2;
    for (int i = 0; i < nq; ++i) {
      ifs >> q2;
      q2grid[i] = q2;
      for (int j = 0; j < nx; ++j) {
        ifs >> xx >> dv >> uv >> ubar >> dbar >> sbar >> cbar >> glue >> f2;
        xgrid[j] = xx;
        rgrid[i * nx + j] = glue;
      }
    }
    fSpline = gsl_spline2d_alloc(interpType, nx, nq);
    fXacc = gsl_interp_accel_alloc();
    fYacc = gsl_interp_accel_alloc();
    gsl_spline2d_init(fSpline, xgrid, q2grid, rgrid, nx, nq);
    isInit = true;
  }

  double Rg;
  if (partPDG == 443) {
    if (type == 0)
      Rg = fRgJpsiLtaWeak->Eval(x, fSp3JpsiLtaWeak);
    if (type == 1)
      Rg = fRgJpsiLtaStrong->Eval(x, fSp3JpsiLtaStrong);
  } else {
    // grid limitations
    if (x < 9.99999975e-6)
      x = 9.99999975e-6;
    if (x > 0.95)
      x = 0.95;
    Rg = gsl_spline2d_eval(fSpline, x,fMu2, fXacc, fYacc);
  }
  return Rg;
}

double UpcPhotoNuclearVM::getRgLtaVG(double x)
{
  static bool isInit = false;

  if (!isInit) {
    std::string csDir = CROSS_SEC_DIR;
    std::string fname = "";
    if (partPDG == 443) {
      fname = Form("%s/vm/lta/LT2013_pb208_cteq6l1_m12_Q2_3.dat", CROSS_SEC_DIR);
    } else if (partPDG == 100443) {
      fname = Form("%s/vm/lta/LT2013_pb208_cteq6l1_m12_Q2_4.dat", CROSS_SEC_DIR);
    } else {
      PLOG_ERROR << "PDG " << partPDG << " is not available for chosen shadowing option";
      std::exit(1);
    }
    fRgJpsiLta = new TGraph();
    std::ifstream ifs(fname);
    if (!ifs) {
      PLOG_ERROR << "Missing file: " << fname;
      std::exit(1);
    }
    double xx, rltas, rltaw;
    const int n = 37;
    for (int i = 0; i < n; ++i) {
      ifs >> xx >> rltas >> rltaw;
      fRgJpsiLta->SetPoint(i, xx, rltaw);
    }
    fSp3JpsiLta = new TSpline3("fSp3JpsiLta", fRgJpsiLta);
    fSp3JpsiLta->SetBit(TGraph::kIsSortedX);
    isInit = true;
  }

  double Rg;
  if (x > 1e-5 && x < 1e-1)
    Rg = fRgJpsiLta->Eval(x, fSp3JpsiLta);
  else
    Rg = fRgJpsiLta->Eval(x, 0, "");
  return Rg;
}


double UpcPhotoNuclearVM::calcCrossSectionY(double y)
{
  // exclusive photonuclear cross section off proton
  double w = mPart / 2. * std::exp(y); // photon energy
  double beamE = 0.5 * UpcCrossSection::sqrts; // beam energy in lab cms
  double Wgp2 = 4. * w * beamE; // photon energy in nucleon cms
  double Wgp = std::sqrt(Wgp2);
  double csGammaP = fDsDt0->Eval(Wgp);

  // integrated squared form factor
  double m2 = mPart * mPart;
  double x = m2 / Wgp2;
  double tmin = x * x * UpcCrossSection::mNucl * UpcCrossSection::mNucl;
  double tmax = tmin + 1.;
  double PhiA = fFormFactorSq->Integral(tmin, tmax);

  // ratio of scaling constants for photonuclear production off nuclei and off proton
  // V. Guzey and M. Zhalov, arXiv: 1307.4526
  double cAcP2 = 1.; // no shadowing

  // shadowing factor
  double Rg = 1.; // IA
  if (fShadowing == 1) { // EPS09 best fit
    cAcP2 = 0.97 * 0.97;
    Rg = getRgEps09(1, 1, x, fMu2);
  }
  // FGS10 parametrization works for psi(2s) and heavier
  if (fShadowing == 2 && fMu2 > 4. - 1e-6) { // LTA weak, FGS10
    cAcP2 = 0.9 * 0.9;
    Rg = getRgLta(0, x);
  }
  if (fShadowing == 3 && fMu2 > 4. - 1e-6) { // LTA strong, FGS10
    cAcP2 = 0.9 * 0.9;
    Rg = getRgLta(1, x);
  }
  // new parametrization works for jpsi and psi(2s)
  if (fShadowing == 4) { // CTEQ6 fits: see Frankfurt, Guzey, Strikman, Phys. Rept. 512 (2012) 255
    cAcP2 = 0.9 * 0.9;
    Rg = getRgLtaVG(x);
  }
  return cAcP2 * csGammaP * Rg * Rg * PhiA * 1e-6;
}
