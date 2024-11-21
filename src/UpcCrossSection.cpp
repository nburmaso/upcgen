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

#include "UpcCrossSection.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace std;

// out-of-line initialization for static members
double UpcCrossSection::rho0 = 0;
double UpcCrossSection::R = 6.68;
double UpcCrossSection::a = 0.447;
int UpcCrossSection::Z = 82;
double UpcCrossSection::sqrts = 5020;
double UpcCrossSection::g1 = sqrts / (2. * phys_consts::mProt);
double UpcCrossSection::g2 = sqrts / (2. * phys_consts::mProt);
int UpcCrossSection::debug = 0;
double* UpcCrossSection::vCachedFormFac = new double[UpcCrossSection::nQ2];

UpcCrossSection::UpcCrossSection()
{
  constexpr int nbc = 10000000;
  vCachedBreakup = new double[nbc];
  gtot = TMath::CosH((TMath::ACosH(g1) + TMath::ACosH(g2)) / 2.);
}

UpcCrossSection::~UpcCrossSection()
{
  delete[] vCachedBreakup;
  delete elemProcess;
}

void UpcCrossSection::setElemProcess(int procID)
{
  switch (procID) {
    case 22: { // light-by-light
      elemProcess = new UpcTwoPhotonLbyL(doMassCut, lowMCut, hiMCut);
      break;
    }
    case 11: { // dielectron photoproduction
      int pdg = 11;
      elemProcess = new UpcTwoPhotonDilep(pdg);
      break;
    }
    case 13: { // dimuon photoproduction
      int pdg = 13;
      elemProcess = new UpcTwoPhotonDilep(pdg);
      break;
    }
    case 15: { // ditau photoproduction
      int pdg = 15;
      elemProcess = new UpcTwoPhotonDilep(pdg);
      break;
    }
    case 51: { // spin-0 ALP production
      elemProcess = new UpcTwoPhotonALP(alpMass, alpWidth); // todo: spin-dependent cross section?
      break;
    }
    case 111: { // pi0pi0 meson photoproduction
      elemProcess = new UpcTwoPhotonDipion(doMassCut, lowMCut, hiMCut);
      break;
    }
    default: {
      PLOG_FATAL << "Unknown process ID! Check manual and enter a correct ID! Exiting...";
      std::_Exit(-1);
    }
  }
}

void UpcCrossSection::init()
{
  PLOG_INFO << "Initializing caches ...";
  // update scaling factor
  factor = Z * Z * phys_consts::alpha / M_PI / M_PI / phys_consts::hc / phys_consts::hc;

  // calculate Woods-Saxon rho0 from R and a
  rho0 = calcWSRho();

  // prepare caches
  prepareGAA();     // G_AA and Fourier-transformed G_AA
  prepareFormFac(); // nuclear form factor
  if (breakupMode > 1) {
    prepareBreakupProb();
  }
  prepareTwoPhotonLumi(); // calculate two-photon luminosity and save into a file
}

template <typename ArrayType>
double UpcCrossSection::simpson(int n, ArrayType* v, double h)
{
  double sum = v[0] + v[n - 1];
  for (int i = 1; i < n - 1; i += 2) {
    sum += 4 * v[i];
  }
  for (int i = 2; i < n - 1; i += 2) {
    sum += 2 * v[i];
  }
  return sum * h / 3;
}

double UpcCrossSection::calcWSRho()
{
  double bmax = 20;
  double db = bmax / (nb - 1);
  for (int ib = 0; ib < nb; ib++) {
    double r = ib * db;
    vRho[ib] = r * r / (1 + exp((r - R) / a));
  }
  double wsRho0 = A / simpson(nb, vRho, db) / 4 / M_PI;
  return wsRho0;
}

double UpcCrossSection::fluxPoint(const double b, const double k)
{
  // flux divided by k
  double g = g1;
  double x = b * k / g / phys_consts::hc;
  double K0 = x > 1e-10 ? TMath::BesselK0(x) : 0;
  double K1 = x > 1e-10 ? TMath::BesselK1(x) : 0;
  double result = factor * k / g / g * (K1 * K1 + K0 * K0 / g / g);
  if (debug > 1) {
    PLOG_DEBUG << "result = " << result;
  }
  return result;
}

double UpcCrossSection::fluxFormInt(double* x, double* par)
{
  double k = x[0];
  double b = par[0];
  double w = par[1];
  double g = g1;

  if (debug > 1) {
    PLOG_DEBUG << "b = " << b << " w = " << w << " g = " << g;
  }

  double t = k * k + w * w / g / g;
  double ff = getCachedFormFac(t);
  double result = k * k * ff / t * TMath::BesselJ1(b * k / phys_consts::hc);
  if (debug > 1) {
    PLOG_DEBUG << "result = " << result;
  }
  return result;
}

double UpcCrossSection::fluxForm(const double b, const double k, TF1* fFluxFormInt)
{
  // flux divided by k
  if (isPoint) {
    return fluxPoint(b, k);
  }

  if (b > 2 * R) {
    return fluxPoint(b, k);
  }

  fFluxFormInt->SetParameter(0, b);
  fFluxFormInt->SetParameter(1, k);
  fFluxFormInt->SetParameter(2, g1);
  double Q = fFluxFormInt->Integral(0, 10, 1e-5) / A;
  double result = factor * Q * Q / k;

  if (debug > 1) {
    PLOG_DEBUG << "result = " << result;
  }

  return result;
}

double UpcCrossSection::calcTwoPhotonLumi(double M, double Y, TF1* fFluxForm, const TGraph* gGAA)
{
  // double differential luminosity
  double k1 = M / 2. * exp(Y);
  double k2 = M / 2. * exp(-Y);

  double b1min = isPoint ? 1 * R : 0.05 * R;
  double b2min = isPoint ? 1 * R : 0.05 * R;
  double b1max = TMath::Max(5. * g1 * phys_consts::hc / k1, 5 * R);
  double b2max = TMath::Max(5. * g2 * phys_consts::hc / k2, 5 * R);
  double log_delta_b1 = (log(b1max) - log(b1min)) / nb1;
  double log_delta_b2 = (log(b2max) - log(b2min)) / nb2;

  vector<double> flux(nb2);
  for (int j = 0; j < nb2; j++) {
    double b2l = b2min * exp(j * log_delta_b2);
    double b2h = b2min * exp((j + 1) * log_delta_b2);
    double b2 = (b2h + b2l) / 2.;
    flux[j] = fluxForm(b2, k2, fFluxForm);
  }

  double sum = 0;
  for (int i = 0; i < nb1; i++) {
    double sum_b2 = 0.;
    double b1l = b1min * exp(i * log_delta_b1);
    double b1h = b1min * exp((i + 1) * log_delta_b1);
    double b1 = (b1h + b1l) / 2.;
    for (int j = 0; j < nb2; j++) {
      double b2l = b2min * exp(j * log_delta_b2);
      double b2h = b2min * exp((j + 1) * log_delta_b2);
      double b2 = (b2h + b2l) / 2.;
      double sum_phi = 0.;
      for (int k = 0; k < ngi10; k++) {
        if (abscissas10[k] < 0) {
          continue;
        }
        double phi = M_PI * abscissas10[k];
        double b = sqrt(b1 * b1 + b2 * b2 + 2. * b1 * b2 * cos(phi));
        double breakup = breakupMode == 1 ? 1 : getCachedBreakupProb(b);
        double gaa = b < 20. ? gGAA->Eval(b) : 1;
        sum_phi += breakup * gaa * weights10[k];
      }
      sum_b2 += flux[j] * sum_phi * b2 * (b2h - b2l);
    }
    sum += fluxForm(b1, k1, fFluxForm) * sum_b2 * b1 * (b1h - b1l);
  }
  double lumi = 2 * M_PI * M_PI * M * sum;
  return lumi;
}

void UpcCrossSection::calcTwoPhotonLumiPol(double& ns, double& np, double M, double Y, TF1* fFluxForm, const TGraph* gGAA)
{
  double k1 = M / 2. * exp(Y);
  double k2 = M / 2. * exp(-Y);

  double b1min = isPoint ? 1 * R : 0.05 * R;
  double b2min = isPoint ? 1 * R : 0.05 * R;
  double b1max = TMath::Max(5. * g1 * phys_consts::hc / k1, 5 * R);
  double b2max = TMath::Max(5. * g2 * phys_consts::hc / k2, 5 * R);
  double log_delta_b1 = (log(b1max) - log(b1min)) / nb1;
  double log_delta_b2 = (log(b2max) - log(b2min)) / nb2;

  double flux[nb2];
  for (int j = 0; j < nb2; ++j) {
    double b2l = b2min * exp(j * log_delta_b2);
    double b2h = b2min * exp((j + 1) * log_delta_b2);
    double b2 = (b2h + b2l) / 2.;
    flux[j] = fluxForm(b2, k2, fFluxForm);
  }

  // calculating two integrals simultaneously
  double sum_b1_s = 0; // scalar part
  double sum_b1_p = 0; // pseudoscalar part
  for (int i = 0; i < nb1; i++) {
    double b1l = b1min * exp(i * log_delta_b1);
    double b1h = b1min * exp((i + 1) * log_delta_b1);
    double b1 = (b1h + b1l) / 2.;
    double ff_b1 = fluxForm(b1, k1, fFluxForm);
    double sum_b2_s = 0;
    double sum_b2_p = 0;
    for (int j = 0; j < nb2; ++j) {
      double b2l = b2min * exp(j * log_delta_b2);
      double b2h = b2min * exp((j + 1) * log_delta_b2);
      double b2 = (b2h + b2l) / 2.;
      double sum_phi_s = 0.;
      double sum_phi_p = 0.;
      for (int k = 0; k < ngi10; k++) {
        if (abscissas10[k] < 0) {
          continue;
        }
        double phi = M_PI * abscissas10[k];
        double cphi = TMath::Cos(phi);
        double sphi = TMath::Sin(phi);
        double b = TMath::Sqrt(b1 * b1 + b2 * b2 - 2. * b1 * b2 * cphi);
        double breakup = breakupMode == 1 ? 1 : getCachedBreakupProb(b);
        double gaa = b < 20 ? gGAA->Eval(b) : 1;
        sum_phi_s += breakup * gaa * weights10[k] * cphi * cphi;
        sum_phi_p += breakup * gaa * weights10[k] * sphi * sphi;
      }
      double ff_b2 = flux[j];
      sum_b2_s += sum_phi_s * ff_b2 * b2 * (b2h - b2l);
      sum_b2_p += sum_phi_p * ff_b2 * b2 * (b2h - b2l);
    }
    sum_b1_s += sum_b2_s * ff_b1 * b1 * (b1h - b1l);
    sum_b1_p += sum_b2_p * ff_b1 * b1 * (b1h - b1l);
  }
  ns = 2 * M_PI * M_PI * M * sum_b1_s;
  np = 2 * M_PI * M_PI * M * sum_b1_p;
}

void UpcCrossSection::fillCrossSectionZM(TH2D* hCrossSectionZM,
                                         double zmin, double zmax, int nz,
                                         double mmin, double mmax, int nm,
                                         int flag)
{
  double m, z;
  double dm = (mmax - mmin) / nm;
  double dz = (zmax - zmin) / nz;
  double cs;
  for (int im = 1; im <= nm; im++) {
    m = mmin + dm * (im - 1);
    for (int iz = 1; iz <= nz; iz++) {
      z = zmin + dz * (iz - 1);
      if (flag == 0) { // the usual unpolarized cross section
        cs = elemProcess->calcCrossSectionZM(z, m);
      }
      if (flag == 1) { // scalar part
        cs = elemProcess->calcCrossSectionZMPolS(z, m);
      }
      if (flag == 2) { // pseudoscalar part
        cs = elemProcess->calcCrossSectionZMPolPS(z, m);
      }
      hCrossSectionZM->SetBinContent(iz, im, cs);
    }
  }
  double scalingFactor = phys_consts::hc * phys_consts::hc * 1e7; // to [nb]
  hCrossSectionZM->Scale(scalingFactor / dm);
}

void UpcCrossSection::prepareGAA()
{
  double bmax = 20;
  double db = bmax / (nb - 1);
  double ssm = pow(sqrts, 2) / pow(2 * phys_consts::mProt + 2.1206, 2);
  double csNN = 0.1 * (34.41 + 0.2720 * pow(log(ssm), 2) + 13.07 * pow(ssm, -0.4473) - 7.394 * pow(ssm, -0.5486)); // PDG 2016
  // calculate rho and TA
  double TAb[nb];
  for (int ib = 0; ib < nb; ib++) {
    double b = ib * db;
    for (int iz = 0; iz < nb; iz++) {
      double z = iz * db;
      double r = TMath::Sqrt(b * b + z * z);
      rho[ib][iz] = rho0 / (1 + exp((r - R) / a));
    }
    TA[ib] = 2 * simpson(nb, rho[ib], db);
    TAb[ib] = TA[ib] * b;
    vb[ib] = b;
  }
  auto* gTA = new TGraph(nb, vb, TA);

  // calculate G_AA
  for (int ib = 0; ib < nb; ib++) {
    double b = ib * db;
    for (int is = 0; is < nb; is++) {
      double s = is * db;
      double sum_phi = 0;
      for (int k = 0; k < ngi10; k++) {
        if (abscissas10[k] < 0)
          continue;
        double r = TMath::Sqrt(b * b + s * s + 2 * b * s * TMath::Cos(M_PI * abscissas10[k]));
        sum_phi += 2 * M_PI * weights10[k] * gTA->Eval(r);
      }
      vs[is] = 2 * s * gTA->Eval(s) * sum_phi;
    }
    vGAA[ib] = exp(-csNN * simpson(nb, vs, db));
  }

  delete gTA;
}

void UpcCrossSection::prepareBreakupProb()
{
  constexpr double bmin = 1e-6;
  constexpr double bmax = 1000;
  constexpr int nbc = 1000000;
  constexpr double db = (bmax - bmin) / nbc;
  for (int i = 0; i < nbc; i++) {
    double b = bmin + db * i;
    double prob = calcBreakupProb(b, breakupMode);
    vCachedBreakup[i] = prob;
  }
}

double UpcCrossSection::getCachedBreakupProb(double b)
{
  constexpr double bmin = 1e-6;
  constexpr double bmax = 1000;
  constexpr int nbc = 1000000;
  constexpr double db = (bmax - bmin) / nbc;
  if (breakupMode == 3 && b > bmax) {
    return 1.;
  }
  if ((breakupMode == 2 || breakupMode == 4) && b > bmax) {
    return 0.;
  }
  double frac = (b - bmin) / db;
  int idx1 = floor(frac);
  double prob1 = vCachedBreakup[idx1];
  double prob2 = vCachedBreakup[idx1 + 1];
  double prob = prob1 + (prob2 - prob1) * (frac - idx1);
  return prob;
}

double UpcCrossSection::calcFormFac(double Q2)
{
  double Q = sqrt(Q2) / phys_consts::hc;
  double coshVal = TMath::CosH(M_PI * Q * a);
  double sinhVal = TMath::SinH(M_PI * Q * a);
  double ff = 4 * M_PI * M_PI * rho0 * a * a * a / (Q * a * Q * a * sinhVal * sinhVal) *
              (M_PI * Q * a * coshVal * TMath::Sin(Q * R) - Q * R * TMath::Cos(Q * R) * sinhVal);
  ff += 8 * M_PI * rho0 * a * a * a * exp(-R / a) / (1 + Q * Q * a * a) / (1 + Q * Q * a * a);
  return ff;
}

void UpcCrossSection::prepareFormFac()
{
  for (int iQ2 = 0; iQ2 < nQ2; iQ2++) {
    double Q2 = Q2min + iQ2 * dQ2;
    double ff = calcFormFac(Q2);
    vCachedFormFac[iQ2] = ff;
  }
}

double UpcCrossSection::getCachedFormFac(double Q2)
{
  if (Q2 > Q2max) {
    return 0;
  }
  double frac = (Q2 - Q2min) / dQ2;
  int idx1 = floor(frac);
  double ff1 = vCachedFormFac[idx1];
  double ff2 = vCachedFormFac[idx1 + 1];
  double ff = ff1 + (ff2 - ff1) * (frac - idx1);
  return ff;
}

void UpcCrossSection::prepareTwoPhotonLumi()
{
  // wait here if an other generator is working on this
  TString fnlock{lumiFileDirectory+"/.lumiIsCalculated"};
  if (usePolarizedCS) {
    fnlock += "Pol";
  }
  bool itIsMyTurn{false};
  while (!itIsMyTurn) {
    auto fd = open(fnlock.Data(), O_CREAT | O_EXCL, 0644);
    if (fd > 0) {
      itIsMyTurn = true;
    } else {
      sleep(1);
    }
  }
  
  bool isFound = false;
  TString fname{lumiFileDirectory+"/twoPhotonLumi"};
  fname += usePolarizedCS ? "Pol.root" : ".root";
  
  if (!gSystem->AccessPathName(fname.Data()) && !usePolarizedCS) {
    PLOG_INFO << "Found pre-calculated unpolarized 2D luminosity";
    isFound = true;
  }
  if (!gSystem->AccessPathName(fname.Data()) && usePolarizedCS) {
    PLOG_INFO << "Found pre-calculated polarized 2D luminosity";
    isFound = true;
  }
  if (!isFound) {
    PLOG_INFO << "Precalculated 2D luminosity is not found. Starting all over...";
    double dy = (ymax - ymin) / ny;
    double dm = (mmax - mmin) / nm;
    auto* f2DLumi = new TFile(fname, "recreate");
    // histograms for unpolarized case
    TH2D* hD2LDMDY = nullptr;
    // histograms for polarized case
    TH2D* hD2LDMDY_s = nullptr;
    TH2D* hD2LDMDY_p = nullptr;
    if (usePolarizedCS) {
      hD2LDMDY_s = new TH2D("hD2LDMDY_s", ";;", nm, mmin, mmax, ny, ymin, ymax);
      hD2LDMDY_p = new TH2D("hD2LDMDY_p", ";;", nm, mmin, mmax, ny, ymin, ymax);
    } else {
      hD2LDMDY = new TH2D("hD2LDMDY", ";;", nm, mmin, mmax, ny, ymin, ymax);
    }
    ROOT::EnableThreadSafety();
    int im, iy;
    int progress = 0;
    int total = ny * nm;
    // using ether parallel or serial implementation
#ifdef USE_OPENMP
    omp_set_num_threads(numThreads);
#pragma omp parallel default(none)                                   \
  shared(hD2LDMDY, hD2LDMDY_s, hD2LDMDY_p, progress) private(im, iy) \
    firstprivate(nb, vb, total, numThreads, dm, dy, ymin, ymax, mmin, mmax, nm, ny, abscissas10, weights10, vGAA)
    {
      auto* fFluxFormInt = new TF1(Form("fFluxFormInt_private_%d", omp_get_thread_num()), fluxFormInt, 0, 10, 3);
      auto* gGAA = new TGraph(nb, vb, vGAA);
      TH2D* hD2LDMDY_private = nullptr;
      TH2D* hD2LDMDY_private_s = nullptr;
      TH2D* hD2LDMDY_private_p = nullptr;
      if (usePolarizedCS) {
        hD2LDMDY_private_s = new TH2D(Form("hD2LDMDY_private_s_%d", omp_get_thread_num()), ";;", nm, mmin, mmax, ny, ymin, ymax);
        hD2LDMDY_private_p = new TH2D(Form("hD2LDMDY_private_p_%d", omp_get_thread_num()), ";;", nm, mmin, mmax, ny, ymin, ymax);
      } else {
        hD2LDMDY_private = new TH2D(Form("hD2LDMDY_private_%d", omp_get_thread_num()), ";;", nm, mmin, mmax, ny, ymin, ymax);
      }
      int threadNum = omp_get_thread_num();
      int lowM = nm * threadNum / numThreads;
      int highM = nm * (threadNum + 1) / numThreads;
      for (im = lowM; im < highM; im++) {
        double m = mmin + dm * im;
        for (iy = 0; iy < ny; iy++) {
          double y = ymin + dy * iy;
          if (usePolarizedCS) {
            double lumi_s, lumi_p;
            calcTwoPhotonLumiPol(lumi_s, lumi_p, m, y, fFluxFormInt, gGAA);
            hD2LDMDY_private_s->SetBinContent(im, iy, lumi_s * dm * dy);
            hD2LDMDY_private_p->SetBinContent(im, iy, lumi_p * dm * dy);
          } else {
            double lumi = calcTwoPhotonLumi(m, y, fFluxFormInt, gGAA);
            hD2LDMDY_private->SetBinContent(im, iy, lumi * dm * dy);
          }
          progress++;
        }
        if (threadNum == 0) {
          double progressBar = 100. * progress / total;
          PLOG_INFO << "Calculating two-photon luminosity: " << fixed << setprecision(2) << progressBar << "%";
        }
      }
#pragma omp critical
      {
        if (usePolarizedCS) {
          hD2LDMDY_s->Add(hD2LDMDY_private_s);
          hD2LDMDY_p->Add(hD2LDMDY_private_p);
        } else {
          hD2LDMDY->Add(hD2LDMDY_private);
        }
        delete hD2LDMDY_private;
        delete hD2LDMDY_private_s;
        delete hD2LDMDY_private_p;
      }
    }
    omp_set_num_threads(1);
#else
    auto* fFluxFormInt = new TF1(Form("fFluxFormInt", omp_get_thread_num()), fluxFormInt, 0, 10, 3);
    auto* gGAA = new TGraph(nb, vb, vGAA);
    for (im = 0; im < nm; im++) {
      double M = mmin + dm * im;
      for (iy = 0; iy < ny; iy++) {
        double y = ymin + dy * iy;
        if (usePolarizedCS) {
          double lumi_s, lumi_p;
          calcTwoPhotonLumiPol(lumi_s, lumi_p, m, y, fFluxFormInt, gGAA);
          hD2LDMDY_s->SetBinContent(im, iy, lumi_s * dm * dy);
          hD2LDMDY_p->SetBinContent(im, iy, lumi_p * dm * dy);
        } else {
          double lumi = calcTwoPhotonLumi(m, y, fFluxFormInt, gGAA);
          hD2LDMDY->SetBinContent(im, iy, lumi * dm * dy);
        }
      }
      double progressBar = 100. * progress / total;
      PLOG_INFO << "Calculating two-photon luminosity: " << fixed << setprecision(2) << progressBar << "%";
    }
#endif
    if (usePolarizedCS) {
      hD2LDMDY_s->Write();
      hD2LDMDY_p->Write();
    } else {
      hD2LDMDY->Write();
    }
    f2DLumi->Close();
    PLOG_INFO << "Two-photon luminosity was written to " << fname;
  }
  
  // remove lock file
  if (remove(fnlock.Data()) != 0) {
    PLOG_WARNING << "Lock file " << fnlock << " was not properly removed!";
  }  
}

void UpcCrossSection::calcNucCrossSectionYM(TH2D* hCrossSectionYM, vector<vector<double>>& hPolCSRatio, double& totCS)
{
  PLOG_INFO << "Calculating nuclear cross section...";

  double dy = (ymax - ymin) / ny;
  double dm = (mmax - mmin) / nm;

  TH2D* hD2LDMDY = nullptr;
  TH2D* hD2LDMDY_s = nullptr;
  TH2D* hD2LDMDY_p = nullptr;

  TString fname{lumiFileDirectory+"/twoPhotonLumi"};
  fname += usePolarizedCS ? "Pol.root" : ".root";
  auto* f2DLumi = new TFile(fname, "r");

  if (usePolarizedCS) {
    hD2LDMDY_s = (TH2D*)f2DLumi->Get("hD2LDMDY_s");
    hD2LDMDY_p = (TH2D*)f2DLumi->Get("hD2LDMDY_p");
  } else { // loading pre-cached two-photon luminosity
    hD2LDMDY = (TH2D*)f2DLumi->Get("hD2LDMDY");
  }

  // calculating nuclear cross section
  ROOT::EnableThreadSafety();
  int im, iy, ib;
  int progress = 0;
  int total = nm * ny;
  std::vector<std::vector<double>> cs(nm, std::vector<double>(ny, 0));
  std::vector<std::vector<double>> cs_rat(nm, std::vector<double>(ny, 0));
  omp_set_num_threads(numThreads);
#pragma omp parallel default(none)                                                   \
  shared(cs, cs_rat, hD2LDMDY, hD2LDMDY_s, hD2LDMDY_p, progress) private(im, iy, ib) \
    firstprivate(elemProcess, nb, vb, total, numThreads, dm, dy, ymin, ymax, mmin, mmax, nm, ny, abscissas10, weights10, vGAA)
  {
    vector<vector<double>> cs_private(nm, vector<double>(ny, 0));
    vector<vector<double>> rat_private(nm, vector<double>(ny, 0));
    TH2D* hD2LDMDY_private = nullptr;
    TH2D* hD2LDMDY_private_s = nullptr;
    TH2D* hD2LDMDY_private_p = nullptr;
    if (usePolarizedCS) {
      hD2LDMDY_private_s = (TH2D*)hD2LDMDY_s->Clone(Form("hD2LDMDY_private_s_%d", omp_get_thread_num()));
      hD2LDMDY_private_p = (TH2D*)hD2LDMDY_p->Clone(Form("hD2LDMDY_private_p_%d", omp_get_thread_num()));
    } else {
      hD2LDMDY_private = (TH2D*)hD2LDMDY->Clone(Form("hD2LDMDY_private_%d", omp_get_thread_num()));
    }
    int threadNum = omp_get_thread_num();
    int lowM = nm * threadNum / numThreads;
    int highM = nm * (threadNum + 1) / numThreads;
    for (im = lowM; im < highM; im++) {
      double m = mmin + dm * im;
      for (iy = 0; iy < ny; iy++) {
        if (!usePolarizedCS) { // unpolarized cross section
          double lumi = hD2LDMDY_private->GetBinContent(im, iy);
          cs_private[im][iy] = elemProcess->calcCrossSectionM(m) * lumi;
        } else { // polarized
          double cs_s = elemProcess->calcCrossSectionMPolS(m);
          double cs_p = elemProcess->calcCrossSectionMPolPS(m);
          double lumi_s = hD2LDMDY_private_s->GetBinContent(im, iy);
          double lumi_p = hD2LDMDY_private_p->GetBinContent(im, iy);
          double nuccs_s = lumi_s * cs_s; // scalar part
          double nuccs_p = lumi_p * cs_p; // psudoscalar part
          double nuccs = nuccs_s + nuccs_p;
          cs_private[im][iy] = nuccs * 1e7; // fm^2 -> nb
          rat_private[im][iy] = nuccs_s / nuccs_p;
        }
        progress++;
      }
      if (threadNum == 0) {
        double progressBar = 100. * progress / total;
        PLOG_INFO << "Calculating nuclear cross section... " << fixed << setprecision(2) << progressBar << "%";
      }
    }
#pragma omp critical
    {
      for (im = lowM; im < highM; im++) {
        for (iy = 0; iy < ny; iy++) {
          cs[im][iy] = cs_private[im][iy];
          cs_rat[im][iy] = rat_private[im][iy];
        }
      }
    }
    delete hD2LDMDY_private;
    delete hD2LDMDY_private_s;
    delete hD2LDMDY_private_p;
  }
  omp_set_num_threads(1);

  delete hD2LDMDY;
  delete hD2LDMDY_s;
  delete hD2LDMDY_p;
  f2DLumi->Close();
  delete f2DLumi;

  PLOG_INFO << "Calculating nuclear cross section...Done!";

  // filling histograms
  double cssum = 0;
  for (int i = 0; i < nm; i++) {
    for (int j = 0; j < ny; j++) {
      double cs_ij = cs[i][j];
      hCrossSectionYM->SetBinContent(j + 1, i + 1, cs_ij);
      if (usePolarizedCS) {
        double rat_ij = cs_rat[i][j];
        hPolCSRatio[j][i] = rat_ij;
      }
    }
  }

  totCS = hCrossSectionYM->Integral() * 1e-6;
  PLOG_INFO << "Total nuclear cross section = " << fixed << setprecision(6) << totCS << " mb";
}

// Function from Starlight
// (by S.R.Klein, J.Nystrand, J.Seger, Y.Gorbunov, J.Butterworth)
double UpcCrossSection::calcBreakupProb(const double impactparameter, const int mode)
{
  static double ee[10001], eee[162], se[10001];

  double _pPhotonBreakup = 0.; // Might default the probability with a different value?
  double b = impactparameter;
  int zp = 82; // What about _beam2? Generic approach?
  int ap = 208;

  // Was initialized at the start of the function originally, been moved inward.
  double pxn = 0.;
  double p1n = 0.;

  double _beamLorentzGamma = UpcCrossSection::g1;
  double hbarcmev = 197.3269718;
  double pi = 3.14159;
  // Used to be done prior to entering the function. Done properly for assymetric?
  double gammatarg = 2. * _beamLorentzGamma * _beamLorentzGamma - 1.;
  double omaxx = 0.;
  // This was done prior entering the function as well
  if (_beamLorentzGamma > 500.) {
    omaxx = 1.E10;
  } else {
    omaxx = 1.E7;
  }

  double e1[23] = {0., 103., 106., 112., 119., 127., 132., 145., 171., 199., 230., 235.,
                   254., 280., 300., 320., 330., 333., 373., 390., 420., 426., 440.};
  double s1[23] = {0., 12.0, 11.5, 12.0, 12.0, 12.0, 15.0, 17.0, 28.0, 33.0,
                   52.0, 60.0, 70.0, 76.0, 85.0, 86.0, 89.0, 89.0, 75.0, 76.0, 69.0, 59.0, 61.0};
  double e2[12] = {0., 2000., 3270., 4100., 4810., 6210., 6600.,
                   7790., 8400., 9510., 13600., 16400.};
  double s2[12] = {0., .1266, .1080, .0805, .1017, .0942, .0844, .0841, .0755, .0827,
                   .0626, .0740};
  double e3[29] = {0., 26., 28., 30., 32., 34., 36., 38., 40., 44., 46., 48., 50., 52., 55.,
                   57., 62., 64., 66., 69., 72., 74., 76., 79., 82., 86., 92., 98., 103.};
  double s3[29] = {0., 30., 21.5, 22.5, 18.5, 17.5, 15., 14.5, 19., 17.5, 16., 14.,
                   20., 16.5, 17.5, 17., 15.5, 18., 15.5, 15.5, 15., 13.5, 18., 14.5, 15.5, 12.5, 13.,
                   13., 12.};
  static double sa[161] = {0., 0., .004, .008, .013, .017, .021, .025, .029, .034, .038, .042, .046,
                           .051, .055, .059, .063, .067, .072, .076, .08, .085, .09, .095, .1, .108, .116,
                           .124, .132, .14, .152, .164, .176, .188, .2, .22, .24, .26, .28, .3, .32, .34,
                           .36, .38, .4, .417, .433, .450, .467, .483, .5, .51, .516, .52, .523, .5245,
                           .525, .5242,
                           .5214, .518, .512, .505, .495, .482, .469, .456, .442, .428, .414, .4, .386,
                           .370, .355, .34, .325, .310, .295, .280, .265, .25, .236, .222, .208, .194,
                           .180, .166,
                           .152, .138, .124, .11, .101, .095, .09, .085, .08, .076, .072, .069, .066,
                           .063, .06, .0575, .055, .0525, .05, .04875, .0475, .04625, .045, .04375,
                           .0425, .04125, .04, .03875, .0375, .03625, .035, .03375, .0325, .03125, .03,
                           .02925, .0285, .02775, .027, .02625, .0255, .02475, .024, .02325, .0225,
                           .02175, .021, .02025, .0195, .01875, .018, .01725, .0165, .01575, .015,
                           .01425, .0135, .01275, .012, .01125, .0105, .00975, .009, .00825, .0075,
                           .00675, .006, .00525, .0045, .00375, .003, .00225, .0015, .00075, 0.};

  double sen[161] = {0., 0., .012, .025, .038, .028, .028, .038, .035, .029, .039, .035,
                     .038, .032, .038, .041, .041, .049, .055, .061, .072, .076, .070, .067,
                     .080, .103, .125, .138, .118, .103, .129, .155, .170, .180, .190, .200,
                     .215, .250, .302, .310, .301, .315, .330, .355, .380, .400, .410, .420,
                     .438, .456, .474, .492, .510, .533, .556, .578, .6, .62, .63, .638,
                     .640, .640, .637, .631, .625, .618, .610, .600, .580, .555, .530, .505,
                     .480, .455, .435, .410, .385, .360, .340, .320, .300, .285, .270, .255,
                     .240, .225, .210, .180, .165, .150, .140, .132, .124, .116, .108, .100,
                     .092, .084, .077, .071, .066, .060, .055, .051, .048, .046, .044, .042,
                     .040, .038, .036, .034, .032, .030, .028, .027, .026, .025, .025, .025,
                     .024, .024, .024, .024, .024, .023, .023, .023, .023, .023, .022, .022,
                     .022, .022, .022, .021, .021, .021, .020, .020,
                     .020, .019, .018, .017, .016, .015, .014, .013, .012, .011, .010, .009,
                     .008, .007, .006, .005, .004, .003, .002, .001, 0.};

  // gammay,p gamma,n of Armstrong begin at 265 incr 25

  double sigt[160] = {0., .4245, .4870, .5269, .4778, .4066, .3341, .2444, .2245, .2005,
                      .1783, .1769, .1869, .1940, .2117, .2226, .2327, .2395, .2646, .2790, .2756,
                      .2607, .2447, .2211, .2063, .2137, .2088, .2017, .2050, .2015, .2121, .2175,
                      .2152, .1917, .1911, .1747, .1650, .1587, .1622, .1496, .1486, .1438, .1556,
                      .1468, .1536, .1544, .1536, .1468, .1535, .1442, .1515, .1559, .1541, .1461,
                      .1388, .1565, .1502, .1503, .1454, .1389, .1445, .1425, .1415, .1424, .1432,
                      .1486, .1539, .1354, .1480, .1443, .1435, .1491, .1435, .1380, .1317, .1445,
                      .1375, .1449, .1359, .1383, .1390, .1361, .1286, .1359, .1395, .1327, .1387,
                      .1431, .1403, .1404, .1389, .1410, .1304, .1363, .1241, .1284, .1299, .1325,
                      .1343, .1387, .1328, .1444, .1334, .1362, .1302, .1338, .1339, .1304, .1314,
                      .1287, .1404, .1383, .1292, .1436, .1280, .1326, .1321, .1268, .1278, .1243,
                      .1239, .1271, .1213, .1338, .1287, .1343, .1231, .1317, .1214, .1370, .1232,
                      .1301, .1348, .1294, .1278, .1227, .1218, .1198, .1193, .1342, .1323, .1248,
                      .1220, .1139, .1271, .1224, .1347, .1249, .1163, .1362, .1236, .1462, .1356,
                      .1198, .1419, .1324, .1288, .1336, .1335, .1266};

  double sigtn[160] = {0., .3125, .3930, .4401, .4582, .3774, .3329, .2996, .2715, .2165,
                       .2297, .1861, .1551, .2020, .2073, .2064, .2193, .2275, .2384, .2150, .2494,
                       .2133, .2023, .1969, .1797, .1693, .1642, .1463, .1280, .1555, .1489, .1435,
                       .1398, .1573, .1479, .1493, .1417, .1403, .1258, .1354, .1394, .1420, .1364,
                       .1325, .1455, .1326, .1397, .1286, .1260, .1314, .1378, .1353, .1264, .1471,
                       .1650, .1311, .1261, .1348, .1277, .1518, .1297, .1452, .1453, .1598, .1323,
                       .1234, .1212, .1333, .1434, .1380, .1330, .12, .12, .12, .12, .12, .12, .12, .12,
                       .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12,
                       .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12,
                       .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12,
                       .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12,
                       .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12, .12};

  static int IFIRSTP = 0;

  double si1 = 0, g1 = 0, o1 = 0;
  int ne = 0, ij = 0;
  double delo = 0, omax = 0, gk1m = 0;
  static double scon = 0., zcon = 0., o0 = 0.;

  double x = 0, y = 0, eps = 0, eta = 0, em = 0, exx = 0, s = 0, ictr = 0, pom = 0, vec = 0, gk1 = 0;

  //  maximum energy for GDR dissocation (in target frame, in MeV)

  double omax1n = 24.01;

  if (IFIRSTP != 0)
    goto L100;

  IFIRSTP = 1;

  // This is dependenant on gold or lead....Might need to expand
  if (zp == 79) {

    ap = 197;
    si1 = 540.;
    g1 = 4.75;

    // peak and minimum energies for GDR excitation (in MeV)
    o1 = 13.70;
    o0 = 8.1;
  } else {
    zp = 82; // assumed to be lead
    ap = 208;
    si1 = 640.;
    g1 = 4.05;
    o1 = 13.42;
    o0 = 7.4;
    for (int j = 1; j <= 160; j++) {

      sa[j] = sen[j];
    }
  }
  // Part II of initialization
  delo = .05;
  //.1 to turn mb into fm^2
  scon = .1 * g1 * g1 * si1;
  zcon = zp / (gammatarg * (pi) * (hbarcmev)) * zp / (gammatarg * (pi) * (hbarcmev)) / 137.04; // alpha?

  // single neutron from GDR, Veyssiere et al. Nucl. Phys. A159, 561 (1970)
  for (int i = 1; i <= 160; i++) {
    eee[i] = o0 + .1 * (i - 1);
    sa[i] = 100. * sa[i];
  }
  // See Baltz, Rhoades-Brown, and Weneser, Phys. Rev. E 54, 4233 (1996)
  // for details of the following photo cross-sections
  eee[161] = 24.1;
  ne = int((25. - o0) / delo) + 1;
  // GDR any number of neutrons, Veyssiere et al., Nucl. Phys. A159, 561 (1970)
  for (int i = 1; i <= ne; i++) {
    ee[i] = o0 + (i - 1) * delo;
    // cout<<" ee 1 "<<ee[i]<<"  "<<i<<endl;

    se[i] = scon * ee[i] * ee[i] / (((o1 * o1 - ee[i] * ee[i]) * (o1 * o1 - ee[i] * ee[i])) + ee[i] * ee[i] * g1 * g1);
  }
  ij = ne; // Risky?
  // 25-103 MeV, Lepretre, et al., Nucl. Phys. A367, 237 (1981)
  for (int j = 1; j <= 27; j++) {
    ij = ij + 1;
    ee[ij] = e3[j];
    // cout<<" ee 2 "<<ee[ij]<<"  "<<ij<<endl;

    se[ij] = .1 * ap * s3[j] / 208.;
  }
  // 103-440 MeV, Carlos, et al., Nucl. Phys. A431, 573 (1984)
  for (int j = 1; j <= 22; j++) {
    ij = ij + 1;
    ee[ij] = e1[j];
    // cout<<" ee 3 "<<ee[ij]<<"  "<<ij<<endl;
    se[ij] = .1 * ap * s1[j] / 208.;
  }
  // 440 MeV-2 GeV Armstrong et al.
  for (int j = 9; j <= 70; j++) {
    ij = ij + 1;
    ee[ij] = ee[ij - 1] + 25.;
    // cout<<" ee 4 "<<ee[ij]<<"  "<<ij<<endl;
    se[ij] = .1 * (zp * sigt[j] + (ap - zp) * sigtn[j]);
  }
  // 2-16.4 GeV Michalowski; Caldwell
  for (int j = 1; j <= 11; j++) {
    ij = ij + 1;
    ee[ij] = e2[j];
    // cout<<" ee 5 "<<ee[ij]<<"   "<<ij<<endl;
    se[ij] = .1 * ap * s2[j];
  }
  // done with initaliation
  // Regge paramteres
  x = .0677;
  y = .129;
  eps = .0808;
  eta = .4525;
  em = .94;
  exx = pow(10, .05);

  // Regge model for high energy
  s = .002 * em * ee[ij];
  // make sure we reach LHC energies
  ictr = 100;
  if (gammatarg > (2. * 150. * 150.))
    ictr = 150;
  for (int j = 1; j <= ictr; j++) {
    ij = ij + 1;
    s = s * exx;
    ee[ij] = 1000. * .5 * (s - em * em) / em;
    // cout<<" ee 6 "<<ee[ij]<<"   "<<ij<<endl;
    pom = x * pow(s, eps);
    vec = y * pow(s, (-eta));
    se[ij] = .1 * .65 * ap * (pom + vec);
  }
  ee[ij + 1] = 99999999999.;

  // clear counters for 1N, XN
L100:

  p1n = 0.;
  pxn = 0.;
  // start XN calculation
  // what's the b-dependent highest energy of interest?

  omax = std::min(omaxx, 4. * gammatarg * (hbarcmev) / b);
  if (omax < o0)
    return _pPhotonBreakup;
  gk1m = TMath::BesselK1(ee[1] * b / ((hbarcmev)*gammatarg));
  int k = 2;
L212:
  if (ee[k] < omax) {
    gk1 = TMath::BesselK1(ee[k] * b / ((hbarcmev)*gammatarg));
    // Eq. 3 of BCW--NIM in Physics Research A 417 (1998) pp1-8:
    pxn = pxn + zcon * (ee[k] - ee[k - 1]) * .5 * (se[k - 1] * ee[k - 1] * gk1m * gk1m + se[k] * ee[k] * gk1 * gk1);
    k = k + 1;
    gk1m = gk1;
    goto L212;
  }
  // one neutron dissociation
  omax = std::min(omax1n, 4. * gammatarg * (hbarcmev) / b);
  gk1m = TMath::BesselK1(eee[1] * b / ((hbarcmev)*gammatarg));
  k = 2;
L102:
  if (eee[k] < omax) {
    gk1 = TMath::BesselK1(eee[k] * b / ((hbarcmev)*gammatarg));
    // Like Eq3 but with only the one neutron out GDR photo cross section input
    p1n = p1n + zcon * (eee[k] - eee[k - 1]) * .5 * (sa[k - 1] * eee[k - 1] * gk1m * gk1m + sa[k] * eee[k] * gk1 * gk1);
    k = k + 1;
    gk1m = gk1;
    goto L102;
  }

  if ((mode) == 1)
    _pPhotonBreakup = 1.;
  if ((mode) == 2)
    _pPhotonBreakup = (1 - exp(-1 * pxn)) * (1 - exp(-1 * pxn));
  if ((mode) == 3)
    _pPhotonBreakup = exp(-2 * pxn);
  if ((mode) == 4)
    _pPhotonBreakup = 2. * exp(-pxn) * (1. - exp(-pxn));

  // cout<<pxn<<" "<<zcon<<" "<<ee[k]<<" "<<se[k-1]<<" "<<gk1m<<"  "<<gk1<<"  "<<k<<"  "<<ee[k+1]<< "  "<<b<< endl;

  return _pPhotonBreakup;
}

double UpcCrossSection::getPhotonPt(double ePhot)
{
  constexpr double pi2x4 = 4 * M_PI * M_PI;
  double ereds = (ePhot * ePhot) / (gtot * gtot);
  int ePhot_key = ePhot * 1e3;
  constexpr int nbins = 5000;
  double pt = 0;
  auto it = photPtDistrMap.find(ePhot_key);
  if (it == photPtDistrMap.end()) {
    TH1D ptDistr(Form("ptDistr%d", ePhot_key), "", nbins, 0., 6. * phys_consts::hc / R);
    ptDistr.SetDirectory(nullptr);
    for (int bin = 1; bin <= nbins; bin++) {
      double pt = 6. * phys_consts::hc / R / nbins * bin;
      double arg = pt * pt + ereds;
      double sFFactPt1 = getCachedFormFac(arg);
      double prob = (sFFactPt1 * sFFactPt1) * pt * pt * pt / (pi2x4 * arg * arg);
      ptDistr.SetBinContent(bin, prob);
    }
    pt = ptDistr.GetRandom();
    photPtDistrMap.emplace(ePhot_key, std::make_pair(1, ptDistr));
  } else {
    pt = it->second.second.GetRandom();
    it->second.first++; // increment frequency
  }
  // clear map if it became too large
  int nElem = photPtDistrMap.size();
  if (nElem > 30000) {
    photPtDistrMap.clear();
  }
  return pt;
}

void UpcCrossSection::getPairMomentum(double mPair, double yPair, TLorentzVector& pPair)
{
  if (!useNonzeroGamPt) {
    double mtPair = mPair; // pairPt = 0
    pPair.SetPxPyPzE(0., 0., mtPair * TMath::SinH(yPair), mtPair * TMath::CosH(yPair));
  }
  if (useNonzeroGamPt) {
    double k1 = mPair / 2 * exp(yPair);
    double k2 = mPair / 2 * exp(-yPair);
    double angle1 = gRandom->Uniform(0, 2 * M_PI);
    double angle2 = gRandom->Uniform(0, 2 * M_PI);
    double pt1 = getPhotonPt(k1);
    double pt2 = getPhotonPt(k2);
    double px = pt1 * TMath::Cos(angle1) + pt2 * TMath::Cos(angle2);
    double py = pt1 * TMath::Sin(angle1) + pt2 * TMath::Sin(angle2);
    double pt = TMath::Sqrt(px * px + py * py);
    double mtPair = TMath::Sqrt(mPair * mPair + pt * pt);
    double pz = mtPair * TMath::SinH(yPair);
    double e = mtPair * TMath::CosH(yPair);
    pPair.SetPxPyPzE(px, py, pz, e);
  }
}
