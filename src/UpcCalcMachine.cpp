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

#include "UpcCalcMachine.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace std;

// out-of-line initialization for static members
double UpcCalcMachine::rho0 = 0;
double UpcCalcMachine::R = 6.68;
double UpcCalcMachine::a = 0.447;
double UpcCalcMachine::Z = 82;
double UpcCalcMachine::sqrts = 5020;
double UpcCalcMachine::g1 = sqrts / (2. * mProt);
double UpcCalcMachine::g2 = sqrts / (2. * mProt);
int UpcCalcMachine::debug = 0;
double* UpcCalcMachine::vCachedFormFac = new double[UpcCalcMachine::nQ2];

UpcCalcMachine::UpcCalcMachine() = default;

UpcCalcMachine::~UpcCalcMachine() = default;

void UpcCalcMachine::init()
{
  // update scaling factor
  factor = Z * Z * alpha / M_PI / M_PI / hc / hc;

  // calculate Woods-Saxon rho0 from R and a
  rho0 = calcWSRho();

  // prepare caches
  prepareGAA();           // G_AA and Fourier-transformed G_AA
  prepareFormFac();       // nuclear form factor
  prepareTwoPhotonLumi(); // calculate two-photon luminosity and save into a file
}

template <typename ArrayType>
double UpcCalcMachine::simpson(int n, ArrayType* v, double h)
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

double UpcCalcMachine::calcWSRho()
{
  for (int ib = 0; ib < nb; ib++) {
    double r = ib * db;
    vRho[ib] = r * r / (1 + exp((r - R) / a));
  }
  double wsRho0 = A / simpson(nb, vRho, db) / 4 / M_PI;
  return wsRho0;
}

double UpcCalcMachine::fluxPoint(const double b, const double k)
{
  // flux divided by k
  double g = g1;
  double x = b * k / g / hc;
  double K0 = x > 1e-10 ? TMath::BesselK0(x) : 0;
  double K1 = x > 1e-10 ? TMath::BesselK1(x) : 0;
  double result = factor * k / g / g * (K1 * K1 + K0 * K0 / g / g);
  if (debug > 1) {
    PLOG_DEBUG << "result = " << result;
  }
  return result;
}

double UpcCalcMachine::fluxFormInt(double* x, double* par)
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
  double result = k * k * ff / t * TMath::BesselJ1(b * k / hc);
  if (debug > 1) {
    PLOG_DEBUG << "result = " << result;
  }
  return result;
}

double UpcCalcMachine::fluxForm(const double b, const double k, TF1* fFluxFormInt)
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

double UpcCalcMachine::calcTwoPhotonLumi(double M, double Y, TF1* fFluxForm, const TGraph* gGAA)
{
  // double differential luminosity
  double k1 = M / 2. * exp(Y);
  double k2 = M / 2. * exp(-Y);

  double b1min = isPoint ? 1 * R : 0.05 * R;
  double b2min = isPoint ? 1 * R : 0.05 * R;
  double b1max = TMath::Max(5. * g1 * hc / k1, 5 * R);
  double b2max = TMath::Max(5. * g2 * hc / k2, 5 * R);
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
        double gaa = b < 20. ? gGAA->Eval(b) : 1;
        sum_phi += gaa * weights10[k];
      }
      sum_b2 += flux[j] * sum_phi * b2 * (b2h - b2l);
    }
    sum += fluxForm(b1, k1, fFluxForm) * sum_b2 * b1 * (b1h - b1l);
  }
  double lumi = 2 * M_PI * M_PI * M * sum;
  return lumi;
}

double UpcCalcMachine::calcCrossSectionMPolS(double m)
{
  double r = 2 * mLep / m;
  double r2 = r * r;
  double r4 = r2 * r2;
  if (r > 1) {
    return 0;
  }
  double cs_s = 4 * M_PI * alpha * alpha * hc * hc / m / m * ((1 + r * r - 3. / 4. * r * r * r * r) * 2 * log(1 / r + sqrt(1 / r / r - 1)) - (1 + 3. / 2. * r * r) * sqrt(1 - r * r));
  return cs_s;
  // fm^2 //
}

double UpcCalcMachine::calcCrossSectionMZPolS(double m, double z)
{
  double mLep2 = mLep * mLep;
  double m2 = m * m;
  double z2 = z * z;
  double cs_s = 2 * M_PI * alpha * alpha;
  cs_s *= m2 - 4 * mLep2;
  cs_s *= sqrt(m2 - 4 * mLep2);
  cs_s *= (4 * mLep2 * (3 - 2 * z2 + z2 * z2)) + m2 * (1 - z2 * z2);
  cs_s /= m2 * m * (m2 * (1 - z2) + 4 * mLep2 * z2) * (m2 * (1 - z2) + 4 * mLep2 * z2);
  return cs_s;
  // GeV^-2 //
}

double UpcCalcMachine::calcCrossSectionMPolPS(double m)
{
  double r = 2 * mLep / m;
  double r2 = r * r;
  double r4 = r2 * r2;
  if (r > 1) {
    return 0;
  }
  return 4 * M_PI * alpha * alpha * hc * hc / m / m * ((1 + r * r - 1. / 4. * r * r * r * r) * 2 * log(1 / r + sqrt(1 / r / r - 1)) - (1 + 1. / 2. * r * r) * sqrt(1 - r * r));
  // fm^2 //
}

double UpcCalcMachine::calcCrossSectionMZPolPS(double m, double z)
{
  double mLep2 = mLep * mLep;
  double m2 = m * m;
  double z2 = z * z;
  double cs_s = 2 * M_PI * alpha * alpha;
  cs_s *= sqrt(m2 - 4 * mLep2);
  cs_s *= m2 * m2 * (1 - z2 * z2) + 8 * m2 * mLep2 * (1 - z2 + z2 * z2) - 16 * mLep2 * mLep2 * (1 - z2) * (1 - z2);
  cs_s /= m2 * m * (m2 * (1 - z2) + 4 * mLep2 * z2) * (m2 * (1 - z2) + 4 * mLep2 * z2);
  return cs_s;
  // GeV^-2 //
}

void UpcCalcMachine::calcTwoPhotonLumiPol(double& ns, double& np, double M, double Y, TF1* fFluxForm, const TGraph* gGAA)
{
  double k1 = M / 2. * exp(Y);
  double k2 = M / 2. * exp(-Y);

  double b1min = isPoint ? 1 * R : 0.05 * R;
  double b2min = isPoint ? 1 * R : 0.05 * R;
  double b1max = TMath::Max(5. * g1 * hc / k1, 5 * R);
  double b2max = TMath::Max(5. * g2 * hc / k2, 5 * R);
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
        double gaa = b < 20 ? gGAA->Eval(b) : 1;
        sum_phi_s += gaa * weights10[k] * cphi * cphi;
        sum_phi_p += gaa * weights10[k] * sphi * sphi;
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

double UpcCalcMachine::calcCrossSectionMZ(double m, double z)
{
  double s = m * m;
  double k = TMath::Sqrt(s) / 2.;              // photon/lepton energy in cm system in GeV
  double p = TMath::Sqrt(k * k - mLep * mLep); // outgoing lepton momentum in GeV
  double norm = 2 * M_PI * alpha * alpha / s * p / k;
  double kt = -2 * k * (k - z * p) / mLep / mLep;
  double ku = -2 * k * (k + z * p) / mLep / mLep;
  double ks = kt + ku;
  double kp = kt * ku;
  double kq = 1. / kt + 1. / ku;
  double kr = ku / kt + kt / ku;
  double cs = 0;
  cs += -8. * (4. * kq * kq + 4. * kq - kr);
  cs += 16. * (2. + kr) * aLep;
  cs += 4. * (2. - 4. * ks + kr) * aLep * aLep;
  cs += -8. * (2. + 2. * ks + kr) * aLep * aLep * aLep;
  cs += -4. * (4. + 2. * ks + 2. * kr - kp) * aLep * aLep * aLep * aLep;
  cs *= norm;
  return cs; // [GeV^-2]
}

double UpcCalcMachine::calcCrossSectionM(double m)
{
  double s = m * m;               // cms invariant mass squared
  double x = 4 * mLep * mLep / s; // inverse lepton gamma-factor squared = 1/g^2 in cms
  double b = TMath::Sqrt(1 - x);  // lepton relativistic velocity in cms
  double y = atanh(b);            // lepton rapidity in cms
  double cs = 0;
  cs += (2 + 2 * x - x * x) * y - b * (1 + x);
  cs += 4 * y * aLep;
  cs += (4 * b / x + y) * aLep * aLep;
  cs += (4 * b / x - 2 * y) * aLep * aLep * aLep;
  cs += ((7. / 12.) * b / x + (1. / 6.) * b / x / x - 0.5 * y) * aLep * aLep * aLep * aLep;
  cs *= 4 * hc * hc * 1e7 * alpha * alpha * M_PI / s;
  return cs;
  // [nb] //
}

void UpcCalcMachine::fillCrossSectionMZ(TH2D* hCrossSectionMZ,
                                        double mmin, double mmax, int nm,
                                        double zmin, double zmax, int nz,
                                        int flag)
{
  double m, z;
  double dm = (mmax - mmin) / nm;
  double dz = (zmax - zmin) / nz;
  double cs;
  for (int im = 0; im <= nm; im++) {
    m = mmin + dm * im;
    for (int iz = 0; iz <= nz; iz++) {
      z = zmin + dz * iz;
      if (flag == 0) { // the usual unpolarized cross section
        cs = calcCrossSectionMZ(m, z);
      }
      if (flag == 1) { // scalar part
        cs = calcCrossSectionMZPolS(m, z);
      }
      if (flag == 2) { // pseudoscalar part
        cs = calcCrossSectionMZPolPS(m, z);
      }
      hCrossSectionMZ->SetBinContent(im, iz, cs);
    }
  }
  double scalingFactor = hc * hc * 1e7; // to [nb]
  hCrossSectionMZ->Scale(scalingFactor / dm);
}

void UpcCalcMachine::fillCrossSectionM(TH1D* hCrossSectionM,
                                       double mmin, double mmax, int nm)
{
  double m;
  double dm = (mmax - mmin) / nm;
  double cs;
  for (int im = 0; im <= nm; im++) {
    m = mmin + dm * im;
    cs = calcCrossSectionM(m);
    hCrossSectionM->SetBinContent(im, cs);
  }
  double scalingFactor = hc * hc * 1e7; // to [nb]
  hCrossSectionM->Scale(scalingFactor);
}

void UpcCalcMachine::prepareGAA()
{
  double ssm = pow(sqrts, 2) / pow(2 * mProt + 2.1206, 2);
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
}

double UpcCalcMachine::calcFormFac(double Q2)
{
  double Q = sqrt(Q2) / hc;
  double coshVal = TMath::CosH(M_PI * Q * a);
  double sinhVal = TMath::SinH(M_PI * Q * a);
  double ff = 4 * M_PI * M_PI * rho0 * a * a * a / (Q * a * Q * a * sinhVal * sinhVal) *
              (M_PI * Q * a * coshVal * TMath::Sin(Q * R) - Q * R * TMath::Cos(Q * R) * sinhVal);
  ff += 8 * M_PI * rho0 * a * a * a * exp(-R / a) / (1 + Q * Q * a * a) / (1 + Q * Q * a * a);
  return ff;
}

void UpcCalcMachine::prepareFormFac()
{
  for (int iQ2 = 0; iQ2 < nQ2; iQ2++) {
    double Q2 = Q2min + iQ2 * dQ2;
    double ff = calcFormFac(Q2);
    vCachedFormFac[iQ2] = ff;
  }
}

double UpcCalcMachine::getCachedFormFac(double Q2)
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

void UpcCalcMachine::prepareTwoPhotonLumi()
{
  bool isFound = false;
  if (!gSystem->AccessPathName("twoPhotonLumi.root") && !usePolarizedCS) {
    PLOG_INFO << "Found pre-calculated unpolarized 2D luminosity";
    isFound = true;
  }
  if (!gSystem->AccessPathName("twoPhotonLumiPol.root") && usePolarizedCS) {
    PLOG_INFO << "Found pre-calculated polarized 2D luminosity";
    isFound = true;
  }
  if (!isFound) {
    PLOG_INFO << "Precalculated 2D luminosity is not found. Starting all over...";
    double dy = (ymax - ymin) / (ny - 1);
    double dm = (mmax - mmin) / (nm - 1);
    TString fname = usePolarizedCS ? "twoPhotonLumiPol.root" : "twoPhotonLumi.root";
    auto* f2DLumi = new TFile(fname, "recreate");
    // histograms for unpolarized case
    TH2D* hD2LDMDY;
    // histograms for polarized case
    TH2D* hD2LDMDY_s;
    TH2D* hD2LDMDY_p;
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
      TH2D* hD2LDMDY_private;
      TH2D* hD2LDMDY_private_s;
      TH2D* hD2LDMDY_private_p;
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
}

void UpcCalcMachine::calcNucCrossSectionYM(TH2D* hCrossSectionYM, TH2D* hPolCSRatio)
{
  PLOG_INFO << "Calculating nuclear cross section for a_lep = " << aLep;

  double dy = (ymax - ymin) / (ny - 1);
  double dm = (mmax - mmin) / (nm - 1);

  TH2D* hD2LDMDY;
  TH2D* hD2LDMDY_s;
  TH2D* hD2LDMDY_p;

  TString fname = usePolarizedCS ? "twoPhotonLumiPol.root" : "twoPhotonLumi.root";
  auto* f2DLumi = new TFile(fname, "r");

  if (usePolarizedCS) {
    if (!hPolCSRatio) {
      PLOG_FATAL << "hPolRatio is not initialized!";
      std::_Exit(-1);
    }
    hD2LDMDY_s = (TH2D*)f2DLumi->Get("hD2LDMDY_s");
    hD2LDMDY_p = (TH2D*)f2DLumi->Get("hD2LDMDY_p");
  } else { // loading pre-cached two-photon luminosity
    hD2LDMDY = (TH2D*)f2DLumi->Get("hD2LDMDY");
  }

  // calculating nuclear cross section
  PLOG_INFO << "Calculating nuclear cross section...";
  ROOT::EnableThreadSafety();
  int im, iy, ib;
  int progress = 0;
  int total = nm * ny;
  double cs[nm][ny];
  double cs_rat[nm][ny];
  omp_set_num_threads(numThreads);
#pragma omp parallel default(none)                                                   \
  shared(cs, cs_rat, hD2LDMDY, hD2LDMDY_s, hD2LDMDY_p, progress) private(im, iy, ib) \
    firstprivate(nb, vb, total, numThreads, dm, dy, ymin, ymax, mmin, mmax, nm, ny, abscissas10, weights10, vGAA)
  {
    vector<vector<double>> cs_private(nm, vector<double>(ny, 0));
    vector<vector<double>> rat_private(nm, vector<double>(ny, 0));
    TH2D* hD2LDMDY_private;
    TH2D* hD2LDMDY_private_s;
    TH2D* hD2LDMDY_private_p;
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
          cs_private[im][iy] = calcCrossSectionM(m) * lumi;
        } else { // polarized
          double cs_s = calcCrossSectionMPolS(m);
          double cs_p = calcCrossSectionMPolPS(m);
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
  }
  omp_set_num_threads(1);
  PLOG_INFO << "Calculating nuclear cross section...Done!";

  // filling histograms
  double cssum = 0;
  for (int i = 0; i < nm - 1; i++) {
    for (int j = 0; j < ny - 1; j++) {
      double cs_ij = (cs[i][j] + cs[i + 1][j] + cs[i][j + 1] + cs[i + 1][j + 1]) / 4.;
      cssum += cs_ij;
      hCrossSectionYM->SetBinContent(j + 1, i + 1, cs_ij);
      if (usePolarizedCS) {
        double rat_ij = (cs_rat[i][j] + cs_rat[i + 1][j] + cs_rat[i][j + 1] + cs_rat[i + 1][j + 1]) / 4.;
        hPolCSRatio->SetBinContent(j + 1, i + 1, rat_ij);
      }
    }
  }

  hCrossSectionYM->Scale(1e-6); // nb -> mb

  PLOG_INFO << "Total nuclear cross section = " << cssum * 1e-6 << " mb";
}

// Ref.: S.R.Klein, J.Nystrand, PRC 60 014903, 1999
double UpcCalcMachine::calcNucFormFactor(double t)
{
  double ffactor;
  if (Z < 7) {
    ffactor = exp(-R * R * t * t / (6 * hc * hc));
  }
  if (Z >= 7) {
    double q = TMath::Sqrt(t);
    double qR = q * R / hc;
    double invqR = hc / (q * R);
    ffactor = (TMath::Sin(qR) - qR * TMath::Cos(qR)) * 3. * invqR * invqR * invqR;
    const double a0 = 0.7; // [fm]
    ffactor = ffactor / (1. + (a0 * a0 * t) / (hc * hc));
  }
  return ffactor;
}

// Function from Starlight
// (by S.R.Klein, J.Nystrand, J.Seger, Y.Gorbunov, J.Butterworth)
double UpcCalcMachine::getPhotonPt(double ePhot)
{
  constexpr double pi2x4 = 4 * M_PI * M_PI;
  double y1 = TMath::ACosH(g1);
  double y2 = -TMath::ACosH(g2);
  double gtot = TMath::CosH((y1 - y2) / 2.);

  double ereds = (ePhot / gtot) * (ePhot / gtot);
  double Cm = TMath::Sqrt(3.) * ePhot / gtot;
  double arg = Cm * Cm + ereds;
  double sFFactCM = calcNucFormFactor(arg);
  double Coef = 3. * (sFFactCM * sFFactCM * Cm * Cm * Cm) / (pi2x4 * arg * arg);

  double x = gRandom->Uniform(0, 1);
  double pp = x * 5. * hc / R;
  arg = pp * pp + ereds;
  double sFFactPt1 = calcNucFormFactor(arg);
  double test = (sFFactPt1 * sFFactPt1) * pp * pp * pp / (pi2x4 * arg * arg);

  bool satisfy = false;
  while (!satisfy) {
    double u = gRandom->Uniform(0, 1);
    if (u * Coef < test) {
      satisfy = true;
    } else {
      x = gRandom->Uniform(0, 1);
      pp = x * 5 * hc / R;
      arg = pp * pp + ereds;
      double sFFactPt2 = calcNucFormFactor(arg);
      test = (sFFactPt2 * sFFactPt2) * pp * pp * pp / (pi2x4 * arg * arg);
    }
  }

  return pp;
}

void UpcCalcMachine::getPairMomentum(double mPair, double yPair, TLorentzVector& pPair)
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
