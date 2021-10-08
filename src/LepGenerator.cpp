//
// Created by nburmaso on 6/25/21.
//

#include "LepGenerator.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"
#include "TRandomGen.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TSystem.h"
#include <fstream>
#include <omp.h>
#include <plog/Log.h>
#include <plog/Init.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Appenders/ColorConsoleAppender.h>

using namespace std;

// out-of-line initialization for static members
double LepGenerator::rho0 = 0.159538;
double LepGenerator::R = 6.68;
double LepGenerator::a = 0.447;
int LepGenerator::debug = 0;

LepGenerator::LepGenerator()
{
  // get parameters from file
  initGeneratorFromFile();

  // initialize the MT64 random number generator
  gRandom = new TRandomMT64();
  gRandom->SetSeed(time(nullptr));

  // initialize pythia for decays
  PLOG_INFO << "Initializing Pythia-based decayer";
  if (pythiaVersion == 8) {
    decayer = new TPythia8Decayer();
  } else if (pythiaVersion == 6) {
    decayer = new TPythia6Decayer();
  } else {
    PLOG_WARNING << "Wrong Pythia version! Please choose either 8 or 6";
    PLOG_WARNING << "Using Pythia8 by default...";
    decayer = new TPythia8Decayer();
  }

  decayer->Init();
}

LepGenerator::~LepGenerator() = default;

double LepGenerator::simpson(int n, double* v, double h)
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

double LepGenerator::fluxPoint(const double b, const double w, const double g)
{
  // flux divided by w
  double x = b * w / g / hc;
  double K0 = x > 1e-10 ? TMath::BesselK0(x) : 0;
  double K1 = x > 1e-10 ? TMath::BesselK1(x) : 0;
  double result = factor * w / g / g * (K1 * K1 + K0 * K0 / g / g);
  if (debug > 1) {
    PLOG_DEBUG << "result = " << result;
  }
  return result;
}

double LepGenerator::fluxFormInt(double* x, double* par)
{
  double k = x[0];
  double b = par[0];
  double w = par[1];
  double g = par[2];

  if (debug > 1) {
    PLOG_DEBUG << "b = " << b << " w = " << w << " g = " << g;
  }

  double t = k * k + w * w / g / g;
  double q = sqrt(t) / hc;
  double sinhVal = sinh(M_PI * q * a);
  double coshVal = cosh(M_PI * q * a);
  double ff = 4 * M_PI * M_PI * rho0 * a * a * a / (q * a * q * a * sinhVal * sinhVal) *
              (M_PI * q * a * coshVal * sin(q * R) - q * R * cos(q * R) * sinhVal);
  for (int n = 1; n < 2; n++) {
    ff += 8 * M_PI * rho0 * a * a * a * (n % 2 ? 1 : -1) * n * exp(-n * R / a) / (n * n + q * q * a * a) / (n * n + q * q * a * a);
  }
  double result = k * k * ff / t * TMath::BesselJ1(b * k / hc);
  if (debug > 1) {
    PLOG_DEBUG << "result = " << result;
  }
  return result;
}

// wrapper for gsl gauss integrator
double LepGenerator::fluxForm(const double b, const double w, const double g)
{
  // flux divided by w
  if (isPoint) {
    return fluxPoint(b, w, g);
  }

  if (b > 2 * R) {
    return fluxPoint(b, w, g);
  }

  fFluxFormInt->SetParameter(0, b);
  fFluxFormInt->SetParameter(1, w);
  fFluxFormInt->SetParameter(2, g);
  double Q = fFluxFormInt->Integral(0, 10, 1e-5) / A;
  double result = factor * Q * Q / w;

  if (debug > 1) {
    PLOG_DEBUG << "result = " << result;
  }

  return result;
}

double LepGenerator::D2LDMDY(double M, double Y)
{
  // double differential luminosity
  double D2LDMDYx = 0.;
  double w1 = M / 2. * exp(Y);
  double w2 = M / 2. * exp(-Y);

  double b1min = 0.05 * R;
  double b2min = 0.05 * R;

  if (isPoint) {
    b1min = 1 * R;
  }
  if (isPoint) {
    b2min = 1 * R;
  }

  double b1max = TMath::Max(5. * g1 * hc / w1, 5 * R);
  double b2max = TMath::Max(5. * g2 * hc / w2, 5 * R);
  double log_delta_b1 = (log(b1max) - log(b1min)) / nb1;
  double log_delta_b2 = (log(b2max) - log(b2min)) / nb2;

  double flux[nb2];
  for (int j = 0; j < nb2; ++j) {
    double b2l = b2min * exp(j * log_delta_b2);
    double b2h = b2min * exp((j + 1) * log_delta_b2);
    double b2 = (b2h + b2l) / 2.;
    flux[j] = fluxForm(b2, w2, g2);
  }

  double sum = 0;
  double sum_b2, b1l, b1h, b1, b2l, b2h, b2, sum_phi, b;
  int i, j, k;

  for (i = 0; i < nb1; ++i) {
    sum_b2 = 0.;
    b1l = b1min * exp(i * log_delta_b1);
    b1h = b1min * exp((i + 1) * log_delta_b1);
    b1 = (b1h + b1l) / 2.;
    for (j = 0; j < nb2; ++j) {
      b2l = b2min * exp(j * log_delta_b2);
      b2h = b2min * exp((j + 1) * log_delta_b2);
      b2 = (b2h + b2l) / 2.;
      sum_phi = 0.;
      for (k = 0; k < ngi; k++) {
        if (abscissas[k] < 0) {
          continue;
        }
        b = sqrt(b1 * b1 + b2 * b2 + 2. * b1 * b2 * cos(M_PI * abscissas[k]));
        sum_phi += (b < 20.) ? (weights[k] * gGAA->Eval(b) * 2) : (weights[k] * 2.);
      }
      sum_b2 += flux[j] * M_PI * sum_phi * b2 * (b2h - b2l);
    }
    sum += fluxForm(b1, w1, g1) * sum_b2 * b1 * (b1h - b1l);
  }
  D2LDMDYx = M_PI * M * sum;
  return D2LDMDYx;
}

double LepGenerator::crossSectionWZ(double s, double z)
{
  double k = sqrt(s / 2.);              // photon/lepton energy in cm system in GeV
  double p = sqrt(k * k - mLep * mLep); // outgoing lepton momentum in GeV
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

double LepGenerator::crossSectionW(double w)
{
  double s = w * w;               // cms invariant mass squared
  double x = 4 * mLep * mLep / s; // inverse lepton gamma-factor squared = 1/g^2 in cms
  double b = sqrt(1 - x);         // lepton relativistic velocity in cms
  double y = atanh(b);            // lepton rapidity in cms
  double cs = 0;
  cs += (2 + 2 * x - x * x) * y - b * (1 + x);
  cs += 4 * y * aLep;
  cs += (4 * b / x + y) * aLep * aLep;
  cs += (4 * b / x - 2 * y) * aLep * aLep * aLep;
  cs += ((7. / 12.) * b / x + (1. / 6.) * b / x / x - 0.5 * y) * aLep * aLep * aLep * aLep;
  cs *= 4 * alpha * alpha * M_PI / s;
  return cs; // [GeV^-2]
}

void LepGenerator::fillCrossSectionWZ(TH2D* hCrossSectionWZ,
                                      double wmin, double wmax, int nw,
                                      double zmin, double zmax, int nz)
{
  double w, z;
  double dw = (wmax - wmin) / nw;
  double dz = (zmax - zmin) / nz;
  double cs;
  for (int iw = 0; iw <= nw; iw++) {
    w = wmin + dw * iw;
    for (int iz = 0; iz <= nz; iz++) {
      z = zmin + dz * iz;
      cs = crossSectionWZ(w * w, z);
      hCrossSectionWZ->SetBinContent(iw, iz, cs);
    }
  }
  double scalingFactor = hc * hc * 1e7; // to [nb]
  hCrossSectionWZ->Scale(scalingFactor / dw);
}

void LepGenerator::fillCrossSectionW(TH1D* hCrossSectionW,
                                     double wmin, double wmax, int nw)
{
  double w;
  double dw = (wmax - wmin) / nw;
  double cs;
  for (int iw = 0; iw <= nw; iw++) {
    w = wmin + dw * iw;
    cs = crossSectionW(w);
    hCrossSectionW->SetBinContent(iw, cs);
  }
  double scalingFactor = hc * hc * 1e7; // to [nb]
  hCrossSectionW->Scale(scalingFactor);
}

void LepGenerator::nuclearCrossSectionYM(TH2D* hCrossSectionYM)
{
  PLOG_INFO << "Calculating nuclear cross section for a_lep = " << aLep;

  double dy = (ymax - ymin) / (ny - 1);
  double dw = (wmax - wmin) / (nw - 1);

  // calculating two-photon luminosity (if needed)
  // -----------------------------------------------------------------------

  if (!gSystem->AccessPathName("hD2LDMDY.root")) {
    PLOG_INFO << "Found precalculated 2D luminosity";
  } else {
    PLOG_INFO << "Precalculated 2D luminosity is not found. Starting all over...";

    // -----------------------------------------------------------------------

    double ssm = pow(sqrts, 2) / pow(2 * mProt + 2.1206, 2);
    double csNN = 0.1 * (34.41 + 0.2720 * pow(log(ssm), 2) + 13.07 * pow(ssm, -0.4473) - 7.394 * pow(ssm, -0.5486)); // PDG 2016

    // calculate rho and TA
    double b, z, r, s, sum_phi;
    int ib, is, iz, k;
    double TAb[nb];
    for (ib = 0; ib < nb; ib++) {
      b = ib * db;
      for (iz = 0; iz < nb; iz++) {
        z = iz * db;
        r = sqrt(b * b + z * z);
        rho[ib][iz] = rho0 / (1 + exp((r - R) / a));
      }
      TA[ib] = 2 * simpson(nb, rho[ib], db);
      TAb[ib] = TA[ib] * b;
      vb[ib] = b;
    }
    auto* gTA = new TGraph(nb, vb, TA);

    // calculate GAA
    for (ib = 0; ib < nb; ib++) {
      b = ib * db;
      for (is = 0; is < nb; is++) {
        s = is * db;
        sum_phi = 0;
        for (k = 0; k < ngi; k++) {
          if (abscissas[k] < 0)
            continue;
          r = sqrt(b * b + s * s + 2 * b * s * cos(M_PI * abscissas[k]));
          sum_phi += 2 * M_PI * weights[k] * gTA->Eval(r);
        }
        vs[is] = 2 * s * gTA->Eval(s) * sum_phi;
      }
      vGAA[ib] = exp(-csNN * simpson(nb, vs, db));
    }

    fFluxFormInt = new TF1("fFluxFormInt", fluxFormInt, 0, 10, 3);
    gGAA = new TGraph(nb, vb, vGAA);

    auto* f2DLumi = new TFile("hD2LDMDY.root", "recreate");
    auto* hD2LDMDY = new TH2D("hD2LDMDY", ";;", nw, wmin, wmax, ny, ymin, ymax);
    PLOG_INFO << "Calculating 2D luminosity grid..." << 0 << "%";
    for (Int_t iw = 0; iw < nw; iw++) {
      double M = wmin + dw * iw;
      for (Int_t iy = 0; iy < ny; iy++) {
        double y = ymin + dy * iy;
        hD2LDMDY->SetBinContent(iw, iy, D2LDMDY(M, y) * dw * dy);
      }
      double progressBar = 100. * (double)(iw + 1) / nw;
      PLOG_INFO << "Calculating 2D luminosity grid..." << fixed << setprecision(2) << progressBar << "%";
    }
    hD2LDMDY->Write();
    f2DLumi->Close();

    PLOG_INFO << "Two-photon luminosity was written to 'hD2LDMDY.root'";
  }

  // calculating nuclear cross section
  // -----------------------------------------------------------------------

  // loading pre-cached two-photon luminosity
  auto* f2DLumi = new TFile("hD2LDMDY.root", "r");
  auto* hD2LDMDY = (TH2D*)f2DLumi->Get("hD2LDMDY");

  // calculating total elementary cross section
  PLOG_INFO << "Calculating total elementary cross section...";
  auto* hCrossSectionW = new TH1D("hCrossSectionW", ";w [gev]; cs [nb];",
                                  nw, wmin, wmax);
  fillCrossSectionW(hCrossSectionW, wmin, wmax, nw);

  PLOG_INFO << "Calculating total elementary cross section...Done!";

  auto* hCrossSectionW_ax = hCrossSectionW->GetXaxis();

  // calculating nuclear cross section
  PLOG_INFO << "Calculating nuclear cross section...";
  double cs[nw][ny];
  for (int iw = 0; iw < nw; iw++) {
    double M = wmin + dw * iw;
    for (int iy = 0; iy < ny; iy++) {
      double csItem = hCrossSectionW->GetBinContent(hCrossSectionW_ax->FindBin(M));
      cs[iw][iy] = csItem * 1e-6 * hD2LDMDY->GetBinContent(iw, iy);
    }
  }
  PLOG_INFO << "Calculating nuclear cross section...Done!";

  // filling a histogram
  double cssum = 0;
  for (int i = 0; i < nw - 1; i++) {
    for (int j = 0; j < ny - 1; j++) {
      double cs_ij = (cs[i][j] + cs[i + 1][j] + cs[i][j + 1] + cs[i + 1][j + 1]) / 4.;
      cssum += cs_ij;
      hCrossSectionYM->SetBinContent(j + 1, i + 1, cs_ij);
    }
  }

  PLOG_INFO << "Total nuclear cross section = " << cssum << " mb";
}

// Ref.: S.R.Klein, J.Nystrand, PRC 60 014903, 1999
double LepGenerator::nucFormFactor(double q)
{
  double ffactor;
  if (Z < 7) {
    ffactor = exp(-R * R * q * q / (6 * hc * hc));
  }
  if (Z >= 7) {
    double qR = q * R / hc;
    double invqR = hc / (q * R);
    ffactor = (sin(qR) - qR * cos(qR)) * 3. * invqR * invqR * invqR;
    ffactor = ffactor / (1. + (a * a * q * q) / (hc * hc));
  }
  return ffactor;
}

double LepGenerator::getPhotonPt(double ePhot)
{
  double y1 = acosh(g1);
  double y2 = -acosh(g2);
  double gtot = cosh((y1 - y2) / 2.);

  double ereds = (ePhot / gtot) * (ePhot / gtot);
  // sqrt(3) * E / gamma_em is p_t where the distribution is a maximum
  double Cm = sqrt(3.) * ePhot / gtot;
  // the amplitude of the p_t spectrum at the maximum
  double sFFactCM = nucFormFactor(Cm * Cm + ereds);
  double Coef = 3. * (sFFactCM * sFFactCM * Cm * Cm * Cm) / (4. * M_PI * (ereds + Cm * Cm) * M_PI * (ereds + Cm * Cm));

  // pick a test value pp, and find the amplitude there
  double x = gRandom->Uniform(0, 1);
  double pp = x * 5. * hc / R;
  double sFFactPt1 = nucFormFactor(pp * pp + ereds);
  double test = (sFFactPt1 * sFFactPt1) * pp * pp * pp / (4. * M_PI * (ereds + pp * pp) * M_PI * (ereds + pp * pp));

  bool satisfy = false;
  while (!satisfy) {
    double u = gRandom->Uniform(0, 1);
    if (u * Coef < test) {
      satisfy = true;
    } else {
      x = gRandom->Uniform(0, 1);
      pp = 5 * hc / R * x;
      double sFFactPt2 = nucFormFactor(pp * pp + ereds);
      test = (sFFactPt2 * sFFactPt2) * pp * pp * pp / (4. * M_PI * (ereds + pp * pp) * M_PI * (ereds + pp * pp));
    }
  }

  return pp;
}

void LepGenerator::getPairMomentum(double mPair, double yPair, TLorentzVector& pPair)
{
  if (!useNonzeroGamPt) {
    double mtPair = sqrt(mPair * mPair); // pairPt = 0
    pPair.SetPxPyPzE(0., 0., mtPair * sinh(yPair), mtPair * cosh(yPair));
  }
  if (useNonzeroGamPt) {
    double w1 = mPair / 2 * exp(yPair);
    double w2 = mPair / 2 * exp(-yPair);
    double angle1 = gRandom->Uniform(0, 2 * M_PI);
    double angle2 = gRandom->Uniform(0, 2 * M_PI);
    double pt1 = getPhotonPt(w1);
    double pt2 = getPhotonPt(w2);
    double px = pt1 * cos(angle1) + pt2 * cos(angle2);
    double py = pt1 * sin(angle1) + pt2 * sin(angle2);
    double pt = sqrt(px * px + py * py);
    double mtPair = sqrt(mPair * mPair + pt * pt);
    double pz = mtPair * sinh(yPair);
    double e = mtPair * cosh(yPair);
    pPair.SetPxPyPzE(px, py, pz, e);
  }
}

void LepGenerator::initGeneratorFromFile()
{
  // todo: use <any> from c++17 for a neat parsing???
  if (!gSystem->AccessPathName("parameters.in")) {
    PLOG_INFO << "Reading parameters from parameters.in ...";
    InputPars parDict;
    ifstream fInputs("parameters.in");
    string line;
    string parameter;
    string parValue;
    while (getline(fInputs, line)) {
      istringstream iss(line);
      iss >> parameter >> parValue;
      if (parameter == parDict.inNEvents) {
        nEvents = stoi(parValue);
      }
      if (parameter == parDict.inCMSqrtS) {
        sqrts = stod(parValue);
        g1 = sqrts / (2. * mProt);
        g2 = sqrts / (2. * mProt);
      }
      if (parameter == parDict.inLepPDG) {
        lepPDG = stoi(parValue);
      }
      if (parameter == parDict.inLepM) {
        mLep = stod(parValue);
      }
      if (parameter == parDict.inLepA) {
        aLep = stod(parValue);
      }
      if (parameter == parDict.inDoPtCut) {
        doPtCut = stoi(parValue);
      }
      if (parameter == parDict.inLowPt) {
        minPt = stod(parValue);
      }
      if (parameter == parDict.inLowZ) {
        zmin = stod(parValue);
      }
      if (parameter == parDict.inHiZ) {
        zmax = stod(parValue);
      }
      if (parameter == parDict.inLowW) {
        wmin = stod(parValue);
      }
      if (parameter == parDict.inHiW) {
        wmax = stod(parValue);
      }
      if (parameter == parDict.inLowY) {
        ymin = stod(parValue);
      }
      if (parameter == parDict.inHiY) {
        ymax = stod(parValue);
      }
      if (parameter == parDict.inBinsZ) {
        nz = stoi(parValue);
      }
      if (parameter == parDict.inBinsW) {
        nw = stoi(parValue);
      }
      if (parameter == parDict.inBinsY) {
        ny = stoi(parValue);
      }
      if (parameter == parDict.inWSRho0) {
        rho0 = stod(parValue);
      }
      if (parameter == parDict.inWSRadius) {
        R = stod(parValue);
      }
      if (parameter == parDict.inWSA) {
        a = stod(parValue);
      }
      if (parameter == parDict.inNucZ) {
        Z = stod(parValue);
      }
      if (parameter == parDict.inNucA) {
        A = stod(parValue);
      }
      if (parameter == parDict.inFluxPoint) {
        isPoint = stoi(parValue);
      }
      if (parameter == parDict.inPythiaVer) {
        pythiaVersion = stoi(parValue);
      }
      if (parameter == parDict.inNonzeroGamPt) {
        useNonzeroGamPt = stoi(parValue);
      }
    }
    fInputs.close();
  } else {
    PLOG_WARNING << "Input file not found! Using default parameters...";
    printParameters();
  }
  // update scaling factor
  factor = Z * Z * alpha / M_PI / M_PI / hc / hc;
}

void LepGenerator::printParameters()
{
  PLOG_WARNING << "NEVENTS " << nEvents;
  PLOG_WARNING << "SQRTS " << sqrts;
  PLOG_WARNING << "LEP_PDG " << lepPDG;
  PLOG_WARNING << "LEP_MASS " << mLep;
  PLOG_WARNING << "LEP_A " << aLep;
  PLOG_WARNING << "DO_PT_CUT " << doPtCut;
  PLOG_WARNING << "PT_MIN " << minPt;
  PLOG_WARNING << "ZMIN " << zmin;
  PLOG_WARNING << "ZMAX " << zmax;
  PLOG_WARNING << "WMIN " << wmin;
  PLOG_WARNING << "WMAX " << wmax;
  PLOG_WARNING << "YMIN " << ymin;
  PLOG_WARNING << "YMAX " << ymax;
  PLOG_WARNING << "BINS_Z " << nz;
  PLOG_WARNING << "BINS_W " << nw;
  PLOG_WARNING << "BINS_Y " << ny;
  PLOG_WARNING << "WS_RHO0 " << rho0;
  PLOG_WARNING << "WS_RAD " << R;
  PLOG_WARNING << "WS_A " << a;
  PLOG_WARNING << "NUCLEUS_Z " << Z;
  PLOG_WARNING << "NUCLUES_A " << A;
  PLOG_WARNING << "FLUX_POINT" << isPoint;
  PLOG_WARNING << "NON_ZERO_GAM_PT" << useNonzeroGamPt;
  PLOG_WARNING << "PYTHIA_VERSION" << pythiaVersion;
}

void LepGenerator::generateEvents()
{
  // calculating elementary cross section in WZ space
  // -----------------------------------------------------------------------
  PLOG_INFO << "Calculating elementary cross section...";

  auto* hCrossSectionWZ = new TH2D("hCrossSectionWZ", ";w [gev]; z; cs [nb/gev]",
                                   nw, wmin, wmax,
                                   nz, zmin, zmax);

  fillCrossSectionWZ(hCrossSectionWZ, wmin, wmax, nw, zmin, zmax, nz);

  // calculating nuclear cross section in YM space
  // -----------------------------------------------------------------------
  auto* hNucCSYM = new TH2D("hNucCSYM", "", ny - 1, ymin, ymax, nw - 1, wmin, wmax);
  nuclearCrossSectionYM(hNucCSYM);

  if (debug > 0) {
    PLOG_DEBUG << "a_lep = " << aLep << ", min. pt = " << minPt;
  }

  double mPair = 0.;
  double yPair = 0.;
  double ptPair = 0.;
  double mtPair;
  double pMag;
  double cost, theta, phi;
  int binW;

  TLorentzVector pPair;
  vector<TLorentzVector> pLeps(2);
  TVector3 vLep;
  TVector3 boost;

  vector<double> cutsZ(hNucCSYM->GetNbinsY());
  if (doPtCut) {
    for (int mBin = 1; mBin <= hNucCSYM->GetNbinsY(); mBin++) {
      double mass = hNucCSYM->GetYaxis()->GetBinLowEdge(mBin);
      double sqt = sqrt(mass * mass * 0.25 - mLep * mLep);
      if (sqt <= minPt) {
        cutsZ[mBin - 1] = 0;
        for (int yBin = 1; yBin <= hNucCSYM->GetNbinsX(); yBin++) {
          hNucCSYM->SetBinContent(yBin, mBin, 0);
        }
        continue;
      }
      double minZ = sqrt(1 - minPt * minPt / (sqt * sqt));
      cutsZ[mBin - 1] = minZ;
    }
  }

  // initialize helper structure
  Particle particle{};

  // initialize file output
  PLOG_INFO << "Events will be written to " << Form("events_%.3f_%.0f.root", aLep, minPt);
  TFile outFile(Form("events_%.3f_%.0f.root", aLep, minPt), "recreate", "", 4 * 100 + 5); // using LZ4 with level 5 compression
  TTree outTree("particles", "Generated particles");
  outTree.Branch("eventNumber", &particle.eventNumber, "eventNumber/I");
  outTree.Branch("pdgCode", &particle.pdgCode, "pdgCode/I");
  outTree.Branch("particleID", &particle.particleID, "particleID/I");
  outTree.Branch("motherID", &particle.motherID, "motherID/I");
  outTree.Branch("px", &particle.px, "px/D");
  outTree.Branch("py", &particle.py, "py/D");
  outTree.Branch("pz", &particle.pz, "pz/D");
  outTree.Branch("e", &particle.e, "e/D");
  outTree.SetAutoSave(0);

  TClonesArray decayParticles("TParticle");
  TLorentzVector tlVector;
  vector<int> pdgLeps(2);

  PLOG_INFO << "Generating " << nEvents << " events...";

  for (int evt = 0; evt < nEvents; evt++) {
    if (debug > 1) {
      PLOG_DEBUG << "Event number: " << evt + 1;
    } else if ((evt + 1) % 100000 == 0) {
      PLOG_INFO << "Event number: " << evt + 1;
    }

    hNucCSYM->GetRandom2(yPair, mPair);

    getPairMomentum(mPair, yPair, pPair);
    pMag = sqrt(pPair.Pz() * pPair.Pz() + pPair.Pt() * pPair.Pt());

    binW = hCrossSectionWZ->GetXaxis()->FindBin(mPair);
    TH1D* hCSSliceAtW = hCrossSectionWZ->ProjectionY("sliceW", binW, binW);

    if (doPtCut) {
      int mBin = hNucCSYM->GetYaxis()->FindBin(mPair);
      double zCut = cutsZ[mBin - 1];
      int zCutBinUp = hCSSliceAtW->GetXaxis()->FindBin(zCut);
      int zCutBinLow = hCSSliceAtW->GetXaxis()->FindBin(-zCut);
      for (int zBin = zCutBinUp + 1; zBin <= hCSSliceAtW->GetNbinsX(); zBin++) {
        hCSSliceAtW->SetBinContent(zBin, 0);
      }
      for (int zBin = 1; zBin < zCutBinLow; zBin++) {
        hCSSliceAtW->SetBinContent(zBin, 0);
      }
    }

    cost = hCSSliceAtW->GetRandom();
    theta = acos(cost);
    phi = gRandom->Uniform(0., 2. * M_PI);

    vLep.SetMagThetaPhi(pMag, theta, phi);
    pLeps[0].SetVectM(vLep, mLep);
    pLeps[1].SetVectM(-vLep, mLep);

    boost = pPair.BoostVector();
    pLeps[0].Boost(boost);
    pLeps[1].Boost(boost);

    int sign = gRandom->Uniform(-1, 1) > 0 ? 1 : -1;
    pdgLeps[0] = sign * lepPDG;
    pdgLeps[1] = -sign * lepPDG;

    // lepton decays
    for (int iLep = 0; iLep < 2; iLep++) {
      int partID = 0;
      decayer->Decay(pdgLeps[iLep], &pLeps[iLep]);
      decayer->ImportParticles(&decayParticles);
      for (int ip = 0; ip < decayParticles.GetEntriesFast(); ip++) {
        auto* part = (TParticle*)decayParticles.At(ip);
        if (debug > 1) {
          PLOG_DEBUG << "Particle info:";
          part->Print();
        }
        int pdg = part->GetPdgCode();
        part->Momentum(tlVector);
        // filling output tree
        particle.eventNumber = evt;
        particle.pdgCode = pdg;
        particle.particleID = partID;
        particle.motherID = part->GetFirstMother();
        particle.px = tlVector[0];
        particle.py = tlVector[1];
        particle.pz = tlVector[2];
        particle.e = tlVector[3];
        outTree.Fill();
        partID++;
      }
    }
    delete hCSSliceAtW;
  }

  outFile.Write();
  outFile.Close();
}