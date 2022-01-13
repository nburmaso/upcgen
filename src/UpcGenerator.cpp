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

#include "UpcGenerator.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace std;

// out-of-line initialization for static members
double UpcGenerator::rho0 = 0;
double UpcGenerator::R = 6.68;
double UpcGenerator::a = 0.447;
double UpcGenerator::Z = 82;
double UpcGenerator::sqrts = 5020;
double UpcGenerator::g1 = sqrts / (2. * mProt);
double UpcGenerator::g2 = sqrts / (2. * mProt);
int UpcGenerator::debug = 0;
std::map<int, double> UpcGenerator::lepMassMap;
double* UpcGenerator::vCachedFormFac = new double[UpcGenerator::nQ2];

UpcGenerator::UpcGenerator()
{
  // populating lepton mass map
  // data from PDG
  lepMassMap[11] = 0.000510998946;
  lepMassMap[13] = 0.1056583745;
  lepMassMap[15] = 1.77686;
}

void UpcGenerator::init()
{
  // get parameters from file
  initGeneratorFromFile();

  mLep = lepMassMap[lepPDG];

  // update scaling factor
  factor = Z * Z * alpha / M_PI / M_PI / hc / hc;

  // calculate Woods-Saxon rho0 from R and a
  rho0 = calcWSRho();

  // initialize the MT64 random number generator
  gRandom = new TRandomMT64();
  gRandom->SetSeed(seed == -1 ? time(nullptr) : seed);

  // initialize pythia for decays
#ifdef USE_PYTHIA8
  PLOG_INFO << "Initializing Pythia-based decayer";
  if (pythiaVersion == 8) {
    isPythiaUsed = true;
    decayer = new UpcPythia8Helper();
    decayer->setFSR(doFSR);
    decayer->setDecays(doDecays);
    decayer->init();
  }
#endif
#ifdef USE_PYTHIA6
  if (pythiaVersion == 6) {
    isPythiaUsed = true;
    decayer = new UpcPythia6Helper();
  }
#endif

  if (!isPythiaUsed) {
    PLOG_WARNING << "Decays with Pythia are not used!";
  }

  PLOG_WARNING << "Check inputs:";
  printParameters();

  // prepare caches
  prepareGAA();     // G_AA and Fourier-transformed G_AA
  prepareFormFac(); // nuclear form factor
  if (!usePolarizedCS) {
    prepareTwoPhotonLumi(); // calculate two-photon luminosity and save into a file
  }
}

UpcGenerator::~UpcGenerator() = default;

void UpcGenerator::initGeneratorFromFile()
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
      if (parameter == parDict.inLowM) {
        mmin = stod(parValue);
      }
      if (parameter == parDict.inHiM) {
        mmax = stod(parValue);
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
      if (parameter == parDict.inBinsM) {
        nm = stoi(parValue);
      }
      if (parameter == parDict.inBinsY) {
        ny = stoi(parValue);
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
      if (parameter == parDict.inPythia8FSR) {
        doFSR = stoi(parValue);
      }
      if (parameter == parDict.inPythia8Decays) {
        doDecays = stoi(parValue);
      }
      if (parameter == parDict.inNonzeroGamPt) {
        useNonzeroGamPt = stoi(parValue);
      }
      if (parameter == parDict.inPolarized) {
        usePolarizedCS = stoi(parValue);
      }
      if (parameter == parDict.inSeed) {
        seed = stol(parValue);
      }
    }
    fInputs.close();
  } else {
    PLOG_WARNING << "Input file not found! Using default parameters...";
    printParameters();
  }
}

void UpcGenerator::printParameters()
{
  PLOG_WARNING << "OMP_NTHREADS " << numThreads;
  PLOG_WARNING << "NUCLEUS_Z " << Z;
  PLOG_WARNING << "NUCLEUS_A " << A;
  PLOG_WARNING << "WS_R " << R;
  PLOG_WARNING << "WS_A " << a;
  PLOG_WARNING << "(CALCULATED) WS_RHO0 " << rho0;
  PLOG_WARNING << "SQRTS " << sqrts;
  PLOG_WARNING << "LEP_PDG " << lepPDG;
  PLOG_WARNING << "LEP_A " << aLep;
  PLOG_WARNING << "NEVENTS " << nEvents;
  PLOG_WARNING << "DO_PT_CUT " << doPtCut;
  PLOG_WARNING << "PT_MIN " << minPt;
  PLOG_WARNING << "ZMIN " << zmin;
  PLOG_WARNING << "ZMAX " << zmax;
  PLOG_WARNING << "MMIN " << mmin;
  PLOG_WARNING << "MMAX " << mmax;
  PLOG_WARNING << "YMIN " << ymin;
  PLOG_WARNING << "YMAX " << ymax;
  PLOG_WARNING << "BINS_Z " << nz;
  PLOG_WARNING << "BINS_M " << nm;
  PLOG_WARNING << "BINS_Y " << ny;
  PLOG_WARNING << "FLUX_POINT " << isPoint;
  PLOG_WARNING << "NON_ZERO_GAM_PT " << useNonzeroGamPt;
  PLOG_WARNING << "PYTHIA_VERSION " << pythiaVersion;
  PLOG_WARNING << "PYTHIA8_FSR " << doFSR;
  PLOG_WARNING << "PYTHIA8_DECAYS " << doDecays;
  PLOG_WARNING << "SEED " << seed;
}

void UpcGenerator::processInPythia(vector<int>& pdgs,
                                   vector<int>& statuses,
                                   vector<int>& mothers,
                                   vector<TLorentzVector>& particles)
{
#if defined(USE_PYTHIA6) || defined(USE_PYTHIA8)
  TClonesArray processedParts("TParticle");
  TLorentzVector tlVector;
  decayer->process(pdgs, statuses, particles);
  decayer->import(&processedParts);
  // todo: not very optimal, to be cleaned up
  pdgs.clear();
  statuses.clear();
  mothers.clear();
  particles.clear();
  for (int ip = 0; ip < processedParts.GetEntriesFast(); ip++) {
    auto* part = (TParticle*)processedParts.At(ip);
    if (debug > 1) {
      PLOG_DEBUG << "Particle info:";
      part->Print();
    }
    int pdg = part->GetPdgCode();
    int status = part->GetStatusCode();
    int mother = part->GetFirstMother();
    part->Momentum(tlVector);
    pdgs.emplace_back(pdg);
    statuses.emplace_back(status);
    mothers.emplace_back(mother);
    particles.emplace_back(tlVector);
  }
#endif
}

void UpcGenerator::writeEvent(int evt,
                              const vector<int>& pdgs,
                              const vector<int>& statuses,
                              const vector<int>& mothers,
                              const vector<TLorentzVector>& particles)
{
#ifndef USE_HEPMC
  for (int i = 0; i < particles.size(); i++) {
    particle.eventNumber = evt;
    particle.pdgCode = pdgs[i];
    particle.particleID = i;
    particle.statusID = statuses[i];
    particle.motherID = i <= 1 ? -1 : mothers[i]; // first two particles = primary leptons
    particle.px = particles[i].Px();
    particle.py = particles[i].Py();
    particle.pz = particles[i].Pz();
    particle.e = particles[i].E();
    mOutTree->Fill();
  }
#endif

#ifdef USE_HEPMC
  // todo: add incoming photons to vertices
  HepMC3::GenEvent eventHepMC;
  eventHepMC.set_event_number(evt);
  HepMC3::GenVertexPtr vertex = std::make_shared<HepMC3::GenVertex>();
  for (int i = 0; i < particles.size(); i++) {
    HepMC3::FourVector p;
    p.setPx(particles[i].Px());
    p.setPy(particles[i].Py());
    p.setPz(particles[i].Pz());
    p.setE(particles[i].E());
    int status = i <= 1 ? 1 : 2;
    HepMC3::GenParticlePtr part = std::make_shared<HepMC3::GenParticle>(p, pdgs[i], status);
    vertex->add_particle_out(part);
  }
  // dummy in-particle for correct output
  vertex->add_particle_in(std::make_shared<HepMC3::GenParticle>());
  eventHepMC.add_vertex(vertex);
  writerHepMC->write_event(eventHepMC);
#endif
}

template <typename ArrayType>
double UpcGenerator::simpson(int n, ArrayType* v, double h)
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

double UpcGenerator::calcWSRho()
{
  for (int ib = 0; ib < nb; ib++) {
    double r = ib * db;
    vRho[ib] = r * r / (1 + exp((r - R) / a));
  }
  double wsRho0 = A / simpson(nb, vRho, db) / 4 / M_PI;
  return wsRho0;
}

double UpcGenerator::fluxPoint(const double b, const double k)
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

double UpcGenerator::fluxFormInt(double* x, double* par)
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

double UpcGenerator::fluxForm(const double b, const double k, TF1* fFluxFormInt)
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

double UpcGenerator::calcTwoPhotonLumi(double M, double Y, TF1* fFluxForm, const TGraph* gGAA)
{
  // double differential luminosity
  double D2LDMDYx = 0.;
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
      for (int k = 0; k < ngi16; k++) {
        if (abscissas10[k] < 0) {
          continue;
        }
        double b = TMath::Sqrt(b1 * b1 + b2 * b2 + 2. * b1 * b2 * TMath::Cos(M_PI * abscissas10[k]));
        sum_phi += (b < 20.) ? (weights10[k] * gGAA->Eval(b) * 2) : (weights10[k] * 2.);
      }
      sum_b2 += flux[j] * M_PI * sum_phi * b2 * (b2h - b2l);
    }
    sum += fluxForm(b1, k1, fFluxForm) * sum_b2 * b1 * (b1h - b1l);
  }
  D2LDMDYx = M_PI * M * sum;
  return D2LDMDYx;
}

double UpcGenerator::crossSectionMPolS(double m)
{
  double r = 2 * mLep / m;
  double r2 = r * r;
  double r4 = r2 * r2;
  if (r > 1) {
    return 0;
  }
  //  double cs_s = M_PI * r2 * alpha * alpha * hc * hc / mLep / mLep * (1 + 1.5 * r2) *
  //                ((1 - 0.5 * r2) * log(-1 + 2 / r * (1 + sqrt(1 - r2))) - sqrt(1 - r2));
  double cs_s = 4 * M_PI * alpha * alpha * hc * hc / m / m * ((1 + r * r - 3. / 4. * r * r * r * r) * 2 * log(1 / r + sqrt(1 / r / r - 1)) - (1 + 3. / 2. * r * r) * sqrt(1 - r * r));
  return cs_s;
  // fm^2 //
}

double UpcGenerator::crossSectionMZPolS(double m, double z)
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

double UpcGenerator::crossSectionMPolPS(double m)
{
  double r = 2 * mLep / m;
  double r2 = r * r;
  double r4 = r2 * r2;
  if (r > 1) {
    return 0;
  }
  //  double cs_ps = M_PI * r2 * alpha * alpha * hc * hc / mLep / mLep *
  //                 ((1 + r2 - 0.25 * r4) * log(-1 + 2 / r * (1 + sqrt(1 - r2))) - sqrt(1 - r2) * (1 + 0.5 * r2));
  //  return cs_ps;
  return 4 * M_PI * alpha * alpha * hc * hc / m / m * ((1 + r * r - 1. / 4. * r * r * r * r) * 2 * log(1 / r + sqrt(1 / r / r - 1)) - (1 + 1. / 2. * r * r) * sqrt(1 - r * r));
  // fm^2 //
}

double UpcGenerator::crossSectionMZPolPS(double m, double z)
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

void UpcGenerator::calcTwoPhotonLumiPol(double& ns, double& np, double M, double Y, TF1* fFluxForm, const TGraph* gGAA)
{
  double k1 = M / 2. * exp(Y);
  double k2 = M / 2. * exp(-Y);

  double xmin = isPoint ? 1 * R : 0.05 * R;
  double xmax = max(5. * g1 * hc / k1, 5 * R);
  double log_delta_x = (log(xmax) - log(xmin)) / nb1;

  double b1min = isPoint ? 1 * R : 0.05 * R;
  double b1max = max(5. * g1 * hc / k2, 5 * R);
  double log_delta_b = (log(b1max) - log(b1min)) / nb1;

  constexpr int ncached = 1000;
  vector<double> flux(ncached, 0);
  vector<double> crdGrid(ncached, 0);
  double bcmin = isPoint ? 1 * R : 1e-7;
  double bcmax = max(xmax, b1max);
  double log_delta_bc = (log(bcmax) - log(bcmin)) / ncached;
  for (int j = 0; j < ncached; j++) {
    double bl = bcmin * exp(j * log_delta_bc);
    double bh = bcmin * exp((j + 1) * log_delta_bc);
    double b = (bh + bl) / 2.;
    flux[j] = fluxForm(b, k1, fFluxForm);
    crdGrid[j] = b;
  }

  // simple linear fit for flux
  auto approxFF = [&](double x) {
    auto it_upp = lower_bound(crdGrid.begin(), crdGrid.end(), x);
    int upp = it_upp - crdGrid.begin();
    int low = upp - 1;
    double c = (flux[low] - flux[upp]) / (crdGrid[low] - crdGrid[upp]);
    double d = flux[low] - c * crdGrid[low];
    double fit = c * x + d;
    return fit;
  };

  double sum_b_s = 0; // scalar part
  double sum_b_p = 0; // pseudoscalar part
  for (int ib = 0; ib < nb1; ib++) {
    double bl = b1min * exp(ib * log_delta_b);
    double bh = b1min * exp((ib + 1) * log_delta_b);
    double b = (bh + bl) / 2.;
    double sum_x_s = 0;
    double sum_x_p = 0;
    for (int ix = 0; ix < nb1; ix++) {
      double xl = xmin * exp(ix * log_delta_x);
      double xh = xmin * exp((ix + 1) * log_delta_x);
      double x = (xh + xl) / 2.;
      double sum_phi_s = 0.;
      double sum_phi_p = 0.;
      double ff_x = fluxForm(x, k2, fFluxForm);
      for (int k = 0; k < ngi16; k++) {
        double phi = M_PI * (1 + abscissas10[k]);
        double cphi = cos(phi);
        double sphi = sin(phi);
        double xmb = sqrt(x * x + b * b - 2. * x * b * cphi);
        double ff_xmb = approxFF(xmb);
        double xs = (x - b * cphi) / xmb;
        double xp = b * sphi / xmb;
        double item = weights10[k] * ff_xmb * ff_x;
        sum_phi_s += item * xs * xs;
        sum_phi_p += item * xp * xp;
      }
      sum_x_s += sum_phi_s * x * (xh - xl);
      sum_x_p += sum_phi_p * x * (xh - xl);
    }
    double gaa = b < 20 ? gGAA->Eval(b) : 1; // |\vec{x_t - b} - \vec{x_t}| - |b|
    sum_b_s += sum_x_s * b * (bh - bl);
    sum_b_p += sum_x_p * b * (bh - bl);
  }
  ns = M_PI * M_PI * M * sum_b_s;
  np = M_PI * M_PI * M * sum_b_p;

  flux.clear();
  crdGrid.clear();
}

double UpcGenerator::crossSectionMZ(double m, double z)
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

double UpcGenerator::crossSectionM(double m)
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

void UpcGenerator::fillCrossSectionMZ(TH2D* hCrossSectionMZ,
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
        cs = crossSectionMZ(m, z);
      }
      if (flag == 1) { // scalar part
        cs = crossSectionMZPolS(m, z);
      }
      if (flag == 2) { // pseudoscalar part
        cs = crossSectionMZPolPS(m, z);
      }
      hCrossSectionMZ->SetBinContent(im, iz, cs);
    }
  }
  double scalingFactor = hc * hc * 1e7; // to [nb]
  hCrossSectionMZ->Scale(scalingFactor / dm);
}

void UpcGenerator::fillCrossSectionM(TH1D* hCrossSectionM,
                                     double mmin, double mmax, int nm)
{
  double m;
  double dm = (mmax - mmin) / nm;
  double cs;
  for (int im = 0; im <= nm; im++) {
    m = mmin + dm * im;
    cs = crossSectionM(m);
    hCrossSectionM->SetBinContent(im, cs);
  }
  double scalingFactor = hc * hc * 1e7; // to [nb]
  hCrossSectionM->Scale(scalingFactor);
}

void UpcGenerator::prepareGAA()
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

double UpcGenerator::calcFormFac(double Q2)
{
  double Q = sqrt(Q2) / hc;
  double coshVal = TMath::CosH(M_PI * Q * a);
  double sinhVal = TMath::SinH(M_PI * Q * a);
  double ff = 4 * M_PI * M_PI * rho0 * a * a * a / (Q * a * Q * a * sinhVal * sinhVal) *
              (M_PI * Q * a * coshVal * TMath::Sin(Q * R) - Q * R * TMath::Cos(Q * R) * sinhVal);
  ff += 8 * M_PI * rho0 * a * a * a * exp(-R / a) / (1 + Q * Q * a * a) / (1 + Q * Q * a * a);
  return ff;
}

void UpcGenerator::prepareFormFac()
{
  for (int iQ2 = 0; iQ2 < nQ2; iQ2++) {
    double Q2 = Q2min + iQ2 * dQ2;
    double ff = calcFormFac(Q2);
    vCachedFormFac[iQ2] = ff;
  }
}

double UpcGenerator::getCachedFormFac(double Q2)
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

void UpcGenerator::prepareTwoPhotonLumi()
{
  // calculating two-photon luminosity (if needed)
  // -----------------------------------------------------------------------
  if (!gSystem->AccessPathName("hD2LDMDY.root")) {
    PLOG_INFO << "Found precalculated 2D luminosity";
  } else {
    PLOG_INFO << "Precalculated 2D luminosity is not found. Starting all over...";
    double dy = (ymax - ymin) / (ny - 1);
    double dm = (mmax - mmin) / (nm - 1);
    // using ether parallel or serial implementation
#ifdef USE_OPENMP
    auto* f2DLumi = new TFile("hD2LDMDY.root", "recreate");
    auto* hD2LDMDY = new TH2D("hD2LDMDY", ";;", nm, mmin, mmax, ny, ymin, ymax);
    ROOT::EnableThreadSafety();
    int im, iy;
    int progress = 0;
    int total = ny * nm;
    omp_set_num_threads(numThreads);
#pragma omp parallel default(none)           \
  shared(hD2LDMDY, progress) private(im, iy) \
    firstprivate(nb, vb, total, numThreads, dm, dy, ymin, ymax, mmin, mmax, nm, ny, abscissas10, weights10, vGAA)
    {
      auto* fFluxFormInt = new TF1(Form("fFluxFormInt_private_%d", omp_get_thread_num()), fluxFormInt, 0, 10, 3);
      auto* gGAA = new TGraph(nb, vb, vGAA);
      auto* hD2LDMDY_private = new TH2D(Form("hD2LDMDY_private_%d", omp_get_thread_num()), ";;", nm, mmin, mmax, ny, ymin, ymax);
      int threadNum = omp_get_thread_num();
      int lowM = nm * threadNum / numThreads;
      int highM = nm * (threadNum + 1) / numThreads;
      for (im = lowM; im < highM; im++) {
        double m = mmin + dm * im;
        for (iy = 0; iy < ny; iy++) {
          double y = ymin + dy * iy;
          double item = calcTwoPhotonLumi(m, y, fFluxFormInt, gGAA);
          hD2LDMDY_private->SetBinContent(im, iy, item * dm * dy);
          progress++;
        }
        if (threadNum == 0) {
          double progressBar = 100. * progress / total;
          PLOG_INFO << "Calculating two-photon luminosity: " << fixed << setprecision(2) << progressBar << "%";
        }
      }
#pragma omp critical
      {
        hD2LDMDY->Add(hD2LDMDY_private);
      }
    }
    omp_set_num_threads(1);
#else
    auto* fFluxFormInt = new TF1("fFluxFormInt", fluxFormInt, 0, 10, 3);
    auto* gGAA = new TGraph(nb, vb, vGAA);
    auto* f2DLumi = new TFile("hD2LDMDY.root", "recreate");
    auto* hD2LDMDY = new TH2D("hD2LDMDY", ";;", nm, mmin, mmax, ny, ymin, ymax);
    PLOG_INFO << "Calculating 2D luminosity grid..." << 0 << "%";
    for (int im = 0; im < nm; im++) {
      double M = mmin + dm * im;
      for (int iy = 0; iy < ny; iy++) {
        double y = ymin + dy * iy;
        hD2LDMDY->SetBinContent(im, iy, D2LDMDY(M, y, fFluxFormInt, gGAA) * dm * dy);
      }
      double progressBar = 100. * (double)(im + 1) / nm;
      PLOG_INFO << "Calculating 2D luminosity grid..." << fixed << setprecision(2) << progressBar << "%";
    }
#endif
    hD2LDMDY->Write();
    f2DLumi->Close();

    PLOG_INFO << "Two-photon luminosity was written to 'hD2LDMDY.root'";
  }
}

void UpcGenerator::nuclearCrossSectionYM(TH2D* hCrossSectionYM, TH2D* hPolCSRatio)
{
  PLOG_INFO << "Calculating nuclear cross section for a_lep = " << aLep;

  double dy = (ymax - ymin) / (ny - 1);
  double dm = (mmax - mmin) / (nm - 1);

  TH2D* hD2LDMDY;

  if (usePolarizedCS) {
    if (!hPolCSRatio) {
      PLOG_FATAL << "hPolRatio is not initialized!";
      std::_Exit(-1);
    }
  } else { // loading pre-cached two-photon luminosity
    auto* f2DLumi = new TFile("hD2LDMDY.root", "r");
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
#pragma omp parallel default(none)                           \
  shared(cs, cs_rat, hD2LDMDY, progress) private(im, iy, ib) \
    firstprivate(nb, vb, total, numThreads, dm, dy, ymin, ymax, mmin, mmax, nm, ny, abscissas10, weights10, vGAA)
  {
    vector<vector<double>> cs_private(nm, vector<double>(ny, 0));
    vector<vector<double>> rat_private(nm, vector<double>(ny, 0));
    TH2D* hD2LDMDY_private;
    TF1* fFluxFormInt_private;
    TGraph* gGAA_private;
    if (usePolarizedCS) {
      fFluxFormInt_private = new TF1(Form("fFluxFormInt_private_%d", omp_get_thread_num()), fluxFormInt, 0, 10, 3);
      gGAA_private = new TGraph(nb, vb, vGAA);
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
          cs_private[im][iy] = crossSectionM(m) * hD2LDMDY_private->GetBinContent(im, iy);
        } else { // polarized
          double y = ymin + dy * iy;
          double cs_s = crossSectionMPolS(m);
          double cs_p = crossSectionMPolPS(m);
          double ns, np;
          calcTwoPhotonLumiPol(ns, np, m, y, fFluxFormInt_private, gGAA_private);
          double ncs_s = ns * cs_s; // scalar part
          double ncs_p = np * cs_p; // psudoscalar part
          double ncs = ncs_s + ncs_p;
          ncs *= dm * dy;
          cs_private[im][iy] = ncs * 1e7; // fm^2 -> nb
          rat_private[im][iy] = ncs_s / ncs_p;
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
double UpcGenerator::nucFormFactor(double t)
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
double UpcGenerator::getPhotonPt(double ePhot)
{
  constexpr double pi2x4 = 4 * M_PI * M_PI;
  double y1 = TMath::ACosH(g1);
  double y2 = -TMath::ACosH(g2);
  double gtot = TMath::CosH((y1 - y2) / 2.);

  double ereds = (ePhot / gtot) * (ePhot / gtot);
  double Cm = TMath::Sqrt(3.) * ePhot / gtot;
  double arg = Cm * Cm + ereds;
  double sFFactCM = nucFormFactor(arg);
  double Coef = 3. * (sFFactCM * sFFactCM * Cm * Cm * Cm) / (pi2x4 * arg * arg);

  double x = gRandom->Uniform(0, 1);
  double pp = x * 5. * hc / R;
  arg = pp * pp + ereds;
  double sFFactPt1 = nucFormFactor(arg);
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
      double sFFactPt2 = nucFormFactor(arg);
      test = (sFFactPt2 * sFFactPt2) * pp * pp * pp / (pi2x4 * arg * arg);
    }
  }

  return pp;
}

void UpcGenerator::getPairMomentum(double mPair, double yPair, TLorentzVector& pPair)
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

void UpcGenerator::generateEvents()
{
  // calculating elementary cross section in WZ space
  // -----------------------------------------------------------------------
  PLOG_INFO << "Calculating elementary cross section...";

  TH2D* hCrossSectionMZ;
  TH2D* hCrossSectionMZPolS;
  TH2D* hCrossSectionMZPolPS;

  if (usePolarizedCS) {
    hCrossSectionMZPolS = new TH2D("hCrossSectionMZPolS", ";m [gev]; z; cs [nb/gev]",
                                   nm, mmin, mmax,
                                   nz, zmin, zmax);

    fillCrossSectionMZ(hCrossSectionMZPolS, mmin, mmax, nm, zmin, zmax, nz, 1);

    hCrossSectionMZPolPS = new TH2D("hCrossSectionMZPolPs", ";m [gev]; z; cs [nb/gev]",
                                    nm, mmin, mmax,
                                    nz, zmin, zmax);

    fillCrossSectionMZ(hCrossSectionMZPolPS, mmin, mmax, nm, zmin, zmax, nz, 2);
  } else {
    hCrossSectionMZ = new TH2D("hCrossSectionMZ", ";m [gev]; z; cs [nb/gev]",
                               nm, mmin, mmax,
                               nz, zmin, zmax);

    fillCrossSectionMZ(hCrossSectionMZ, mmin, mmax, nm, zmin, zmax, nz, 0);
  }

  // calculating nuclear cross section in YM space
  // -----------------------------------------------------------------------
  auto* hNucCSYM = new TH2D("hNucCSYM", "", ny - 1, ymin, ymax, nm - 1, mmin, mmax);
  TH2D* hPolCSRatio;
  if (usePolarizedCS) {
    hPolCSRatio = new TH2D("hPolCSRatio", "", ny - 1, ymin, ymax, nm - 1, mmin, mmax);
  }
  nuclearCrossSectionYM(hNucCSYM, hPolCSRatio);

  if (debug > 0) {
    PLOG_DEBUG << "a_lep = " << aLep << ", min. pt = " << minPt;
  }

  vector<double> cutsZ(hNucCSYM->GetNbinsY());
  if (doPtCut) {
    for (int mBin = 1; mBin <= hNucCSYM->GetNbinsY(); mBin++) {
      double mass = hNucCSYM->GetYaxis()->GetBinLowEdge(mBin);
      double sqt = TMath::Sqrt(mass * mass * 0.25 - mLep * mLep);
      if (sqt <= minPt) {
        cutsZ[mBin - 1] = 0;
        for (int yBin = 1; yBin <= hNucCSYM->GetNbinsX(); yBin++) {
          hNucCSYM->SetBinContent(yBin, mBin, 0);
        }
        continue;
      }
      double minZ = TMath::Sqrt(1 - minPt * minPt / (sqt * sqt));
      cutsZ[mBin - 1] = minZ;
    }
  }

  vector<int> pdgs;
  vector<int> statuses;
  vector<int> mothers;
  vector<TLorentzVector> particles;

#ifndef USE_HEPMC
  // initialize file output
  PLOG_WARNING << "Using ROOT tree for output!";
  PLOG_INFO << "Events will be written to " << Form("events_%.3f_%.0f.root", aLep, minPt);
  mOutFile = new TFile(Form("events_%.3f_%.0f.root", aLep, minPt), "recreate", "", 4 * 100 + 5); // using LZ4 with level 5 compression
  mOutTree = new TTree("particles", "Generated particles");
  mOutTree->Branch("eventNumber", &particle.eventNumber, "eventNumber/I");
  mOutTree->Branch("pdgCode", &particle.pdgCode, "pdgCode/I");
  mOutTree->Branch("particleID", &particle.particleID, "particleID/I");
  mOutTree->Branch("statusID", &particle.statusID, "statusID/I");
  mOutTree->Branch("motherID", &particle.motherID, "motherID/I");
  mOutTree->Branch("px", &particle.px, "px/D");
  mOutTree->Branch("py", &particle.py, "py/D");
  mOutTree->Branch("pz", &particle.pz, "pz/D");
  mOutTree->Branch("e", &particle.e, "e/D");
  mOutTree->SetAutoSave(0);
#else
  PLOG_WARNING << "Using HepMC format for output!";
  PLOG_INFO << "Events will be written to " << Form("events_%.3f_%.0f.hepmc", aLep, minPt);
  writerHepMC = new HepMC3::WriterAscii(Form("events_%.3f_%.0f.hepmc", aLep, minPt));
#endif

  PLOG_INFO << "Generating " << nEvents << " events...";

  for (int evt = 0; evt < nEvents; evt++) {
    if (debug > 1) {
      PLOG_DEBUG << "Event number: " << evt + 1;
    } else if ((evt + 1) % 10000 == 0) {
      PLOG_INFO << "Event number: " << evt + 1;
    }

    double yPair;
    double mPair;

    hNucCSYM->GetRandom2(yPair, mPair);

    TH1D* hCSSliceAtM;
    if (!usePolarizedCS) {
      int binM = hCrossSectionMZ->GetXaxis()->FindBin(mPair);
      hCSSliceAtM = hCrossSectionMZ->ProjectionY("sliceM", binM, binM);
    } else {
      int frBinY = hPolCSRatio->GetXaxis()->FindBin(yPair);
      int frBinM = hPolCSRatio->GetYaxis()->FindBin(mPair);
      double frac = hPolCSRatio->GetBinContent(frBinM, frBinY);
      bool pickScalar = gRandom->Uniform(0, 1) < frac;
      if (pickScalar) {
        int binM = hCrossSectionMZPolS->GetXaxis()->FindBin(mPair);
        hCSSliceAtM = hCrossSectionMZPolS->ProjectionY("sliceM", binM, binM);
      } else {
        int binM = hCrossSectionMZPolPS->GetXaxis()->FindBin(mPair);
        hCSSliceAtM = hCrossSectionMZPolPS->ProjectionY("sliceM", binM, binM);
      }
    }

    if (doPtCut) {
      int mBin = hNucCSYM->GetYaxis()->FindBin(mPair);
      double zCut = cutsZ[mBin - 1];
      int zCutBinUp = hCSSliceAtM->GetXaxis()->FindBin(zCut);
      int zCutBinLow = hCSSliceAtM->GetXaxis()->FindBin(-zCut);
      for (int zBin = zCutBinUp + 1; zBin <= hCSSliceAtM->GetNbinsX(); zBin++) {
        hCSSliceAtM->SetBinContent(zBin, 0);
      }
      for (int zBin = 1; zBin < zCutBinLow; zBin++) {
        hCSSliceAtM->SetBinContent(zBin, 0);
      }
    }

    double cost = hCSSliceAtM->GetRandom();
    double theta = TMath::ACos(cost);
    double phi = gRandom->Uniform(0., 2. * M_PI);

    TLorentzVector pPair;
    getPairMomentum(mPair, yPair, pPair);
    double pMag = TMath::Sqrt(pPair.Mag() * pPair.Mag() / 4 - mLep * mLep);
    TVector3 vLep;
    vLep.SetMagThetaPhi(pMag, theta, phi);
    TLorentzVector tlLep1, tlLep2;
    tlLep1.SetVectM(vLep, mLep);
    tlLep2.SetVectM(-vLep, mLep);
    particles.emplace_back(tlLep1);
    particles.emplace_back(tlLep2);

    TVector3 boost = pPair.BoostVector();
    particles[0].Boost(boost);
    particles[1].Boost(boost);

    int sign = gRandom->Uniform(-1, 1) > 0 ? 1 : -1;
    pdgs.emplace_back(sign * lepPDG);
    pdgs.emplace_back(-sign * lepPDG);
    mothers.emplace_back(-1);
    mothers.emplace_back(-1);
    statuses.emplace_back(23);
    statuses.emplace_back(23);

    // lepton decays for taus
    // todo: at the moment "fsr" and "decays" flags
    //       are only really meaningful for pythia8
    if ((doFSR || doDecays) && isPythiaUsed) {
      processInPythia(pdgs, statuses, mothers, particles);
    }

    writeEvent(evt, pdgs, statuses, mothers, particles);

    pdgs.clear();
    statuses.clear();
    mothers.clear();
    particles.clear();
    delete hCSSliceAtM;
  }

#ifndef USE_HEPMC
  hNucCSYM->Write();
  mOutFile->Write();
  mOutFile->Close();
#endif

#ifdef USE_HEPMC
  writerHepMC->close();
#endif
}