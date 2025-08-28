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

/// helper class for calculations

#pragma once

#include <cmath>
#include <fstream>
#include <map>
#include <random>
#include <unordered_map>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include "TClonesArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TParticle.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TRandomGen.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "plog/Appenders/ColorConsoleAppender.h"
#include "plog/Formatters/TxtFormatter.h"
#include "plog/Init.h"
#include "plog/Initializers/RollingFileInitializer.h"
#include "plog/Log.h"

#include "UpcElemProcess.h"
#include "UpcPhysConstants.h"
#include "UpcTwoPhotonALP.h"
#include "UpcTwoPhotonDilep.h"
#include "UpcTwoPhotonDipion.h"
#include "UpcTwoPhotonLbyL.h"
#include "UpcPhotoNuclearVM.h"

class UpcCrossSection
{

 public:
  UpcCrossSection();
  ~UpcCrossSection();

  // Woods-Saxon parameters
  inline static double rho0{0.}; // fm^-3
  inline static double R{6.68};  // fm
  inline static double a{0.447}; // fm

  // parameters of the nucleus
  inline static int Z{82};
  inline static int A{208};
  inline static double mNucl{(Z * phc::mProt + (A - Z) * phc::mNeut) / A};

  // beam parameters
  inline static double sqrts{5020.};
  inline static double g1{sqrts / (2. * phc::mProt)};
  inline static double g2{sqrts / (2. * phc::mProt)};

  // Gaussian integration n = 10
  // since cos is symmetric around 0 we only need 5
  // of the points in the gaussian integration.
  static const int ngi10 = 10;
  double weights10[ngi10]{0.0666713443086881,
                          0.1494513491505806,
                          0.2190863625159820,
                          0.2692667193099963,
                          0.2955242247147529,
                          0.2955242247147529,
                          0.2692667193099963,
                          0.2190863625159820,
                          0.1494513491505806,
                          0.0666713443086881};

  double abscissas10[ngi10]{-0.9739065285171717,
                            -0.8650633666889845,
                            -0.6794095682990244,
                            -0.4333953941292472,
                            -0.1488743389816312,
                            0.1488743389816312,
                            0.4333953941292472,
                            0.6794095682990244,
                            0.8650633666889845,
                            0.9739065285171717};

  // photon luminosity calculation parameters
  TString lumiFileDirectory{"."};
  void setLumiFileDirectory(TString directory) { lumiFileDirectory = directory; };
  const int nb1{120};
  const int nb2{120};

  // cross sections binning
  double zmin{-1};   // min z
  double zmax{1};    // max z
  int nz{100};       // nbins in z = cos(theta)
  double mmin{3.56}; // min M in GeV
  double mmax{50.};  // max M in GeV
  int nm{1001};      // n bins in M
  double ymin{-6.};  // min pair rapidity
  double ymax{6.};   // max pair rapidity
  int ny{121};       // n bins in Y

  bool doMassCut{false};
  double lowMCut{0};
  double hiMCut{9999};

  // scaling factor
  double factor{Z * Z * phc::alpha / M_PI / M_PI / phc::hc / phc::hc};

  // helper containers for calculations
  // ----------------------------------------------------------------------
  static const int nb{200};

  // 'constant' for photon pt calculation
  // initialized once in the constructor
  double gtot;

  // map with photon pt distributions as function of e_photon
  std::unordered_map<int, std::pair<int, TH1D>> photPtDistrMap;

  // simulation parameters
  bool isPoint{true}; // flux calculation parameter
  bool useNonzeroGamPt{true};
  bool usePolarizedCS{false};
  // 1 -- no breakup
  // 2 -- XNXN
  // 3 -- 0N0N
  // 4 -- 0NXN (+ XN0N)
  int breakupMode{1};
  int shadowingOption{0};
  int dghtPDG{};

  // debug level
  inline static int debug{0};

  // openmp: n threads
  int numThreads{1};

  // elementary process setter
  void setElemProcess(int procID);
  UpcElemProcess* elemProcess{nullptr};

  // process-specific parameters (with some default values)
  double alpMass{1.}; // 1 GeV
  double alpWidth{0.010}; // 10 MeV

  int64_t seed{0};

  void init();
  void setSeed(int64_t _seed) { seed = _seed; }

  // methods for cross section calculation
  // ----------------------------------------------------------------------

  // Simpson integrator
  template <typename ArrayType>
  static double simpson(int n, ArrayType* v, double h);

  // Woods-Saxon rho0 from normalization
  double calcWSRho();

  // photon fluxes
  double fluxPoint(double b, double k);

  static double calcFormFac(double Q2);

  double fluxForm(double b, double k);

  // two-photon luminosity
  double calcTwoPhotonLumi(double M, double Y);

  // two-photon luminosity for scalar part
  void calcTwoPhotonLumiPol(double& ns, double& np, double M, double Y);

  double calcPhotonFlux(double M, double Y);

  // histogram filler for MZ-cross section
  void fillCrossSectionZM(std::vector<std::vector<double>>& crossSectionZM,
                          double mmin, double mmax, int nm,
                          double zmin, double zmax, int nz,
                          int flag);

  // function to calculate nuclear cross section
  // using 2D elementary cross section and two-photon luminosity
  void calcNucCrossSectionYM(std::vector<std::vector<double>>& crossSectionYM,
                             std::vector<std::vector<double>>& polCSRatio,
                             double& totCS);

  double vegasNucCrossSectionYM();

  void calcNucCrossSectionY(std::vector<std::vector<double>>& crossSectionY,
                            std::vector<std::vector<double>>& csYRatio,
                            double& totCS);

  double calcBreakupProb(double b, int mode);

  // functions for calculating pair momentum
  // accounting for non-zero photon pt
  double getPhotonPt(double ePhot);
  void getPairMomentum(double mPair, double yPair, TLorentzVector& pPair);
  void getMomentumVM(double m, double y, int target, TLorentzVector& pPair);

  // various cachers-getters for lookup tables
  // ----------------------------------------------------------------------

  // G_AA helpers
  void prepareGAA();

  // nuclear breakup probability helpers
  void prepareBreakupProb();

  // form factor helpers and some workarounds
  void prepareFormFac();

  // prepare two photon luminosity, cache to file
  void prepareTwoPhotonLumi();
};
