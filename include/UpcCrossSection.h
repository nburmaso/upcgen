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

/// helper class for calculations

#ifndef UPCGENERATOR_INCLUDE_UPCCROSSSECTION_H_
#define UPCGENERATOR_INCLUDE_UPCCROSSSECTION_H_

#include "UpcPhysConstants.h"
#include "UpcElemProcess.h"
#include "UpcTwoPhotonALP.h"
#include "UpcTwoPhotonDilep.h"
#include "UpcTwoPhotonDipion.h"
#include "UpcTwoPhotonLbyL.h"

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

#include <fstream>
#include <map>
#include <random>
#include <unordered_map>

class UpcCrossSection
{

 public:
  UpcCrossSection();
  ~UpcCrossSection();

  // Woods-Saxon parameters
  static double rho0; // fm-3
  static double R;    // fm
  static double a;    // fm

  // parameters of the nucleus
  static int Z;
  int A{208};

  // beam parameters
  static double sqrts;
  static double g1;
  static double g2;

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
  double factor{Z * Z * phys_consts::alpha / M_PI / M_PI / phys_consts::hc / phys_consts::hc};

  // helper containers for calculations
  // ----------------------------------------------------------------------
  static const int nb{200};
  double vb[nb]{};
  double vs[nb]{};
  double TA[nb]{};
  double rho[nb][nb]{};
  double vGAA[nb]{};
  double vRho[nb]{};

  // 'constant' for photon pt calculation
  // initialized once in the constructor
  double gtot;

  // lookup tables and maps for cached functions etc.
  static constexpr double Q2min{1e-9};
  static constexpr double Q2max{100};
  static constexpr int nQ2{10000000};
  static constexpr double dQ2{(Q2max - Q2min) / nQ2};
  static double* vCachedFormFac;                                // Q^2-grid for possible form factor values
  double* vCachedBreakup{nullptr};                              // b-grid for possible breakup probabilities
  std::unordered_map<int, std::pair<int, TH1D>> photPtDistrMap; // map with photon pt distributions in dependence of e_photon
  double avgFreq{-1};                                           // calculate average frequency only once

  // simulation parameters
  bool isPoint{true}; // flux calculation parameter
  bool useNonzeroGamPt{true};
  bool usePolarizedCS{false};
  // 1 -- no breakup
  // 2 -- XNXN
  // 3 -- 0N0N
  // 4 -- 0NXN (+ XN0N)
  int breakupMode{1};

  // debug level
  static int debug;

  // openmp: n threads
  int numThreads{1};

  // elementary process setter
  void setElemProcess(int procID);
  UpcElemProcess* elemProcess{nullptr};

  // process-specific parameters (with some default values)
  double alpMass{1.}; // 1 GeV
  double alpWidth{0.010}; // 10 MeV

  void init();

  // methods for cross section calculation
  // ----------------------------------------------------------------------

  // Simpson integrator
  template <typename ArrayType>
  static double simpson(int n, ArrayType* v, double h);

  // Woods-Saxon rho0 from normalization
  double calcWSRho();

  // photon fluxes
  double fluxPoint(double b, double k);

  static double fluxFormInt(double* x, double* par);

  static double calcFormFac(double Q2);

  double fluxForm(double b, double k, TF1* fFluxForm);

  // two-photon luminosity
  double calcTwoPhotonLumi(double M, double Y, TF1* fFluxForm, const TGraph* gGAA);

  // two-photon luminosity for scalar part
  void calcTwoPhotonLumiPol(double& ns, double& np, double M, double Y, TF1* fFluxForm, const TGraph* gGAA);

  // histogram filler for MZ-cross section
  void fillCrossSectionZM(TH2D* hCrossSectionZM,
                          double mmin, double mmax, int nm,
                          double zmin, double zmax, int nz,
                          int flag);

  // function to calculate nuclear cross section
  // using 2D elementary cross section and two-photon luminosity
  void calcNucCrossSectionYM(TH2D* hCrossSectionYM, std::vector<std::vector<double>>& hPolCSRatio, double& totCS);

  double calcBreakupProb(double b, int mode);

  // functions for calculating pair momentum
  // accounting for non-zero photon pt
  double getPhotonPt(double ePhot);
  void getPairMomentum(double mPair, double yPair, TLorentzVector& pPair);

  // various cachers-getters for lookup tables
  // ----------------------------------------------------------------------

  // prepare G_AA
  void prepareGAA();

  // prepare nuclear breakup probability
  void prepareBreakupProb();
  double getCachedBreakupProb(double b);

  // prepare form factor
  void prepareFormFac();
  static double getCachedFormFac(double Q2);

  // prepare two photon luminosity, cache to file
  void prepareTwoPhotonLumi();
};

#endif // UPCGENERATOR_INCLUDE_UPCCROSSSECTION_H_
