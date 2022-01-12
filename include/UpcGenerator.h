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

// main class of the generator

#ifndef UPCGENERATOR__UPCGENERATOR_H_
#define UPCGENERATOR__UPCGENERATOR_H_

#include "UpcPythia8Helper.h"
#include "UpcPythia6Helper.h"

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

#ifdef USE_HEPMC
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#endif

using namespace std;

class UpcGenerator
{
 public:
  UpcGenerator();
  ~UpcGenerator();

  // parse inputs, set flags, prepare caches...
  void init();

  // debug level:
  //  0  -- no debug output
  //  >0 -- enable debug info
  void setDebugLevel(int level) { debug = level; }

  // number of threads for two-photon luminosity calculation
  void setNumThreads(int n) { numThreads = n; }

  // ----------------------------------------------------------------------

  // file parser
  void initGeneratorFromFile();

  // print parameters
  void printParameters();

  // the main method
  void generateEvents();

 private:
  // internal methods for event treating
  // ----------------------------------------------------------------------
  void processInPythia(vector<int>& pdgs,
                       vector<int>& statuses,
                       vector<int>& mothers,
                       vector<TLorentzVector>& particles);

  // helper struct for file output
  struct {
    int eventNumber;
    int pdgCode;
    int particleID;
    int statusID;
    int motherID;
    double px;
    double py;
    double pz;
    double e;
  } particle;

  TFile* mOutFile;
  TTree* mOutTree;

  void writeEvent(int evt,
                  const vector<int>& pdgs,
                  const vector<int>& statuses,
                  const vector<int>& mothers,
                  const vector<TLorentzVector>& particles);

  // pythia helper & decayer parameters
  int pythiaVersion{-1}; // not using Pythia at all by default
  bool isPythiaUsed{false};
#if defined(USE_PYTHIA6) || defined(USE_PYTHIA8)
  UpcPythiaBase* decayer;
#endif
  bool doFSR{false};
  bool doDecays{false};

#ifdef USE_HEPMC
  // helper for HepMC output format
  HepMC3::WriterAscii* writerHepMC;
#endif

  // number of worker threads for OpenMP
  int numThreads{1};

  // internal methods for calculations
  // ----------------------------------------------------------------------

  // Simpson integrator
  static double simpson(int n, double* v, double h);

  // Woods-Saxon rho0 from normalization
  double calcWSRho();

  // photon fluxes
  double fluxPoint(const double b, const double k);

  static double fluxFormInt(double* x, double* par);

  static double calcFormFac(double Q2);

  double fluxForm(const double b, const double k, TF1* fFluxForm);

  // two-photon luminosity
  double calcTwoPhotonLumi(double M, double Y, TF1* fFluxForm, const TGraph* gGAA);

  // polarized elementary cross sections
  double crossSectionMZPolS(double m, double z);

  double crossSectionMPolS(double m);

  double crossSectionMZPolPS(double m, double z);

  double crossSectionMPolPS(double m);

  // two-photon luminosity for scalar part
  void calcTwoPhotonLumiPol(double& ns, double& np, double M, double Y, TF1* fFluxForm, const TGraph* gGAA);

  // elementary cross section for dilepton production in MZ space
  double crossSectionMZ(double m, double z);

  // elementary cross section for dilepton production in M space
  double crossSectionM(double m);

  // histogram filler for MZ-cross section
  void fillCrossSectionMZ(TH2D* hCrossSectionMZ,
                          double mmin, double mmax, int nm,
                          double zmin, double zmax, int nz,
                          int flag);

  // histogram filler for M-cross section
  void fillCrossSectionM(TH1D* hCrossSectionM,
                         double mmin, double mmax, int nm);

  // function to calculate nuclear cross section
  // using 2D elementary cross section and two-photon luminosity
  void nuclearCrossSectionYM(TH2D* hCrossSectionYM, TH2D* hPolCSRatio);

  // nuclear form factor for momentum transfer q
  static double nucFormFactor(double t);

  // functions for calculating pair momentum
  // accounting for non-zero photon pt
  double getPhotonPt(double ePhot);
  void getPairMomentum(double mPair, double yPair, TLorentzVector& pPair);

  // various cachers-getters for lookup tables
  // ----------------------------------------------------------------------

  // prepare G_AA
  void prepareGAA();

  // prepare form factor
  void prepareFormFac();
  static double getCachedFormFac(double Q2);

  // prepare two photon luminosity, cache to file
  void prepareTwoPhotonLumi();

  // simulation & calculation parameters
  // ----------------------------------------------------------------------
  long seed{-1}; // seed for random numbers generator

  int lepPDG{15};       // tau by default
  double mLep{1.77682}; // tau by default
  double aLep{0};       // lepton anomalous magnetic moment

  // physics constants
  constexpr static double alpha{1.0 / 137.035999074}; // fine structure constant
  constexpr static double hc{0.1973269718};           // scaling factor
  constexpr static double mProt{0.9382720813};        // proton mass

  // Woods-Saxon parameters
  static double rho0; // fm-3
  static double R;    // fm
  static double a;    // fm

  // parameters of the nucleus
  static double Z;
  double A{208};

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

  // Gaussian integration n = 16
  static const int ngi16 = 16;
  double abscissas16[ngi16]{-0.9894009349916499,
                            -0.9445750230732326,
                            -0.8656312023878318,
                            -0.7554044083550030,
                            -0.6178762444026438,
                            -0.4580167776572274,
                            -0.2816035507792589,
                            -0.0950125098376374,
                            0.0950125098376374,
                            0.2816035507792589,
                            0.4580167776572274,
                            0.6178762444026438,
                            0.7554044083550030,
                            0.8656312023878318,
                            0.9445750230732326,
                            0.9894009349916499};

  double weights16[ngi16]{0.0271524594117541,
                          0.0622535239386479,
                          0.0951585116824928,
                          0.1246289712555339,
                          0.1495959888165767,
                          0.1691565193950025,
                          0.1826034150449236,
                          0.1894506104550685,
                          0.1894506104550685,
                          0.1826034150449236,
                          0.1691565193950025,
                          0.1495959888165767,
                          0.1246289712555339,
                          0.0951585116824928,
                          0.0622535239386479,
                          0.0271524594117541};

  // photon luminosity calculation parameters
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

  // scaling factor
  double factor{Z * Z * alpha / M_PI / M_PI / hc / hc};

  // helper containers for calculations
  static const int nb{200};
  double bmax{20};
  double db{bmax / (nb - 1)};
  double vb[nb];
  double vs[nb];
  double TA[nb];
  double rho[nb][nb];
  double vGAA[nb];
  double vRho[nb];

  // lookup tables
  static constexpr double Q2min{1e-9};
  static constexpr double Q2max{100};
  static const int nQ2{10000000};
  static constexpr double dQ2{(Q2max - Q2min) / nQ2};
  static double* vCachedFormFac; // Q^2-grid for possible form factor values

  // simulation parameters
  bool doPtCut{false};
  double minPt{0};
  int nEvents{1000};
  bool isPoint{true}; // flux calculation parameter
  bool useNonzeroGamPt{true};
  bool usePolarizedCS{false};
  static std::map<int, double> lepMassMap;

  // parameters dictionary
  // todo: use <any> from c++17 for a neat parsing???
  struct InputPars {
    string inNucZ{"NUCLEUS_Z"};
    string inNucA{"NUCLEUS_A"};
    string inWSRadius{"WS_R"};
    string inWSA{"WS_A"};
    string inCMSqrtS{"SQRTS"};
    string inLepPDG{"LEP_PDG"};
    string inLepA{"LEP_A"};
    string inDoPtCut{"DO_PT_CUT"};
    string inNEvents{"NEVENTS"};
    string inLowPt{"PT_MIN"};
    string inLowZ{"ZMIN"};
    string inHiZ{"ZMAX"};
    string inLowM{"MMIN"};
    string inHiM{"MMAX"};
    string inLowY{"YMIN"};
    string inHiY{"YMAX"};
    string inBinsZ{"BINS_Z"};
    string inBinsM{"BINS_M"};
    string inBinsY{"BINS_Y"};
    string inFluxPoint{"FLUX_POINT"};
    string inNonzeroGamPt{"NON_ZERO_GAM_PT"};
    string inPolarized{"USE_POLARIZED_CS"};
    string inPythiaVer{"PYTHIA_VERSION"};
    string inPythia8FSR{"PYTHIA8_FSR"};
    string inPythia8Decays{"PYTHIA8_DECAYS"};
    string inSeed{"SEED"};
  };

  // debug level
  static int debug;
};

#endif // UPCGENERATOR__UPCGENERATOR_H_
