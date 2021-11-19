//
// created by Nazar Burmasov on 6/25/21.
//

#ifndef UPCGENERATOR__UPCGENERATOR_H_
#define UPCGENERATOR__UPCGENERATOR_H_

#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TRandom.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "plog/Log.h"
#include "plog/Initializers/RollingFileInitializer.h"

#ifdef USE_PYTHIA6
#include "TPythia6.h"
#include "TPythia6Decayer.h"
#endif

#ifdef USE_PYTHIA8
#include "Pythia8/Pythia.h"
#include "TPythia8.h"
#include "TPythia8Decayer.h"
#endif

using namespace std;

class UpcGenerator
{
 public:
  UpcGenerator();
  ~UpcGenerator();

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
  int numThreads{1};

  // internal class methods
  // ----------------------------------------------------------------------

  // Simpson integrator
  static double simpson(int n, double* v, double h);

  // photon fluxes
  double fluxPoint(const double b, const double w, const double g);

  static double fluxFormInt(double* x, double* par);

  double fluxForm(const double b, const double w, const double g, TF1* fFluxForm);

  // two-photon luminosity
  double D2LDMDY(double M, double Y, TF1* fFluxForm, const TGraph* gGAA);

  // elementary cross section for dilepton production in WZ space
  double crossSectionWZ(double s, double z);

  // elementary cross section for dilepton production in W space
  double crossSectionW(double w);

  // histogram filler for WZ-cross section
  void fillCrossSectionWZ(TH2D* hCrossSectionWZ,
                          double wmin, double wmax, int nw,
                          double zmin, double zmax, int nz);

  // histogram filler for W-cross section
  void fillCrossSectionW(TH1D* hCrossSectionW,
                         double wmin, double wmax, int nw);

  // function to calculate nuclear cross section
  // using 2D elementary cross section and 2-gamma luminosity
  void nuclearCrossSectionYM(TH2D* hCrossSectionYM);

  // nuclear form factor for momentum transfer q
  static double nucFormFactor(double q);

  // functions for calculating pair momentum
  // accounting for non-zero photon pt
  double getPhotonPt(double ePhot);
  void getPairMomentum(double mPair, double yPair, TLorentzVector& pPair);

  // pythia helper & decayer parameters
  int pythiaVersion{-1}; // not using Pythia at all by default
  bool isPythiaUsed{false};
#if defined(USE_PYTHIA6) || defined(USE_PYTHIA8)
  TVirtualMCDecayer* decayer;
#endif

  // simulation & calculation parameters
  // ----------------------------------------------------------------------
  long seed{-1}; // seed for random numbers generator

  int lepPDG{15};       // tau by default
  double mLep{1.77682}; // tau by default
  double aLep{0};       // lepton anomalous magnetic moment

  // physics constants
  const double alpha{1.0 / 137.035999074};  // fine structure constant
  constexpr static double hc{0.1973269718}; // scaling factor
  const double mProt{0.9382720813};         // proton mass

  // Woods-Saxon parameters
  static double rho0; // fm-3
  static double R;    // fm
  static double a;    // fm

  // parameters of the nucleus
  static double Z;
  double A{208};

  // beam parameters
  double sqrts{5020};
  double g1{sqrts / (2. * mProt)};
  double g2{sqrts / (2. * mProt)};

  // Gaussian integration n = 10
  // since cos is symmetric around 0 we only need 5
  // of the points in the gaussian integration.
  static const int ngi = 5;
  double weights[ngi]{0.2955242247147529, 0.2692667193099963,
                      0.2190863625159820, 0.1494513491505806,
                      0.0666713443086881};
  double abscissas[ngi]{0.1488743389816312, 0.4333953941292472,
                        0.6794095682990244, 0.8650633666889845,
                        0.9739065285171717};

  // photon luminosity calculation parameters
  const int nb1{120};
  const int nb2{120};

  // cross sections binning
  double zmin{-1};   // min z
  double zmax{1};    // max z
  int nz{100};       // nbins in z = cos(theta)
  double wmin{3.56}; // min W in GeV
  double wmax{50.};  // max W in GeV
  int nw{1001};      // n bins in W
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

  // simulation parameters
  bool doPtCut{false};
  double minPt{0};
  int nEvents{1000};
  bool isPoint{true}; // flux calculation parameter
  bool useNonzeroGamPt{true};
  static std::map<int, double> lepMassMap;

  // helper struct for file output
  struct Particle {
    int eventNumber;
    int pdgCode;
    int particleID;
    int motherID;
    double px;
    double py;
    double pz;
    double e;
  };

  // parameters dictionary
  // todo: use <any> from c++17 for a neat parsing???
  struct InputPars {
    string inNucZ{"NUCLEUS_Z"};
    string inNucA{"NUCLEUS_A"};
    string inWSRho0{"WS_RHO0"};
    string inWSRadius{"WS_RAD"};
    string inWSA{"WS_A"};
    string inCMSqrtS{"SQRTS"};
    string inLepPDG{"LEP_PDG"};
    string inLepA{"LEP_A"};
    string inDoPtCut{"DO_PT_CUT"};
    string inNEvents{"NEVENTS"};
    string inLowPt{"PT_MIN"};
    string inLowZ{"ZMIN"};
    string inHiZ{"ZMAX"};
    string inLowW{"WMIN"};
    string inHiW{"WMAX"};
    string inLowY{"YMIN"};
    string inHiY{"YMAX"};
    string inBinsZ{"BINS_Z"};
    string inBinsW{"BINS_W"};
    string inBinsY{"BINS_Y"};
    string inFluxPoint{"FLUX_POINT"};
    string inNonzeroGamPt{"NON_ZERO_GAM_PT"};
    string inPythiaVer{"PYTHIA_VERSION"};
    string inSeed{"SEED"};
  };

  // debug level
  static int debug;
};

#endif // UPCGENERATOR__UPCGENERATOR_H_
