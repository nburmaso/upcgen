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

// main class of the generator

#ifndef UPCGENERATOR__UPCGENERATOR_H_
#define UPCGENERATOR__UPCGENERATOR_H_

#include "UpcCrossSection.h"
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

class UpcGenerator
{
 public:
  UpcGenerator();
  ~UpcGenerator();

  // process-specific parameters
  int procID;
  double aLep{0};

  // simulation parameters
  bool doPtCut{false};
  double minPt{0};
  bool doEtaCut{false};
  double minEta{0};
  double maxEta{0};
  bool usePolarizedCS{false};
  long int seed{0};
  long int nEvents{1000};

  // parameters dictionary
  // todo: use <any> from c++17 for a neat parsing???
  static const int nParameters = 37;
  const std::string parsNames[nParameters]
  {
    // parameters of incoming nuclei
    "NUCLEUS_Z",
    "NUCLEUS_A",
    "WS_R",
    "WS_A",
    "SQRTS",
    "PROC_ID",
    "NEVENTS",
    // kinematic cuts for final-state particles (decay products not included)
    "DO_PT_CUT",
    "PT_MIN",
    "DO_ETA_CUT",
    "ETA_MIN",
    "ETA_MAX",
    // grid sizes and binnings
    "ZMIN",
    "ZMAX",
    "MMIN",
    "MMAX",
    "YMIN",
    "YMAX",
    "BINS_Z",
    "BINS_M",
    "BINS_Y",
    // switches and flags
    "FLUX_POINT",
    "BREAKUP_MODE",
    "NON_ZERO_GAM_PT",
    "USE_POLARIZED_CS",
    "PYTHIA_VERSION",
    "PYTHIA8_FSR",
    "PYTHIA8_DECAYS",
    "SEED",
    "USE_ROOT_OUTPUT",
    "USE_HEPMC_OUTPUT",
    // mass cuts: hack for pre-calculated cross sections (lbyl and pi0 pair production)
    // fixme : to be removed
    "DO_M_CUT",
    "LOW_M_CUT",
    "HIGH_M_CUT",
    // process-specific parameters
    "LEP_A",
    "ALP_MASS",
    "ALP_WIDTH"
  };

  void setParameterValue(std::string parameter, std::string parValue);

  // debug level
  static int debug;

  // set configuration parameters
  void configGeneratorFromFile();

  // perform initial calculations
  void init();

  // set the directory where the lumi file is expected to be
  void setLumiFileDirectory(std::string directory) { nucProcessCS->setLumiFileDirectory(directory); };

  // debug level:
  //  0  -- no debug output
  //  >0 -- enable debug info
  void setDebugLevel(int level) { debug = level; }

  // number of threads for two-photon luminosity calculation
  void setNumThreads(int n) { numThreads = n; }
  
  void setParFile(std::string parfilename) { parFileName = parfilename; }

  // ----------------------------------------------------------------------
  // set beam parameters
  void setCollisionSystem(float sqrts, int nucl_z, int nucl_a);

  // set seed for random generator
  void setSeed(int seedIn) {
    seed = seedIn;
    
    // initialize the MT64 random number generator
    delete gRandom;
    gRandom = new TRandomMT64();

    PLOG_INFO << "<seed> = " << seed;
    gRandom->SetSeed(seed == 0 ? time(nullptr) : seed);
  };

  // print parameters
  void printParameters();

  // compute nuclear cross section
  void computeNuclXsection();
  double totNuclX() { return totCS; }
  double fidNuclX() { return fidCS; }

  // event generation
  long int generateEvent(std::vector<int>& pdgs, std::vector<int>& statuses, std::vector<int>& mothers, std::vector<TLorentzVector>& particles);
  const std::vector<TParticle>& getParticles() const { return genParticles; };

  // the main method in an event loop
  void generateEvents();

 private:
  UpcCrossSection* nucProcessCS{nullptr};

  std::string parFileName{"parameters.in"};

  bool useROOTOut{true};
  bool useHepMCOut{false};

  // global variables needed for event generation
  double totCS;
  double fidCS;
  bool ignoreCSZ;
  TH2D* hNucCSYM;
  std::vector<std::vector<double>> hPolCSRatio;
  std::vector<TH1D*> hCrossSecsZ;
  std::vector<TH1D*> hCrossSecsZ_S;
  std::vector<TH1D*> hCrossSecsZ_PS;

  // helper struct for ROOT file output
  struct outParticle {
    long int eventNumber;
    int pdgCode;
    int particleID;
    int statusID;
    int motherID;
    double px;
    double py;
    double pz;
    double e;
  };

  TFile* mOutFile{nullptr};
  TTree* mOutTree{nullptr};
  outParticle particle{};
  std::vector<TParticle> genParticles;

  void writeEvent(long int evt,
                  const std::vector<int>& pdgs,
                  const std::vector<int>& statuses,
                  const std::vector<int>& mothers,
                  const std::vector<TLorentzVector>& particles);

  // pythia helper & decayer parameters
  int pythiaVersion{-1}; // not using Pythia at all by default
  bool isPythiaUsed{false};
#if defined(USE_PYTHIA6) || defined(USE_PYTHIA8)
  UpcPythiaBase* decayer{nullptr};
#endif
  bool doFSR{false};
  bool doDecays{false};

  // a very simple HepMC writer
  // writing only particles, omitting any vertex information
  class WriterHepMC
  {
   public:
    WriterHepMC(const std::string& fname)
    {
      openFile(fname);
    }

    ~WriterHepMC()
    {
      closeFile();
    }

    std::ofstream outfile;

    // open output file and print preamble
    void openFile(const std::string& fname)
    {
      outfile.open(fname);
      outfile << "HepMC::Version 3.02.04"
              << "\n"
              << "HepMC::Asciiv3-START_EVENT_LISTING"
              << "\n";
    }

    // write end-of-listing message and close output file
    void closeFile()
    {
      outfile << "HepMC::Asciiv3-END_EVENT_LISTING"
              << "\n";
      outfile.close();
    }

    // writing basic event info with default HepMC units
    void writeEventInfo(long int eventID, int nParticles, int nVertices = 0)
    {
      outfile << "E " << eventID << " " << nVertices << " " << nParticles << "\n"
              << "U GEV MM"
              << "\n";
    }

    void writeParticleInfo(int id, int motherID, int pdg,
                           double px, double py, double pz, double e, double m,
                           int status)
    {
      outfile << std::setprecision(9)
              << "P " << id << " " << motherID << " " << pdg << " "
              << px << " " << py << " " << pz << " " << e << " " << m << " "
              << status << "\n";
    }
  };

  WriterHepMC* writerHepMC;

  // number of worker threads for OpenMP
  int numThreads{1};

  // internal methods for event treating
  // ----------------------------------------------------------------------
  bool isPairProduction{false};
  void pairProduction(TLorentzVector& pPair,                  // input lorentz pair-momentum vector
                      TVector3& vec,                          // momentum of outgoing particles
                      double mPart,                           // outgoing particle mass
                      int partPDG,                            // pdg of final-state particles
                      bool isCharged,                         // final-state particles are charged or not
                      std::vector<TLorentzVector>& particles, // vector for final-state particles
                      std::vector<int>& pdgs,
                      std::vector<int>& mothers,
                      std::vector<int>& statuses);

  bool isSingleProduction{false};
  void singleProduction(TLorentzVector& pPair,                  // input lorentz pair-momentum vector
                        int partPDG,                            // pdg of final-state particle
                        std::vector<TLorentzVector>& particles, // vector for final-state particles
                        std::vector<int>& pdgs,
                        std::vector<int>& mothers,
                        std::vector<int>& statuses);

  void processInPythia(std::vector<int>& pdgs,
                       std::vector<int>& statuses,
                       std::vector<int>& mothers,
                       std::vector<std::vector<int>>& daughters,
                       std::vector<TLorentzVector>& particles);

  void twoPartDecayUniform(std::vector<int>& pdgs,
                           std::vector<int>& statuses,
                           std::vector<int>& mothers,
                           std::vector<TLorentzVector>& particles,
                           int id,
                           double decayMass,
                           int decayProdPDG);

  // check kinematic cuts for particles in 'particles' vector
  bool checkKinCuts(std::vector<TLorentzVector>& particles);
};

#endif // UPCGENERATOR__UPCGENERATOR_H_
