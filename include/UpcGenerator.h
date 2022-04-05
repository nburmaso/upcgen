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

#ifdef USE_HEPMC
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#endif

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
  struct InputPars {
    std::string inNucZ{"NUCLEUS_Z"};
    std::string inNucA{"NUCLEUS_A"};
    std::string inWSRadius{"WS_R"};
    std::string inWSA{"WS_A"};
    std::string inCMSqrtS{"SQRTS"};
    std::string inProcID{"PROC_ID"};
    std::string inLepA{"LEP_A"};
    std::string inDoPtCut{"DO_PT_CUT"};
    std::string inLowPt{"PT_MIN"};
    std::string inDoEtaCut{"DO_ETA_CUT"};
    std::string inLowEta{"ETA_MIN"};
    std::string inHiEta{"ETA_MAX"};
    std::string inNEvents{"NEVENTS"};
    std::string inLowZ{"ZMIN"};
    std::string inHiZ{"ZMAX"};
    std::string inLowM{"MMIN"};
    std::string inHiM{"MMAX"};
    std::string inLowY{"YMIN"};
    std::string inHiY{"YMAX"};
    std::string inBinsZ{"BINS_Z"};
    std::string inBinsM{"BINS_M"};
    std::string inBinsY{"BINS_Y"};
    std::string inFluxPoint{"FLUX_POINT"};
    std::string inBreakupMode{"BREAKUP_MODE"};
    std::string inNonzeroGamPt{"NON_ZERO_GAM_PT"};
    std::string inPolarized{"USE_POLARIZED_CS"};
    std::string inPythiaVer{"PYTHIA_VERSION"};
    std::string inPythia8FSR{"PYTHIA8_FSR"};
    std::string inPythia8Decays{"PYTHIA8_DECAYS"};
    std::string inSeed{"SEED"};
    // fixme : to be removed
    std::string inDoMCut{"DO_M_CUT"};
    std::string inLowMCut{"LOW_M_CUT"};
    std::string inHighMCut{"HIGH_M_CUT"};
  };

  // debug level
  static int debug;

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
  UpcCrossSection* nucProcessCS{nullptr};

  // helper struct for file output
  struct outParticle {
    int eventNumber;
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

  void writeEvent(int evt,
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

#ifdef USE_HEPMC
  // helper for HepMC output format
  HepMC3::WriterAscii* writerHepMC;
#endif

  // number of worker threads for OpenMP
  int numThreads{1};

  // internal methods for event treating
  // ----------------------------------------------------------------------
  void processInPythia(std::vector<int>& pdgs,
                       std::vector<int>& statuses,
                       std::vector<int>& mothers,
                       std::vector<TLorentzVector>& particles);

  void processPions(std::vector<int>& pdgs,
                    std::vector<int>& statuses,
                    std::vector<int>& mothers,
                    std::vector<TLorentzVector>& particles);
};

#endif // UPCGENERATOR__UPCGENERATOR_H_
