//////////////////////////////////////////////////////////////////////////
// Copyright (C) 2021-2022, Nazar Burmasov, Evgeny Kryshen
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
#include "UpcPhysConstants.h"
#include "UpcSampler.h"

#include "UpcTwoPhotonDipion.h"

using namespace std;

// out-of-line initialization for static members
int UpcGenerator::debug = 0;

UpcGenerator::UpcGenerator()
{
  nucProcessCS = new UpcCrossSection();
}

void UpcGenerator::init()
{
  if (nucProcessCS == nullptr) {
    PLOG_FATAL << "UpcCrossSection was not initialized! Exiting...";
    std::_Exit(-1);
  }

  TH1::AddDirectory(false);

  // get parameters from file
  initGeneratorFromFile();

  // sanity checks: lbyl
  if (procID == 22) {
    isPairProduction = true;
    // for gamma+gamma -> gamma+gamma grid size is set explicitly
    // because elem. cross sections are already stored in files
    PLOG_WARNING << "For light-by-light scattering grid sizes along Z and M are fixed -- see parameters check";
    nucProcessCS->zmin = -0.99;
    nucProcessCS->zmax = 0.99;
    nucProcessCS->nz = 198;
    nucProcessCS->mmin = 0.05;
    nucProcessCS->mmax = 50.;
    nucProcessCS->nm = 1000;
    if (usePolarizedCS) {
      PLOG_WARNING << "For light-by-light scattering polarized cross section is not available";
      usePolarizedCS = false;
    }
  }

  // sanity checks: pi0pi0
  if (procID == 111) {
    isPairProduction = true;
    // for gamma+gamma -> pi0 pi0 grid size is set explicitly
    // because elem. cross sections are already stored in files
    PLOG_WARNING << "For dipion photoproduction grid sizes along Z and M are fixed -- see parameters check";
    nucProcessCS->zmin = -1;
    nucProcessCS->zmax = 1;
    nucProcessCS->nz = 100;
    nucProcessCS->mmin = 0.275;
    nucProcessCS->mmax = 5.;
    nucProcessCS->nm = 91;
    if (usePolarizedCS) {
      PLOG_WARNING << "For dipion photoproduction polarized cross section is not available";
      usePolarizedCS = false;
    }
  }

  // checks for ALP production
  if (procID == 51) {
    isSingleProduction = true;
    double mass = nucProcessCS->alpMass;
    double width = nucProcessCS->alpWidth;
    // gauss peak, +/- 4 sigmas
    nucProcessCS->mmin = mass - 4. * width;
    nucProcessCS->mmax = mass + 4. * width;
    PLOG_WARNING << "For ALP production angular distribution is ignored!";
    PLOG_WARNING << "ALP will be generated according to chosen ALP_MASS and ALP_WIDTH: check input parameters print-out";
  }

  nucProcessCS->setElemProcess(procID);

  // checks for dilepton production
  if (procID >= 11 && procID <= 15) {
    isPairProduction = true;
    auto* proc = (UpcTwoPhotonDilep*)nucProcessCS->elemProcess;
    proc->aLep = aLep;
    double mPart = proc->mPart;
    if (nucProcessCS->mmin < mPart * 2.) {
      PLOG_WARNING << "MMIN is lower than 2 lepton masses! Setting MMIN to 2 lepton masses...";
      nucProcessCS->mmin = mPart * 2.;
    }
  }

  nucProcessCS->numThreads = numThreads;

  PLOG_WARNING << "Check inputs:";
  printParameters();

  nucProcessCS->init();

  // initialize the MT64 random number generator
  delete gRandom;
  gRandom = new TRandomMT64();
  gRandom->SetSeed(seed == 0 ? time(nullptr) : seed);

  // initialize pythia for decays
#ifdef USE_PYTHIA8
  PLOG_INFO << "Initializing Pythia-based decayer";
  if (pythiaVersion == 8) {
    isPythiaUsed = true;
    decayer = new UpcPythia8Helper();
    decayer->setFSR(doFSR);
    decayer->setDecays(doDecays);
    decayer->setSeed(seed);
    decayer->init();
  }
#endif
#ifdef USE_PYTHIA6
  if (pythiaVersion == 6) {
    isPythiaUsed = true;
    decayer = new UpcPythia6Helper();
    decayer->setSeed(seed);
  }
#endif

  if (!isPythiaUsed) {
    PLOG_WARNING << "Decays with Pythia are not used!";
  }
}

UpcGenerator::~UpcGenerator()
{
#if defined(USE_PYTHIA6) || defined(USE_PYTHIA8)
  delete decayer;
#endif
  delete nucProcessCS;
  delete mOutFile;
}

void UpcGenerator::initGeneratorFromFile()
{
  // todo: use <any> from c++17 for a neat parsing?
  if (!gSystem->AccessPathName("parameters.in")) {
    PLOG_INFO << "Reading parameters from parameters.in ...";
    InputPars parDict;
    ifstream fInputs("parameters.in");
    string line;
    string parameter;
    string parValue;
    while (getline(fInputs, line)) {
      istringstream iss(line);
      // skip comment lines
      if (line[0] == '#') {
        continue;
      }
      line.erase(std::find(line.begin(), line.end(), '#'), line.end()); // trim comments
      iss >> parameter >> parValue;
      if (parameter == parDict.inNEvents) {
        nEvents = stol(parValue);
      }
      if (parameter == parDict.inCMSqrtS) {
        UpcCrossSection::sqrts = stod(parValue);
        UpcCrossSection::g1 = UpcCrossSection::sqrts / (2. * phys_consts::mProt);
        UpcCrossSection::g2 = UpcCrossSection::sqrts / (2. * phys_consts::mProt);
      }
      if (parameter == parDict.inProcID) {
        procID = stoi(parValue);
      }
      if (parameter == parDict.inLepA) {
        aLep = stod(parValue);
      }
      if (parameter == parDict.inALPMass) {
        nucProcessCS->alpMass = stod(parValue);
      }
      if (parameter == parDict.inALPWidth) {
        nucProcessCS->alpWidth = stod(parValue);
      }
      if (parameter == parDict.inDoPtCut) {
        doPtCut = stoi(parValue);
      }
      if (parameter == parDict.inLowPt) {
        minPt = stod(parValue);
      }
      if (parameter == parDict.inDoEtaCut) {
        doEtaCut = stoi(parValue);
      }
      if (parameter == parDict.inLowEta) {
        minEta = stod(parValue);
      }
      if (parameter == parDict.inHiEta) {
        maxEta = stod(parValue);
      }
      if (parameter == parDict.inLowZ) {
        nucProcessCS->zmin = stod(parValue);
      }
      if (parameter == parDict.inHiZ) {
        nucProcessCS->zmax = stod(parValue);
      }
      if (parameter == parDict.inLowM) {
        nucProcessCS->mmin = stod(parValue);
      }
      if (parameter == parDict.inHiM) {
        nucProcessCS->mmax = stod(parValue);
      }
      if (parameter == parDict.inLowY) {
        nucProcessCS->ymin = stod(parValue);
      }
      if (parameter == parDict.inHiY) {
        nucProcessCS->ymax = stod(parValue);
      }
      if (parameter == parDict.inBinsZ) {
        nucProcessCS->nz = stoi(parValue);
      }
      if (parameter == parDict.inBinsM) {
        nucProcessCS->nm = stoi(parValue);
      }
      if (parameter == parDict.inBinsY) {
        nucProcessCS->ny = stoi(parValue);
      }
      if (parameter == parDict.inWSRadius) {
        UpcCrossSection::R = stod(parValue);
      }
      if (parameter == parDict.inWSA) {
        UpcCrossSection::a = stod(parValue);
      }
      if (parameter == parDict.inNucZ) {
        UpcCrossSection::Z = stoi(parValue);
      }
      if (parameter == parDict.inNucA) {
        nucProcessCS->A = stoi(parValue);
      }
      if (parameter == parDict.inFluxPoint) {
        nucProcessCS->isPoint = stoi(parValue);
      }
      if (parameter == parDict.inBreakupMode) {
        nucProcessCS->breakupMode = stoi(parValue);
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
        nucProcessCS->useNonzeroGamPt = stoi(parValue);
      }
      if (parameter == parDict.inPolarized) {
        nucProcessCS->usePolarizedCS = stoi(parValue);
        usePolarizedCS = stoi(parValue);
      }
      if (parameter == parDict.inSeed) {
        seed = stol(parValue);
      }
      if (parameter == parDict.inROOTOut) {
        useROOTOut = stoi(parValue);
      }
      if (parameter == parDict.inHepMCOut) {
        useHepMCOut = stoi(parValue);
      }
      if (parameter == parDict.inDoMCut) {
        nucProcessCS->doMassCut = stoi(parValue);
      }
      if (parameter == parDict.inLowMCut) {
        nucProcessCS->lowMCut = stod(parValue);
      }
      if (parameter == parDict.inHighMCut) {
        nucProcessCS->hiMCut = stod(parValue);
      }
    }
    // do a sanity check for output formats
    if (!useROOTOut && !useHepMCOut) {
      PLOG_FATAL << "Output format not set! Choose ROOT or/and HepMC via flags USE_ROOT_OUTPUT and USE_HEPMC_OUTPUT!";
      std::_Exit(-1);
    }
    fInputs.close();
  } else {
    PLOG_WARNING << "Input file not found! Using default parameters...";
    printParameters();
  }
}

void UpcGenerator::printParameters()
{
  PLOG_INFO << "OMP_NTHREADS " << numThreads;
  PLOG_INFO << "NUCLEUS_Z " << UpcCrossSection::Z;
  PLOG_INFO << "NUCLEUS_A " << nucProcessCS->A;
  PLOG_INFO << "WS_R " << UpcCrossSection::R;
  PLOG_INFO << "WS_A " << UpcCrossSection::a;
  PLOG_INFO << "SQRTS " << UpcCrossSection::sqrts;
  PLOG_INFO << "PROC_ID " << procID;
  PLOG_INFO << "LEP_A " << aLep;
  PLOG_INFO << "NEVENTS " << nEvents;
  PLOG_INFO << "DO_PT_CUT " << doPtCut;
  PLOG_INFO << "PT_MIN " << minPt;
  PLOG_INFO << "DO_ETA_CUT " << doEtaCut;
  PLOG_INFO << "ETA_MIN " << minEta;
  PLOG_INFO << "ETA_MAX " << maxEta;
  PLOG_INFO << "ZMIN " << nucProcessCS->zmin;
  PLOG_INFO << "ZMAX " << nucProcessCS->zmax;
  PLOG_INFO << "MMIN " << nucProcessCS->mmin;
  PLOG_INFO << "MMAX " << nucProcessCS->mmax;
  PLOG_INFO << "YMIN " << nucProcessCS->ymin;
  PLOG_INFO << "YMAX " << nucProcessCS->ymax;
  PLOG_INFO << "BINS_Z " << nucProcessCS->nz;
  PLOG_INFO << "BINS_M " << nucProcessCS->nm;
  PLOG_INFO << "BINS_Y " << nucProcessCS->ny;
  PLOG_INFO << "FLUX_POINT " << nucProcessCS->isPoint;
  PLOG_INFO << "BREAKUP_MODE " << nucProcessCS->breakupMode;
  PLOG_INFO << "NON_ZERO_GAM_PT " << nucProcessCS->useNonzeroGamPt;
  PLOG_INFO << "USE_POLARIZED_CS " << usePolarizedCS;
  PLOG_INFO << "PYTHIA_VERSION " << pythiaVersion;
  PLOG_INFO << "PYTHIA8_FSR " << doFSR;
  PLOG_INFO << "PYTHIA8_DECAYS " << doDecays;
  PLOG_INFO << "SEED " << seed;
  PLOG_INFO << "DO_M_CUT " << nucProcessCS->doMassCut;
  PLOG_INFO << "LOW_M_CUT " << nucProcessCS->lowMCut;
  PLOG_INFO << "HIGH_M_CUT " << nucProcessCS->hiMCut;
}

void UpcGenerator::pairProduction(TLorentzVector& pPair,             // lorentz pair-momentum vector of incoming photons
                                  TVector3& vec,                     // momentum of outgoing particles
                                  double mPart,                      // outgoing particle mass
                                  int partPDG,                       // pdg of final-state particles
                                  bool isCharged,                    // final-state particles are charged or not
                                  vector<TLorentzVector>& particles, // vector for final-state particles
                                  vector<int>& pdgs,
                                  vector<int>& mothers,
                                  vector<int>& statuses)
{
  TLorentzVector tlVec1, tlVec2;
  tlVec1.SetVectM(vec, mPart);
  tlVec2.SetVectM(-vec, mPart);
  particles.emplace_back(tlVec1);
  particles.emplace_back(tlVec2);

  TVector3 boost = pPair.BoostVector();
  particles[0].Boost(boost);
  particles[1].Boost(boost);

  int sign1, sign2;
  if (isCharged) {
    sign1 = gRandom->Uniform(-1, 1) > 0 ? 1 : -1;
    sign2 = -sign1;
  } else {
    sign1 = 1;
    sign2 = 1;
  }

  pdgs.emplace_back(sign1 * partPDG);
  pdgs.emplace_back(sign2 * partPDG);
  mothers.emplace_back(-1);
  mothers.emplace_back(-1);
  statuses.emplace_back(23);
  statuses.emplace_back(23);
}

void UpcGenerator::singleProduction(TLorentzVector& pPair,                  // input lorentz pair-momentum vector
                                    int partPDG,                            // pdg of final-state particle
                                    std::vector<TLorentzVector>& particles, // vector for final-state particles
                                    std::vector<int>& pdgs,
                                    std::vector<int>& mothers,
                                    std::vector<int>& statuses)
{
  particles.emplace_back(pPair);
  // TVector3 boost = pPair.BoostVector();
  // particles[0].Boost(boost);
  pdgs.emplace_back(partPDG);
  mothers.emplace_back(-1);
  statuses.emplace_back(23);
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
    if (ip + 1 > mothers.size())
      mothers.emplace_back(mother);
    particles.emplace_back(tlVector);
  }
#endif
}

// todo: for now valid only for neutral-charged decay products
void UpcGenerator::twoPartDecayUniform(vector<int>& pdgs, // output vectors with pdg codes, status codes, mother IDs and 4-momenta
                                       vector<int>& statuses,
                                       vector<int>& mothers,
                                       vector<TLorentzVector>& particles,
                                       int id,           // index of particle to be decayed
                                       double decayMass, // mass of decay products
                                       int decayProdPDG) // pdg of decay products
{
  const int status = 33; // code from pythia 8: outgoing particles from subsequent subprocesses
  vector<TLorentzVector> decays(2);
  double mDecay = decayMass;
  double ePhot1 = particles[id].Mag() / 2.;
  double pPhot1 = sqrt(ePhot1 * ePhot1 - mDecay * mDecay);
  double phi1 = gRandom->Uniform(0., 2. * M_PI);
  double cost1 = gRandom->Uniform(-1., 1.);
  double theta1 = acos(cost1);
  TVector3 vPhot1;
  vPhot1.SetMagThetaPhi(pPhot1, theta1, phi1);
  decays[0].SetVectM(-vPhot1, mDecay);
  decays[1].SetVectM(vPhot1, mDecay);
  TVector3 boost1 = particles[id].BoostVector();
  TVector3 zAxis1 = particles[id].Vect().Unit();
  decays[0].RotateUz(zAxis1);
  decays[1].RotateUz(zAxis1);
  decays[0].Boost(boost1);
  decays[1].Boost(boost1);
  pdgs.emplace_back(decayProdPDG);
  statuses.emplace_back(status);
  mothers.emplace_back(id);
  particles.emplace_back(decays[0]);
  pdgs.emplace_back(decayProdPDG);
  statuses.emplace_back(status);
  mothers.emplace_back(id);
  particles.emplace_back(decays[1]);
}

bool UpcGenerator::checkKinCuts(std::vector<TLorentzVector>& particles)
{
  bool pass = true;

  for (const auto& tlvec : particles) {
    // check pt cut
    if (doPtCut) {
      double pt = tlvec.Pt();
      if (pt < minPt) {
        pass = false;
        break;
      }
    }
    // check eta cuts
    if (doEtaCut) {
      double eta = tlvec.Eta();
      if (eta < minEta || eta > maxEta) {
        pass = false;
        break;
      }
    }
  }

  return pass;
}

void UpcGenerator::writeEvent(int evt,
                              const vector<int>& pdgs,
                              const vector<int>& statuses,
                              const vector<int>& mothers,
                              const vector<TLorentzVector>& particles)
{
  if (useROOTOut) {
    for (int i = 0; i < particles.size(); i++) {
      particle.eventNumber = evt;
      particle.pdgCode = pdgs[i];
      particle.particleID = i;
      particle.statusID = statuses[i];
      particle.motherID = mothers[i];
      particle.px = particles[i].Px();
      particle.py = particles[i].Py();
      particle.pz = particles[i].Pz();
      particle.e = particles[i].E();
      mOutTree->Fill();
    }
  }

  if (useHepMCOut) {
    writerHepMC->writeEventInfo(evt, particles.size());
    for (int i = 0; i < particles.size(); i++) {
      writerHepMC->writeParticleInfo(i + 1, mothers[i] + 1, pdgs[i],
                                     particles[i].Px(), particles[i].Py(), particles[i].Pz(), particles[i].E(), particles[i].M(),
                                     mothers[i] < 0 ? 1 : 2);
    }
  }
}

void UpcGenerator::generateEvents()
{
  gErrorIgnoreLevel = 6001; // fixme: temporary workaround

  TH2D* hCrossSectionZM = nullptr;
  TH2D* hCrossSectionZMPolS = nullptr;
  TH2D* hCrossSectionZMPolPS = nullptr;

  int nm = nucProcessCS->nm;
  double mmin = nucProcessCS->mmin;
  double mmax = nucProcessCS->mmax;
  double dm = (mmax - mmin) / nm;
  int nz = nucProcessCS->nz;
  double zmin = nucProcessCS->zmin;
  double zmax = nucProcessCS->zmax;
  double dz = (zmax - zmin) / nz;
  int ny = nucProcessCS->ny;
  double ymin = nucProcessCS->ymin;
  double ymax = nucProcessCS->ymax;
  double dy = (ymax - ymin) / ny;

  // mass of a particle in final state
  double mPart = nucProcessCS->elemProcess->mPart;
  int partPDG = nucProcessCS->elemProcess->partPDG;
  bool isCharged = nucProcessCS->elemProcess->isCharged;

  if (usePolarizedCS) {
    hCrossSectionZMPolS = new TH2D("hCrossSectionZMPolS", ";m [gev]; z; cs [nb/gev]",
                                   nz, zmin, zmax,
                                   nm, mmin, mmax);

    nucProcessCS->fillCrossSectionZM(hCrossSectionZMPolS, zmin, zmax, nz, mmin, mmax, nm, 1);

    hCrossSectionZMPolPS = new TH2D("hCrossSectionZMPolPS", ";m [gev]; z; cs [nb/gev]",
                                    nz, zmin, zmax,
                                    nm, mmin, mmax);

    nucProcessCS->fillCrossSectionZM(hCrossSectionZMPolPS, zmin, zmax, nz, mmin, mmax, nm, 2);
  } else {
    hCrossSectionZM = new TH2D("hCrossSectionZM_", ";m [gev]; z; cs [nb/gev]",
                               nz, zmin, zmax,
                               nm, mmin, mmax);

    nucProcessCS->fillCrossSectionZM(hCrossSectionZM, zmin, zmax, nz, mmin, mmax, nm, 0);
    // hCrossSectionZM = ((UpcTwoPhotonDipion*)nucProcessCS->elemProcess)->hCrossSectionZM;
  }

  // calculating nuclear cross section in YM space
  // -----------------------------------------------------------------------
  auto* hNucCSYM = new TH2D("hNucCSYM", "", ny, ymin, ymax, nm, mmin, mmax);
  vector<vector<double>> hPolCSRatio;
  if (usePolarizedCS) {
    hPolCSRatio.resize(ny, vector<double>(nm));
  }
  double totCS;
  nucProcessCS->calcNucCrossSectionYM(hNucCSYM, hPolCSRatio, totCS);

  vector<int> pdgs;
  vector<int> statuses;
  vector<int> mothers;
  vector<TLorentzVector> particles;

  // setting up histograms for sampling
  // -----------------------------------------------------------------------
  auto* hNucCSM = hNucCSYM->ProjectionY();
  vector<TH1D*> hCrossSecsZ(nm, nullptr);
  vector<TH1D*> hCrossSecsZ_S(nm, nullptr);
  vector<TH1D*> hCrossSecsZ_PS(nm, nullptr);
  if (!usePolarizedCS) {
    for (int i = 0; i < nm; i++) {
      hCrossSecsZ[i] = hCrossSectionZM->ProjectionX(Form("hCrossSecsZ_%d", i + 1), i + 1, i + 1);
    }
  } else {
    for (int i = 0; i < nm; i++) {
      hCrossSecsZ_S[i] = hCrossSectionZMPolS->ProjectionX(Form("hCrossSecsZPolS_%d", i + 1), i + 1, i + 1);
      hCrossSecsZ_PS[i] = hCrossSectionZMPolPS->ProjectionX(Form("hCrossSecsZPolPS_%d", i + 1), i + 1, i + 1);
    }
  }

  // hack for ALP: cos(theta) is uniformly distributed
  bool ignoreCSZ = false;
  if (procID == 51) {
    ignoreCSZ = true;
  }

  // generationg events
  // -----------------------------------------------------------------------

  // initialize file output
  if (useROOTOut) {
    mOutFile = new TFile("events.root", "recreate", "", 4 * 100 + 9); // using LZ4 with level 5 compression
    PLOG_WARNING << "Using ROOT tree for output!";
    PLOG_INFO << "Events will be written to events.root";
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
  }

  if (useHepMCOut) {
    writerHepMC = new WriterHepMC("events.hepmc");
  }

  PLOG_INFO << "Generating " << nEvents << " events...";

  long int rejected = 0;
  long int evt = 0;
  while (evt < nEvents) {
    if (debug > 1) {
      PLOG_DEBUG << "Event number: " << evt + 1;
    }
    if (debug <= 1 && ((evt + 1) % 10000 == 0)) {
      PLOG_INFO << "Event number: " << evt + 1;
    }

    // pick pair m and y from nuclear cross section
    double mPair, yPair;
    hNucCSYM->GetRandom2(yPair, mPair);
    int yPairBin = hNucCSYM->GetXaxis()->FindBin(yPair);
    int mPairBin = hNucCSYM->GetYaxis()->FindBin(mPair);

    // pick z = cos(theta) for corresponding m from elem. cross section
    double cost;
    if (!ignoreCSZ) { // ignoring z distribution for narrow resonances
      if (usePolarizedCS) {
        double frac = hPolCSRatio[yPairBin][mPairBin];
        bool pickScalar = gRandom->Uniform(0, 1) < frac;
        if (pickScalar) {
          cost = hCrossSecsZ_S[mPairBin - 1]->GetRandom();
        } else {
          cost = hCrossSecsZ_PS[mPairBin - 1]->GetRandom();
        }
      } else {
        cost = hCrossSecsZ[mPairBin - 1]->GetRandom();
      }
    } else {
      cost = gRandom->Uniform(-1., 1.);
    }

    double theta = acos(cost);
    double phi = gRandom->Uniform(0., 2. * M_PI);

    TLorentzVector pPair;
    nucProcessCS->getPairMomentum(mPair, yPair, pPair);
    double pMag = sqrt(pPair.Mag() * pPair.Mag() / 4 - mPart * mPart);

    TVector3 vec;
    vec.SetMagThetaPhi(pMag, theta, phi);

    if (isPairProduction) {
      pairProduction(pPair, vec, mPart, partPDG, isCharged, particles, pdgs, mothers, statuses);
    }

    if (isSingleProduction) {
      singleProduction(pPair, partPDG, particles, pdgs, mothers, statuses);
    }

    // check kinematic cuts
    if (!checkKinCuts(particles)) {
      rejected++;
      pdgs.clear();
      statuses.clear();
      mothers.clear();
      particles.clear();
      continue;
    }

    // lepton decays for taus
    // todo: at the moment "fsr" and "decays" flags
    //       are only really meaningful for pythia8
    if ((doFSR || doDecays) && isPythiaUsed && (procID >= 11 && procID <= 15)) {
      processInPythia(pdgs, statuses, mothers, particles);
    }

    // uniform pion decays into photon pairs
    if (procID == 111) {
      twoPartDecayUniform(pdgs, statuses, mothers, particles, 0, 0., 22);
      twoPartDecayUniform(pdgs, statuses, mothers, particles, 1, 0., 22);
    }

    // uniform ALP decay into two photons
    if (procID == 51) {
      twoPartDecayUniform(pdgs, statuses, mothers, particles, 0, 0., 22);
    }

    writeEvent(evt, pdgs, statuses, mothers, particles);

    pdgs.clear();
    statuses.clear();
    mothers.clear();
    particles.clear();
    evt++;
  }

  PLOG_INFO << "Event generation is finished!";

  PLOG_INFO << fixed << setprecision(6) << "Total nuclear cross section: " << totCS << " mb";

  bool doAnyCuts = doPtCut || doEtaCut;
  if (doAnyCuts && nEvents > 0) {
    double fidCS = totCS * (double)nEvents / (double)(nEvents + rejected);
    PLOG_INFO << "Kinematic cuts were used";
    PLOG_INFO << "Number of rejected events = " << rejected;
    PLOG_INFO << fixed << setprecision(6) << "Cross section with cuts = " << fidCS << " mb";
  }

  if (useROOTOut) {
    if (debug > 0) {
      hNucCSM->Write();
      hNucCSYM->Write();
    }
    mOutFile->Write();
    mOutFile->Close();
  }

  if (useHepMCOut) {
    delete writerHepMC;
  }

  delete hCrossSectionZM;
  delete hCrossSectionZMPolS;
  delete hCrossSectionZMPolPS;

  delete hNucCSM;
  delete hNucCSYM;

  for (int i = 0; i < nm; i++) {
    delete hCrossSecsZ[i];
    delete hCrossSecsZ_S[i];
    delete hCrossSecsZ_PS[i];
  }
}