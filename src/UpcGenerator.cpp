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
  // get parameters from file
  initGeneratorFromFile();

  // sanity checks: lbyl
  if (procID == 1) {
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
  if (procID == 20) {
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

  nucProcessCS->numThreads = numThreads;
  nucProcessCS->init();
  nucProcessCS->setElemProcess(procID);

  PLOG_WARNING << "Check inputs:";
  printParameters();

  // if not dummy, set for dilepton photoproduction
  if (aLep > -9999 && (procID >= 10 && procID <= 12)) {
    auto* proc = (UpcTwoPhotonDilep*)nucProcessCS->elemProcess;
    proc->aLep = aLep;
  }

  // initialize the MT64 random number generator
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
    decayer->setSeed(seed);
    decayer = new UpcPythia6Helper();
  }
#endif

  if (!isPythiaUsed) {
    PLOG_WARNING << "Decays with Pythia are not used!";
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
      if (parameter == parDict.inDoPtCut) {
        doPtCut = stoi(parValue);
      }
      if (parameter == parDict.inLowPt) {
        minPt = stod(parValue);
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
  PLOG_INFO << "(CALCULATED) WS_RHO0 " << UpcCrossSection::rho0;
  PLOG_INFO << "SQRTS " << UpcCrossSection::sqrts;
  PLOG_INFO << "PROC_ID " << procID;
  PLOG_INFO << "LEP_A " << aLep;
  PLOG_INFO << "NEVENTS " << nEvents;
  PLOG_INFO << "DO_PT_CUT " << doPtCut;
  PLOG_INFO << "PT_MIN " << minPt;
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

void UpcGenerator::processPions(vector<int>& pdgs,
                                vector<int>& statuses,
                                vector<int>& mothers,
                                vector<TLorentzVector>& particles)
{
  const int status = 33; // code from pythia 8: outgoing particles from subsequent subprocesses
  vector<TLorentzVector> photons(4);

  // first pion
  double mPhot1 = 0.;
  double ePhot1 = particles[0].Mag() / 2.;
  double pPhot1 = sqrt(ePhot1 * ePhot1 - mPhot1 * mPhot1);
  double phi1 = gRandom->Uniform(0., 2. * M_PI);
  double cost1 = gRandom->Uniform(-1., 1.);
  double theta1 = acos(cost1);
  TVector3 vPhot1;
  vPhot1.SetMagThetaPhi(pPhot1, theta1, phi1);
  photons[0].SetVectM(-vPhot1, mPhot1);
  photons[1].SetVectM(vPhot1, mPhot1);
  TVector3 boost1 = particles[0].BoostVector();
  TVector3 zAxis1 = particles[0].Vect().Unit();
  photons[0].RotateUz(zAxis1);
  photons[1].RotateUz(zAxis1);
  photons[0].Boost(boost1);
  photons[1].Boost(boost1);
  pdgs.emplace_back(22);
  statuses.emplace_back(status);
  mothers.emplace_back(0);
  particles.emplace_back(photons[0]);
  pdgs.emplace_back(22);
  statuses.emplace_back(status);
  mothers.emplace_back(0);
  particles.emplace_back(photons[1]);

  // second pion
  double mPhot2 = 0.;
  double ePhot2 = particles[1].Mag() / 2.;
  double pPhot2 = sqrt(ePhot2 * ePhot2 - mPhot2 * mPhot2);
  double phi2 = gRandom->Uniform(0., 2. * M_PI);
  double cost2 = gRandom->Uniform(-1., 1.);
  double theta2 = acos(cost2);
  TVector3 vPhot2;
  vPhot2.SetMagThetaPhi(pPhot2, theta2, phi2);
  photons[2].SetVectM(-vPhot2, mPhot2);
  photons[3].SetVectM(vPhot2, mPhot2);
  TVector3 boost2 = particles[1].BoostVector();
  TVector3 zAxis2 = particles[1].Vect().Unit();
  photons[2].RotateUz(zAxis2);
  photons[3].RotateUz(zAxis2);
  photons[2].Boost(boost2);
  photons[3].Boost(boost2);
  pdgs.emplace_back(22);
  statuses.emplace_back(status);
  mothers.emplace_back(1);
  particles.emplace_back(photons[2]);
  pdgs.emplace_back(22);
  statuses.emplace_back(status);
  mothers.emplace_back(1);
  particles.emplace_back(photons[3]);
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
    particle.motherID = i <= 1 ? -1 : mothers[i]; // first two particles are primary
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

void UpcGenerator::generateEvents()
{
  gErrorIgnoreLevel = 6001; // fixme: temporary workaround

  TH2D* hCrossSectionZM;
  TH2D* hCrossSectionZMPolS;
  TH2D* hCrossSectionZMPolPS;

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
  auto* hNucCSYM = new TH2D("hNucCSYM", "", ny - 1, ymin, ymax, nm - 1, mmin, mmax);
  vector<vector<double>> hPolCSRatio;
  if (usePolarizedCS) {
    hPolCSRatio.resize(ny, vector<double>(nm));
  }
  nucProcessCS->calcNucCrossSectionYM(hNucCSYM, hPolCSRatio);

  vector<double> cutsZ(hNucCSYM->GetNbinsY());
  if (doPtCut) {
    for (int mBin = 1; mBin <= nm; mBin++) {
      double mass = hNucCSYM->GetYaxis()->GetBinLowEdge(mBin);
      double sqt = TMath::Sqrt(mass * mass * 0.25 - mPart * mPart);
      if (sqt <= minPt) {
        cutsZ[mBin - 1] = 0;
        for (int yBin = 1; yBin <= ny; yBin++) {
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

  if (doPtCut) {
    for (int mBin = 1; mBin <= nm; mBin++) {
      double zCut = cutsZ[mBin - 1];
      int zCutBinUp = hCrossSectionZM->GetXaxis()->FindBin(zCut);
      int zCutBinLow = hCrossSectionZM->GetXaxis()->FindBin(-zCut);
      for (int zBin = zCutBinUp + 1; zBin <= nz; zBin++) {
        hCrossSectionZM->SetBinContent(zBin, mBin, 0);
      }
      for (int zBin = 1; zBin < zCutBinLow; zBin++) {
        hCrossSectionZM->SetBinContent(zBin, mBin, 0);
      }
      if (usePolarizedCS) {
        for (int zBin = zCutBinUp + 1; zBin <= nz; zBin++) {
          hCrossSectionZMPolS->SetBinContent(zBin, mBin, 0);
          hCrossSectionZMPolPS->SetBinContent(zBin, mBin, 0);
        }
        for (int zBin = 1; zBin < zCutBinLow; zBin++) {
          hCrossSectionZMPolS->SetBinContent(zBin, mBin, 0);
          hCrossSectionZMPolPS->SetBinContent(zBin, mBin, 0);
        }
      }
    }
  }

  // setting up histograms for sampling
  // -----------------------------------------------------------------------
  auto* hNucCSM = hNucCSYM->ProjectionY();
  vector<TH1D*> hCrossSecsZ(nm);
  vector<TH1D*> hCrossSecsZ_S(nm);
  vector<TH1D*> hCrossSecsZ_PS(nm);
  if (!usePolarizedCS) {
    for (int i = 0; i < nm; i++) {
      hCrossSecsZ[i] = hCrossSectionZM->ProjectionX(Form("hCrossSecsZ_%d", i), i, i);
    }
  } else {
    for (int i = 0; i < nm; i++) {
      hCrossSecsZ_S[i] = hCrossSectionZMPolS->ProjectionX(Form("hCrossSecsZPolS_%d", i), i, i);
      hCrossSecsZ_PS[i] = hCrossSectionZMPolPS->ProjectionX(Form("hCrossSecsZPolPS_%d", i), i, i);
    }
  }

  // generationg events
  // -----------------------------------------------------------------------

#ifndef USE_HEPMC
  // initialize file output
  PLOG_WARNING << "Using ROOT tree for output!";
  PLOG_INFO << "Events will be written to "
            << "events.root";
  mOutTree = new TTree("particles", "Generated particles");
  mOutTree->SetDirectory(nullptr);
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
  PLOG_INFO << "Events will be written to "
            << "events.hepmc";
  writerHepMC = new HepMC3::WriterAscii("events.hepmc");
#endif

  PLOG_INFO << "Generating " << nEvents << " events...";

  for (int evt = 0; evt < nEvents; evt++) {
    if (debug > 1) {
      PLOG_DEBUG << "Event number: " << evt + 1;
    } else if ((evt + 1) % 10000 == 0) {
      PLOG_INFO << "Event number: " << evt + 1;
    }

    // pick pair m and y from nuclear cross section
    double mPair, yPair;
    hNucCSYM->GetRandom2(yPair, mPair);
    int yPairBin = hNucCSYM->GetXaxis()->FindBin(yPair);
    int mPairBin = hNucCSYM->GetYaxis()->FindBin(mPair);

    // pick z = cos(theta) for corresponding m from elem. cross section
    double cost;
    if (usePolarizedCS) {
      double frac = hPolCSRatio[yPairBin][mPairBin];
      bool pickScalar = gRandom->Uniform(0, 1) < frac;
      if (pickScalar) {
        cost = hCrossSecsZ_S[mPairBin]->GetRandom();
      } else {
        cost = hCrossSecsZ_PS[mPairBin]->GetRandom();
      }
    } else {
      cost = hCrossSecsZ[mPairBin]->GetRandom();
    }

    double theta1 = acos(cost);
    double theta2 = acos(-cost);
    double phi1 = gRandom->Uniform(0., 2. * M_PI);
    double phi2 = phi1 > M_PI ? phi1 - M_PI : phi1 + M_PI;

    TLorentzVector pPair;
    nucProcessCS->getPairMomentum(mPair, yPair, pPair);
    double pMag = sqrt(pPair.Mag() * pPair.Mag() / 4 - mPart * mPart);

    TVector3 vec1, vec2;
    vec1.SetMagThetaPhi(pMag, theta1, phi1);
    vec2.SetMagThetaPhi(pMag, theta2, phi2);

    TLorentzVector tlVec1, tlVec2;
    tlVec1.SetVectM(vec1, mPart);
    tlVec2.SetVectM(vec2, mPart);
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

    // lepton decays for taus
    // todo: at the moment "fsr" and "decays" flags
    //       are only really meaningful for pythia8
    if ((doFSR || doDecays) && isPythiaUsed && (procID >= 10 && procID <= 12)) {
      processInPythia(pdgs, statuses, mothers, particles);
    }

    // uniform pion decays into photon pairs
    if (procID == 20) {
      processPions(pdgs, statuses, mothers, particles);
    }

    writeEvent(evt, pdgs, statuses, mothers, particles);

    pdgs.clear();
    statuses.clear();
    mothers.clear();
    particles.clear();
  }

#ifndef USE_HEPMC
  mOutFile = new TFile("events.root", "recreate", "", 4 * 100 + 9); // using LZ4 with level 5 compression
  mOutTree->Write();
  if (debug > 0) {
    hNucCSM->Write();
    hNucCSYM->Write();
  }
  mOutFile->Write();
  mOutFile->Close();
#endif

#ifdef USE_HEPMC
  writerHepMC->close();
#endif
}