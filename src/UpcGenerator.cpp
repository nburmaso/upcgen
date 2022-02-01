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

using namespace std;

// out-of-line initialization for static members
int UpcGenerator::debug = 0;

UpcGenerator::UpcGenerator()
{
  calcMachine = new UpcCalcMachine();
}

void UpcGenerator::init()
{
  // get parameters from file
  initGeneratorFromFile();
  PLOG_WARNING << "Check inputs:";
  printParameters();

  calcMachine->numThreads = numThreads;
  calcMachine->init();
  calcMachine->setElemProcess(procID);

  // if not dummy, set for dilepton photoproduction
  if (aLep > -9999 && (procID >= 10 && procID <= 12)) {
    auto* proc = (UpcTwoPhotonDilep*)calcMachine->elemProcess;
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
        UpcCalcMachine::sqrts = stod(parValue);
        UpcCalcMachine::g1 = UpcCalcMachine::sqrts / (2. * phys_consts::mProt);
        UpcCalcMachine::g2 = UpcCalcMachine::sqrts / (2. * phys_consts::mProt);
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
        calcMachine->zmin = stod(parValue);
      }
      if (parameter == parDict.inHiZ) {
        calcMachine->zmax = stod(parValue);
      }
      if (parameter == parDict.inLowM) {
        calcMachine->mmin = stod(parValue);
      }
      if (parameter == parDict.inHiM) {
        calcMachine->mmax = stod(parValue);
      }
      if (parameter == parDict.inLowY) {
        calcMachine->ymin = stod(parValue);
      }
      if (parameter == parDict.inHiY) {
        calcMachine->ymax = stod(parValue);
      }
      if (parameter == parDict.inBinsZ) {
        calcMachine->nz = stoi(parValue);
      }
      if (parameter == parDict.inBinsM) {
        calcMachine->nm = stoi(parValue);
      }
      if (parameter == parDict.inBinsY) {
        calcMachine->ny = stoi(parValue);
      }
      if (parameter == parDict.inWSRadius) {
        UpcCalcMachine::R = stod(parValue);
      }
      if (parameter == parDict.inWSA) {
        UpcCalcMachine::a = stod(parValue);
      }
      if (parameter == parDict.inNucZ) {
        UpcCalcMachine::Z = stoi(parValue);
      }
      if (parameter == parDict.inNucA) {
        calcMachine->A = stoi(parValue);
      }
      if (parameter == parDict.inFluxPoint) {
        calcMachine->isPoint = stoi(parValue);
      }
      if (parameter == parDict.inBreakupMode) {
        calcMachine->breakupMode = stoi(parValue);
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
        calcMachine->useNonzeroGamPt = stoi(parValue);
      }
      if (parameter == parDict.inPolarized) {
        calcMachine->usePolarizedCS = stoi(parValue);
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
  PLOG_WARNING << "NUCLEUS_Z " << UpcCalcMachine::Z;
  PLOG_WARNING << "NUCLEUS_A " << calcMachine->A;
  PLOG_WARNING << "WS_R " << UpcCalcMachine::R;
  PLOG_WARNING << "WS_A " << UpcCalcMachine::a;
  PLOG_WARNING << "(CALCULATED) WS_RHO0 " << UpcCalcMachine::rho0;
  PLOG_WARNING << "SQRTS " << UpcCalcMachine::sqrts;
  PLOG_WARNING << "PROC_ID " << procID;
  PLOG_WARNING << "LEP_A " << aLep;
  PLOG_WARNING << "NEVENTS " << nEvents;
  PLOG_WARNING << "DO_PT_CUT " << doPtCut;
  PLOG_WARNING << "PT_MIN " << minPt;
  PLOG_WARNING << "ZMIN " << calcMachine->zmin;
  PLOG_WARNING << "ZMAX " << calcMachine->zmax;
  PLOG_WARNING << "MMIN " << calcMachine->mmin;
  PLOG_WARNING << "MMAX " << calcMachine->mmax;
  PLOG_WARNING << "YMIN " << calcMachine->ymin;
  PLOG_WARNING << "YMAX " << calcMachine->ymax;
  PLOG_WARNING << "BINS_Z " << calcMachine->nz;
  PLOG_WARNING << "BINS_M " << calcMachine->nm;
  PLOG_WARNING << "BINS_Y " << calcMachine->ny;
  PLOG_WARNING << "FLUX_POINT " << calcMachine->isPoint;
  PLOG_WARNING << "BREAKUP_MODE " << calcMachine->breakupMode;
  PLOG_WARNING << "NON_ZERO_GAM_PT " << calcMachine->useNonzeroGamPt;
  PLOG_WARNING << "USE_POLARIZED_CS " << usePolarizedCS;
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

void UpcGenerator::generateEvents()
{
  TH2D* hCrossSectionZM;
  TH2D* hCrossSectionZMPolS;
  TH2D* hCrossSectionZMPolPS;

  int nm = calcMachine->nm;
  double mmin = calcMachine->mmin;
  double mmax = calcMachine->mmax;
  double dm = (mmax - mmin) / nm;
  int nz = calcMachine->nz;
  double zmin = calcMachine->zmin;
  double zmax = calcMachine->zmax;
  double dz = (zmax - zmin) / nz;
  int ny = calcMachine->ny;
  double ymin = calcMachine->ymin;
  double ymax = calcMachine->ymax;
  double dy = (ymax - ymin) / ny;

  // mass of a particle in final state
  double mPart = calcMachine->elemProcess->mPart;
  int partPDG = calcMachine->elemProcess->partPDG;

  if (usePolarizedCS) {
    hCrossSectionZMPolS = new TH2D("hCrossSectionZMPolS", ";m [gev]; z; cs [nb/gev]",
                                   nz, zmin, zmax,
                                   nm, mmin, mmax);

    calcMachine->fillCrossSectionZM(hCrossSectionZMPolS, zmin, zmax, nz, mmin, mmax, nm, 1);

    hCrossSectionZMPolPS = new TH2D("hCrossSectionZMPolPS", ";m [gev]; z; cs [nb/gev]",
                                    nz, zmin, zmax,
                                    nm, mmin, mmax);

    calcMachine->fillCrossSectionZM(hCrossSectionZMPolPS, zmin, zmax, nz, mmin, mmax, nm, 2);
  } else {
    hCrossSectionZM = new TH2D("hCrossSectionMZ", ";m [gev]; z; cs [nb/gev]",
                               nz, zmin, zmax,
                               nm, mmin, mmax);

    calcMachine->fillCrossSectionZM(hCrossSectionZM, zmin, zmax, nz, mmin, mmax, nm, 0);
  }

  // calculating nuclear cross section in YM space
  // -----------------------------------------------------------------------
  auto* hNucCSYM = new TH2D("hNucCSYM", "", ny - 1, ymin, ymax, nm - 1, mmin, mmax);
  vector<vector<double>> hPolCSRatio;
  if (usePolarizedCS) {
    hPolCSRatio.resize(ny, vector<double>(nm));
  }
  calcMachine->calcNucCrossSectionYM(hNucCSYM, hPolCSRatio);

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

  // setting up samplers
  // -----------------------------------------------------------------------
  auto* hNucCSM = hNucCSYM->ProjectionY();
  vector<double> mGrid(nm, 0.);
  for (int i = 0; i < nm; i++) {
    mGrid[i] = mmin + i * dm;
  }

  vector<double> yGrid(ny, 0.);
  for (int i = 0; i < ny; i++) {
    yGrid[i] = ymin + i * dy;
  }

  vector<double> zGrid(nz, 0.);
  for (int i = 0; i < nz; i++) {
    zGrid[i] = zmin + i * dz;
  }

  // sampler for nuclear cross section
  auto* samplerNuc = new UpcSampler();
  samplerNuc->setSeed(seed);
  samplerNuc->setDistr1D(hNucCSM, nm);
  samplerNuc->setDistr2D(hNucCSYM, ny, nm);
  samplerNuc->setYGrid(mGrid);
  samplerNuc->setXGrid(yGrid);

  // sampler for unpolarized cs
  UpcSampler* samplerElem;

  // samplers for polarized cs
  UpcSampler* samplerElemS;
  UpcSampler* samplerElemPS;

  if (!usePolarizedCS) {
    samplerElem = new UpcSampler();
    samplerElem->setSeed(seed);
    samplerElem->setDistr2D(hCrossSectionZM, nz, nm);
    samplerElem->setXGrid(zGrid);
    samplerElem->setYGrid(mGrid);
  } else {
    samplerElemS = new UpcSampler();
    samplerElemS->setSeed(seed);
    samplerElemS->setDistr2D(hCrossSectionZMPolS, nz, nm);
    samplerElemS->setXGrid(zGrid);
    samplerElemS->setYGrid(mGrid);

    samplerElemPS = new UpcSampler();
    samplerElemPS->setSeed(seed);
    samplerElemPS->setDistr2D(hCrossSectionZMPolPS, nz, nm);
    samplerElemPS->setXGrid(zGrid);
    samplerElemPS->setYGrid(mGrid);
  }

  // generationg events
  // -----------------------------------------------------------------------

#ifndef USE_HEPMC
  // initialize file output
  PLOG_WARNING << "Using ROOT tree for output!";
  PLOG_INFO << "Events will be written to " << "events.root";
  mOutFile = new TFile("events.root", "recreate", "", 4 * 100 + 5); // using LZ4 with level 5 compression
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
  PLOG_INFO << "Events will be written to " << "events.hepmc";
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
    int yPairBin;
    int mPairBin;
    samplerNuc->sample2D(mPairBin, yPairBin);
    double mPair = mGrid[mPairBin];
    double yPair = yGrid[yPairBin];

    // pick z = cos(theta) for corresponding m from elem. cross section
    double cost;
    if (usePolarizedCS) {
      double frac = hPolCSRatio[yPairBin][mPairBin];
      bool pickScalar = gRandom->Uniform(0, 1) < frac;
      if (pickScalar) {
        samplerElemS->sampleXAtY(mPairBin, cost);
      } else {
        samplerElemPS->sampleXAtY(mPairBin, cost);
      }
    } else {
      samplerElem->sampleXAtY(mPairBin, cost);
    }

    double theta = TMath::ACos(cost);
    double phi = gRandom->Uniform(0., 2. * M_PI);

    TLorentzVector pPair;
    calcMachine->getPairMomentum(mPair, yPair, pPair);
    double pMag = TMath::Sqrt(pPair.Mag() * pPair.Mag() / 4 - mPart * mPart);
    TVector3 vLep;
    vLep.SetMagThetaPhi(pMag, theta, phi);
    TLorentzVector tlLep1, tlLep2;
    tlLep1.SetVectM(vLep, mPart);
    tlLep2.SetVectM(-vLep, mPart);
    particles.emplace_back(tlLep1);
    particles.emplace_back(tlLep2);

    TVector3 boost = pPair.BoostVector();
    particles[0].Boost(boost);
    particles[1].Boost(boost);

    int sign = gRandom->Uniform(-1, 1) > 0 ? 1 : -1;
    pdgs.emplace_back(sign * partPDG);
    pdgs.emplace_back(-sign * partPDG);
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
  }

#ifndef USE_HEPMC
  mOutFile->Write();
  mOutFile->Close();
#endif

#ifdef USE_HEPMC
  writerHepMC->close();
#endif
}