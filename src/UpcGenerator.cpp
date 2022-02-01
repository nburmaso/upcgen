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
#include "UpcSampler.h"

using namespace std;

// out-of-line initialization for static members
int UpcGenerator::debug = 0;
std::map<int, double> UpcGenerator::lepMassMap;

UpcGenerator::UpcGenerator()
{
  // populating lepton mass map
  // data from PDG
  lepMassMap[11] = 0.000510998946;
  lepMassMap[13] = 0.1056583745;
  lepMassMap[15] = 1.77686;
  calcMachine = new UpcCalcMachine();
}

void UpcGenerator::init()
{
  // get parameters from file
  initGeneratorFromFile();

  calcMachine->mLep = lepMassMap[lepPDG];
  calcMachine->numThreads = numThreads;
  calcMachine->init();

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

  PLOG_WARNING << "Check inputs:";
  printParameters();
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
        UpcCalcMachine::g1 = UpcCalcMachine::sqrts / (2. * UpcCalcMachine::mProt);
        UpcCalcMachine::g2 = UpcCalcMachine::sqrts / (2. * UpcCalcMachine::mProt);
      }
      if (parameter == parDict.inLepPDG) {
        lepPDG = stoi(parValue);
      }
      if (parameter == parDict.inLepA) {
        calcMachine->aLep = stod(parValue);
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
  PLOG_WARNING << "LEP_PDG " << lepPDG;
  PLOG_WARNING << "LEP_A " << calcMachine->aLep;
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
  TH2D* hCrossSectionMZ;
  TH2D* hCrossSectionMZPolS;
  TH2D* hCrossSectionMZPolPS;

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

  double aLep = calcMachine->aLep;
  double mLep = calcMachine->mLep;

  if (usePolarizedCS) {
    hCrossSectionMZPolS = new TH2D("hCrossSectionMZPolS", ";m [gev]; z; cs [nb/gev]",
                                   nm, mmin, mmax,
                                   nz, zmin, zmax);

    calcMachine->fillCrossSectionMZ(hCrossSectionMZPolS, mmin, mmax, nm, zmin, zmax, nz, 1);

    hCrossSectionMZPolPS = new TH2D("hCrossSectionMZPolPs", ";m [gev]; z; cs [nb/gev]",
                                    nm, mmin, mmax,
                                    nz, zmin, zmax);

    calcMachine->fillCrossSectionMZ(hCrossSectionMZPolPS, mmin, mmax, nm, zmin, zmax, nz, 2);
  } else {
    hCrossSectionMZ = new TH2D("hCrossSectionMZ", ";m [gev]; z; cs [nb/gev]",
                               nm, mmin, mmax,
                               nz, zmin, zmax);

    calcMachine->fillCrossSectionMZ(hCrossSectionMZ, mmin, mmax, nm, zmin, zmax, nz, 0);
  }

  // calculating nuclear cross section in YM space
  // -----------------------------------------------------------------------
  auto* hNucCSYM = new TH2D("hNucCSYM", "", ny - 1, ymin, ymax, nm - 1, mmin, mmax);
  vector<vector<double>> hPolCSRatio;
  if (usePolarizedCS) {
    hPolCSRatio.resize(ny, vector<double>(nm));
  }
  calcMachine->calcNucCrossSectionYM(hNucCSYM, hPolCSRatio);

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

  if (doPtCut) {
    for (int mBin = 1; mBin <= nm; mBin++) {
      double zCut = cutsZ[mBin - 1];
      int zCutBinUp = hCrossSectionMZ->GetXaxis()->FindBin(zCut);
      int zCutBinLow = hCrossSectionMZ->GetXaxis()->FindBin(-zCut);
      for (int zBin = zCutBinUp + 1; zBin <= nz; zBin++) {
        hCrossSectionMZ->SetBinContent(mBin, zBin, 0);
      }
      for (int zBin = 1; zBin < zCutBinLow; zBin++) {
        hCrossSectionMZ->SetBinContent(mBin, zBin, 0);
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
  TH2D* hCrossSectionZM;
  UpcSampler* samplerElem;

  // samplers for polarized cs
  TH2D* hCrossSectionZMPolS;
  TH2D* hCrossSectionZMPolPS;
  UpcSampler* samplerElemS;
  UpcSampler* samplerElemPS;

  if (!usePolarizedCS) {
    hCrossSectionZM = new TH2D("hCrossSectionZM", "",
                               nz, zmin, zmax,
                               nm, mmin, mmax);
    for (int i = 1; i <= nz; i++) {
      for (int j = 1; j <= nm; j++) {
        hCrossSectionZM->SetBinContent(i, j, hCrossSectionMZ->GetBinContent(j, i));
      }
    }
    samplerElem = new UpcSampler();
    samplerElem->setSeed(seed);
    samplerElem->setDistr2D(hCrossSectionZM, nz, nm);
    samplerElem->setXGrid(zGrid);
    samplerElem->setYGrid(mGrid);
  } else {
    hCrossSectionZMPolS = new TH2D("hCrossSectionZMPolS", "",
                               nz, zmin, zmax,
                               nm, mmin, mmax);
    hCrossSectionZMPolPS = new TH2D("hCrossSectionZMPolPS", "",
                                   nz, zmin, zmax,
                                   nm, mmin, mmax);
    for (int i = 1; i <= nz; i++) {
      for (int j = 1; j <= nm; j++) {
        hCrossSectionZMPolS->SetBinContent(i, j, hCrossSectionMZPolS->GetBinContent(j, i));
        hCrossSectionZMPolPS->SetBinContent(i, j, hCrossSectionMZPolPS->GetBinContent(j, i));
      }
    }

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
  }

#ifndef USE_HEPMC
  mOutFile->Write();
  mOutFile->Close();
#endif

#ifdef USE_HEPMC
  writerHepMC->close();
#endif
}