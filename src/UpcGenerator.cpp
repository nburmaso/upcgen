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

#include "UpcGenerator.h"
#include "UpcPhysConstants.h"
#include "UpcSampler.h"
#include "UpcTwoPhotonDipion.h"

// out-of-line initialization for static members
int UpcGenerator::debug = 0;

UpcGenerator::UpcGenerator()
{
  nucProcessCS = new UpcCrossSection();

  // reset nuclear cross section
  totCS = 0.0;
  fidCS = 0.0;

  // set defaults for the collisions system
  setCollisionSystem(5020., 82, 208);
}

UpcGenerator::~UpcGenerator()
{
#if defined(USE_PYTHIA6) || defined(USE_PYTHIA8)
  delete decayer;
#endif
  delete gRandom;
  delete samplerCsYM;
  for (auto* s : samplersCsZ)
    delete s;
  for (auto* s : samplersCsSZ)
    delete s;
  for (auto* s : samplersCsPsZ)
    delete s;
}

void UpcGenerator::init()
{
  TH1::AddDirectory(false);

  if (nucProcessCS == nullptr) {
    PLOG_FATAL << "UpcCrossSection was not initialized! Exiting...";
    std::_Exit(-1);
  }

  ignoreCSZ = false;

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
    ignoreCSZ = true; // hack for ALP: cos(theta) is uniformly distributed
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

  if (procID == 443 || procID == 100443 || procID == 553) {
    isSingleProduction = true; // one particle
    isPairProductionVM = true; // two-part decay
    ignoreCSZ = true;
    auto* proc = (UpcPhotoNuclearVM*)nucProcessCS->elemProcess;
    double mPart = proc->mPart;
    nucProcessCS->mmin = mPart - 1e-6;
    nucProcessCS->mmax = mPart + 1e-6;
    nucProcessCS->nm = 1;
  }

  nucProcessCS->numThreads = numThreads;

  PLOG_WARNING << "Check inputs:";
  printParameters();

  // calculate two-photon luminosity
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

  // compute the nuclear cross section
  computeNuclXsection();
}

void UpcGenerator::setParameterValue(const std::string& parameter, const std::string& parValue)
{
  if (parameter == "NEVENTS") {
    nEvents = stol(parValue);
  }
  if (parameter == "SQRTS") {
    UpcCrossSection::sqrts = stod(parValue);
    UpcCrossSection::g1 = UpcCrossSection::sqrts / (2. * phys_consts::mProt);
    UpcCrossSection::g2 = UpcCrossSection::sqrts / (2. * phys_consts::mProt);
  }
  if (parameter == "PROC_ID") {
    procID = stoi(parValue);
  }
  if (parameter == "LEP_A") {
    aLep = stod(parValue);
  }
  if (parameter == "ALP_MASS") {
    nucProcessCS->alpMass = stod(parValue);
  }
  if (parameter == "ALP_WIDTH") {
    nucProcessCS->alpWidth = stod(parValue);
  }
  if (parameter == "DO_PT_CUT") {
    doPtCut = stoi(parValue);
  }
  if (parameter == "PT_MIN") {
    minPt = stod(parValue);
  }
  if (parameter == "DO_ETA_CUT") {
    doEtaCut = stoi(parValue);
  }
  if (parameter == "ETA_MIN") {
    minEta = stod(parValue);
  }
  if (parameter == "ETA_MAX") {
    maxEta = stod(parValue);
  }
  if (parameter == "ZMIN") {
    nucProcessCS->zmin = stod(parValue);
  }
  if (parameter == "ZMAX") {
    nucProcessCS->zmax = stod(parValue);
  }
  if (parameter == "MMIN") {
    nucProcessCS->mmin = stod(parValue);
  }
  if (parameter == "MMAX") {
    nucProcessCS->mmax = stod(parValue);
  }
  if (parameter == "YMIN") {
    nucProcessCS->ymin = stod(parValue);
  }
  if (parameter == "YMAX") {
    nucProcessCS->ymax = stod(parValue);
  }
  if (parameter == "BINS_Z") {
    nucProcessCS->nz = stoi(parValue);
  }
  if (parameter == "BINS_M") {
    nucProcessCS->nm = stoi(parValue);
  }
  if (parameter == "BINS_Y") {
    nucProcessCS->ny = stoi(parValue);
  }
  if (parameter == "WS_R") {
    UpcCrossSection::R = stod(parValue);
  }
  if (parameter == "WS_A") {
    UpcCrossSection::a = stod(parValue);
  }
  if (parameter == "NUCLEUS_Z") {
    UpcCrossSection::Z = stoi(parValue);
  }
  if (parameter == "NUCLEUS_A") {
    UpcCrossSection::A = stoi(parValue);
  }
  if (parameter == "FLUX_POINT") {
    nucProcessCS->isPoint = stoi(parValue);
  }
  if (parameter == "BREAKUP_MODE") {
    nucProcessCS->breakupMode = stoi(parValue);
  }
  if (parameter == "PYTHIA_VERSION") {
    pythiaVersion = stoi(parValue);
  }
  if (parameter == "PYTHIA8_FSR") {
    doFSR = stoi(parValue);
  }
  if (parameter == "PYTHIA8_DECAYS") {
    doDecays = stoi(parValue);
  }
  if (parameter == "NON_ZERO_GAM_PT") {
    nucProcessCS->useNonzeroGamPt = stoi(parValue);
  }
  if (parameter == "USE_POLARIZED_CS") {
    nucProcessCS->usePolarizedCS = stoi(parValue);
    usePolarizedCS = stoi(parValue);
  }
  if (parameter == "SEED") {
    seed = stol(parValue);
  }
  if (parameter == "USE_ROOT_OUTPUT") {
    useROOTOut = stoi(parValue);
  }
  if (parameter == "USE_HEPMC_OUTPUT") {
    useHepMCOut = stoi(parValue);
  }
  if (parameter == "DO_M_CUT") {
    nucProcessCS->doMassCut = stoi(parValue);
  }
  if (parameter == "LOW_M_CUT") {
    nucProcessCS->lowMCut = stod(parValue);
  }
  if (parameter == "HIGH_M_CUT") {
    nucProcessCS->hiMCut = stod(parValue);
  }
  if (parameter == "SHADOWING") {
    nucProcessCS->shadowingOption = stoi(parValue);
  }
  if (parameter == "DECAY_PDG") {
    nucProcessCS->dghtPDG = stoi(parValue);
  }
}

void UpcGenerator::setCollisionSystem(float sqrts, int nucl_z, int nucl_a)
{
  UpcCrossSection::sqrts = sqrts;
  UpcCrossSection::g1 = UpcCrossSection::sqrts / (2. * phys_consts::mProt);
  UpcCrossSection::g2 = UpcCrossSection::sqrts / (2. * phys_consts::mProt);

  UpcCrossSection::Z = nucl_z;
  UpcCrossSection::A = nucl_a;
  UpcCrossSection::mNucl = (nucl_z * phys_consts::mProt + (nucl_a - nucl_z) * phys_consts::mNeut) / nucl_a;
}

void UpcGenerator::configGeneratorFromFile()
{
  // todo: use <any> from c++17 for a neat parsing?
  if (!gSystem->AccessPathName(parFileName.c_str())) {
    PLOG_INFO << "Reading parameters from " << parFileName << " ...";
    std::ifstream fInputs(parFileName);
    std::string line;
    std::string parameter;
    std::string parValue;
    while (getline(fInputs, line)) {
      std::istringstream iss(line);
      // skip comment lines
      if (line[0] == '#') {
        continue;
      }
      line.erase(std::find(line.begin(), line.end(), '#'), line.end()); // trim comments
      iss >> parameter >> parValue;

      setParameterValue(parameter, parValue);
    }
    // do a sanity check for output formats
    if (!useROOTOut && !useHepMCOut) {
      PLOG_WARNING << "Output format not set! Choose ROOT or/and HepMC via flags USE_ROOT_OUTPUT and USE_HEPMC_OUTPUT!";
    }
    fInputs.close();
  } else {
    PLOG_WARNING << "Input file not found! Using default parameters...";
  }
}

void UpcGenerator::printParameters()
{
  PLOG_INFO << "OMP_NTHREADS " << numThreads;
  PLOG_INFO << "NUCLEUS_Z " << UpcCrossSection::Z;
  PLOG_INFO << "NUCLEUS_A " << UpcCrossSection::A;
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
  PLOG_INFO << "SHADOWING " << nucProcessCS->shadowingOption;
  PLOG_INFO << "DECAY_PDG " << nucProcessCS->dghtPDG;
}

void UpcGenerator::pairProduction(TLorentzVector& pPair,                  // lorentz pair-momentum vector of incoming photons
                                  TVector3& vec,                          // momentum of outgoing particles
                                  double mPart,                           // outgoing particle mass
                                  int partPDG,                            // pdg of final-state particles
                                  bool isCharged,                         // final-state particles are charged or not
                                  std::vector<TLorentzVector>& particles, // vector for final-state particles
                                  std::vector<int>& pdgs,
                                  std::vector<int>& mothers,
                                  std::vector<int>& statuses)
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
  mothers.emplace_back(0);
  mothers.emplace_back(0);
  statuses.emplace_back(23);
  statuses.emplace_back(23);
}

void UpcGenerator::twoPartDecayVM(std::vector<int>& pdgs, // output vectors with pdg codes, status codes, mother IDs and 4-momenta
                                  std::vector<int>& statuses,
                                  std::vector<int>& mothers,
                                  std::vector<TLorentzVector>& particles,
                                  int id) // index of particle to be decayed
{
  int decayProdPDG = nucProcessCS->elemProcess->dghtPDG;
  // J. Breitweg et al., Eur. Phys. J. C2, 247 (1998)
  double theta;
  double dndtheta;
  while (true) {
    theta = M_PI * gRandom->Uniform(0., 1.);
    double test = gRandom->Uniform(0., 1.);
    if (decayProdPDG == 11 || decayProdPDG == 13) {
      dndtheta = std::sin(theta) * (1. + (std::cos(theta) * std::cos(theta)));
    } else if (decayProdPDG == 2212) {
      dndtheta = std::sin(theta) * (1. + (0.605 * std::cos(theta) * std::cos(theta)));
    }
    if (test < dndtheta)
      break;
  }
  int sign1 = gRandom->Uniform(-1, 1) > 0 ? 1 : -1;
  int sign2 = -sign1;
  const int status = 33; // code from pythia 8: outgoing particles from subsequent subprocesses
  std::vector<TLorentzVector> decays(2);
  const auto& particle = particles[id-1];
  double mDecay = nucProcessCS->elemProcess->mDght;
  double pMag = std::sqrt(particle.Mag2() / 4. - mDecay * mDecay);
  double phi = gRandom->Uniform(0., 2. * M_PI);
  TVector3 vec;
  vec.SetMagThetaPhi(pMag, theta, phi);
  decays[0].SetVectM(-vec, mDecay);
  decays[1].SetVectM(vec, mDecay);
  TVector3 boost = particle.BoostVector();
  TVector3 zAxis1 = particle.Vect().Unit();
  decays[0].RotateUz(zAxis1);
  decays[1].RotateUz(zAxis1);
  decays[0].Boost(boost);
  decays[1].Boost(boost);
  pdgs.emplace_back(sign1 * decayProdPDG);
  statuses.emplace_back(status);
  mothers.emplace_back(id);
  particles.emplace_back(decays[0]);
  pdgs.emplace_back(sign2 * decayProdPDG);
  statuses.emplace_back(status);
  mothers.emplace_back(id);
  particles.emplace_back(decays[1]);
}

void UpcGenerator::singleProduction(TLorentzVector& pPair,                  // input lorentz pair-momentum vector
                                    int partPDG,                            // pdg of final-state particle
                                    std::vector<TLorentzVector>& particles, // vector for final-state particles
                                    std::vector<int>& pdgs,
                                    std::vector<int>& mothers,
                                    std::vector<int>& statuses)
{
  particles.emplace_back(pPair);
  pdgs.emplace_back(partPDG);
  mothers.emplace_back(0);
  statuses.emplace_back(23);
}

void UpcGenerator::processInPythia(std::vector<int>& pdgs,
                                   std::vector<int>& statuses,
                                   std::vector<int>& mothers,
                                   std::vector<std::vector<int>>& daughters,
                                   std::vector<TLorentzVector>& particles)
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
  daughters.clear();
  particles.clear();
  for (int ip = 0; ip < processedParts.GetEntriesFast(); ip++) {
    auto* part = (TParticle*)processedParts.At(ip);
    if (debug > 1) {
      PLOG_DEBUG << "Particle info:";
      part->Print();
    }
    int pdg = part->GetPdgCode();
    int status = part->GetStatusCode();
    int mother = abs(status) == 23 ? 0 : part->GetFirstMother() + 1;
    int fDaughter = part->GetFirstDaughter();
    int lDaughter = part->GetLastDaughter();
    part->Momentum(tlVector);
    pdgs.emplace_back(pdg);
    statuses.emplace_back(status);
    mothers.emplace_back(mother);
    daughters.emplace_back(std::vector<int>{fDaughter, lDaughter});
    particles.emplace_back(tlVector);
  }
#endif
}

// todo: for now valid only for neutral-charged decay products
void UpcGenerator::twoPartDecayUniform(std::vector<int>& pdgs, // output vectors with pdg codes, status codes, mother IDs and 4-momenta
                                       std::vector<int>& statuses,
                                       std::vector<int>& mothers,
                                       std::vector<TLorentzVector>& particles,
                                       int id,           // index of particle to be decayed
                                       double decayMass, // mass of decay products
                                       int decayProdPDG) // pdg of decay products
{
  const int status = 33; // code from pythia 8: outgoing particles from subsequent subprocesses
  std::vector<TLorentzVector> decays(2);
  const auto& particle = particles[id-1];
  double mDecay = decayMass;
  double ePhot1 = particle.Mag() / 2.;
  double pPhot1 = sqrt(ePhot1 * ePhot1 - mDecay * mDecay);
  double phi1 = gRandom->Uniform(0., 2. * M_PI);
  double cost1 = gRandom->Uniform(-1., 1.);
  double theta1 = acos(cost1);
  TVector3 vPhot1;
  vPhot1.SetMagThetaPhi(pPhot1, theta1, phi1);
  decays[0].SetVectM(-vPhot1, mDecay);
  decays[1].SetVectM(vPhot1, mDecay);
  TVector3 boost1 = particle.BoostVector();
  TVector3 zAxis1 = particle.Vect().Unit();
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

void UpcGenerator::writeEvent(long int evt,
                              const std::vector<int>& pdgs,
                              const std::vector<int>& statuses,
                              const std::vector<int>& mothers,
                              const std::vector<TLorentzVector>& particles)
{
  if (useROOTOut) {
    for (int i = 0; i < particles.size(); i++) {
      particle.eventNumber = evt;
      particle.pdgCode = pdgs[i];
      particle.particleID = i + 1;
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
    // compute number of vertices
    int nVertices = -1;
    int lastMotherId = -1;
    for (auto mother : mothers) {
      if (mother != lastMotherId) {
        nVertices++;
        lastMotherId = mother;
      }
    }
    writerHepMC->writeEventInfo(evt, static_cast<int>(particles.size()), nVertices);

    for (int i = 0; i < particles.size(); ++i) {
      writerHepMC->writeParticleInfo(i + 1, mothers[i], pdgs[i],
                                     particles[i].Px(), particles[i].Py(), particles[i].Pz(), particles[i].E(), particles[i].M(),
                                     statuses[i]);
    }
  }
}

void UpcGenerator::computeNuclXsection()
{
  bool isVM = false;
  if (procID == 443 || procID == 100443 || procID == 553)
    isVM = true;

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

  // setup for two-photon processes
  if (!isVM) {
    std::vector<std::vector<double>> crossSectionZM;
    std::vector<std::vector<double>> crossSectionZMPolS;
    std::vector<std::vector<double>> crossSectionZMPolPS;

    if (usePolarizedCS) {
      crossSectionZMPolS.resize(nm, std::vector<double>(nz, 0.));
      crossSectionZMPolPS.resize(nm, std::vector<double>(nz, 0.));
      nucProcessCS->fillCrossSectionZM(crossSectionZMPolS, zmin, zmax, nz, mmin, mmax, nm, 1);
      nucProcessCS->fillCrossSectionZM(crossSectionZMPolPS, zmin, zmax, nz, mmin, mmax, nm, 2);
    } else {
      crossSectionZM.resize(nm, std::vector<double>(nz, 0.));
      nucProcessCS->fillCrossSectionZM(crossSectionZM, zmin, zmax, nz, mmin, mmax, nm, 0);
    }

    // calculating upc cross section in YM space
    nucCSYM.resize(ny, std::vector<double>(nm, 0.));
    if (usePolarizedCS) {
      polCSRatio.resize(ny, std::vector<double>(nm));
    }
    nucProcessCS->calcNucCrossSectionYM(nucCSYM, polCSRatio, totCS);

    // setup random samplers
    binEdgesM.resize(nm + 1, 0.);
    binEdgesZ.resize(nz + 1, 0.);
    binEdgesY.resize(ny + 1, 0.);

    for (int i = 0; i < nm + 1; ++i)
      binEdgesM[i] = mmin + dm * i;

    for (int i = 0; i < nz + 1; ++i)
      binEdgesZ[i] = zmin + dz * i;

    for (int i = 0; i < ny + 1; ++i)
      binEdgesY[i] = ymin + dy * i;

    if (usePolarizedCS) {
      for (int i = 0; i < nm; ++i) {
        auto* samplerSZ = new UpcSampler1D(crossSectionZMPolS[i], binEdgesZ, seed);
        auto* samplerPsZ = new UpcSampler1D(crossSectionZMPolPS[i], binEdgesZ, seed);
        samplersCsSZ.push_back(samplerSZ);
        samplersCsPsZ.push_back(samplerPsZ);
      }
    } else {
      for (int i = 0; i < nm; ++i) {
        auto* samplerZ = new UpcSampler1D(crossSectionZM[i], binEdgesZ, seed);
        samplersCsZ.push_back(samplerZ);
      }
    }

    samplerCsYM = new UpcSampler2D(nucCSYM, binEdgesY, binEdgesM, seed);
  } else {
    nucCSYM.resize(ny, std::vector<double>(1, 0.)); // ignoring nm, nz etc.
    nucTargRatioCSYM.resize(ny, std::vector<double>(1, 0.)); // ignoring nm, nz etc.
    nucProcessCS->calcNucCrossSectionY(nucCSYM, nucTargRatioCSYM, totCS);
    binEdgesM.resize(1 + 1, 0.); // force 1 mass bin
    binEdgesM[0] = mmin;
    binEdgesM[1] = mmax;
    binEdgesY.resize(ny + 1, 0.);
    for (int i = 0; i < ny + 1; ++i)
      binEdgesY[i] = ymin + dy * i;
    samplerCsYM = new UpcSampler2D(nucCSYM, binEdgesY, binEdgesM, seed);
  }
}

long int UpcGenerator::generateEvent(std::vector<int>& pdgs,
                                     std::vector<int>& statuses,
                                     std::vector<int>& mothers,
                                     std::vector<TLorentzVector>& particles)
{
  bool isVM = false;
  if (procID == 443 || procID == 100443 || procID == 553)
    isVM = true;

  // clear the particle vectors
  pdgs.clear();
  statuses.clear();
  mothers.clear();
  particles.clear();

  // mass of a particle in final state
  double mPart = nucProcessCS->elemProcess->mPart;
  int partPDG = nucProcessCS->elemProcess->partPDG;
  bool isCharged = nucProcessCS->elemProcess->isCharged;

  // pick pair m and y from nuclear cross section
  double mPair, yPair;
  (*samplerCsYM)(yPair, mPair);
  int yPairBin = samplerCsYM->getBinX(yPair);
  int mPairBin = samplerCsYM->getBinY(mPair);

  // pick z = cos(theta) for corresponding m from elem. cross section
  double cost;
  if (!ignoreCSZ) { // ignoring z distribution for narrow resonances
    if (usePolarizedCS) {
      double frac = polCSRatio[yPairBin][mPairBin];
      bool pickScalar = gRandom->Uniform(0, 1) < frac;
      if (pickScalar) {
        cost = (*samplersCsSZ[mPairBin])();
      } else {
        cost = (*samplersCsPsZ[mPairBin])();
      }
    } else {
      cost = (*samplersCsZ[mPairBin])();
    }
  } else {
    cost = gRandom->Uniform(-1., 1.);
  }

  TLorentzVector pPair;

  if (isVM) {
    double ratio = nucTargRatioCSYM[yPairBin][mPairBin];
    bool target = gRandom->Uniform(0, 1) < ratio;
    nucProcessCS->getMomentumVM(mPair, yPair, target, pPair);
  } else {
    nucProcessCS->getPairMomentum(mPair, yPair, pPair);
  }

  if (isPairProduction) {
    double pMag = std::sqrt(pPair.Mag2() / 4 - mPart * mPart);
    double theta = acos(cost);
    double phi = gRandom->Uniform(0., 2. * M_PI);
    TVector3 vec;
    vec.SetMagThetaPhi(pMag, theta, phi);
    pairProduction(pPair, vec, mPart, partPDG, isCharged, particles, pdgs, mothers, statuses);
  }

  if (isSingleProduction) {
    singleProduction(pPair, partPDG, particles, pdgs, mothers, statuses);
  }

  // check kinematic cuts
  if (!checkKinCuts(particles)) {
    pdgs.clear();
    statuses.clear();
    mothers.clear();
    particles.clear();
    return 0;
  }

  // lepton decays for taus
  // todo: at the moment "fsr" and "decays" flags
  //       are only really meaningful for pythia8
  std::vector<std::vector<int>> daughters;
  if ((doFSR || doDecays) && isPythiaUsed && (procID >= 11 && procID <= 15)) {
    processInPythia(pdgs, statuses, mothers, daughters, particles);
  }

  // uniform pion decays into photon pairs
  if (procID == 111) {
    twoPartDecayUniform(pdgs, statuses, mothers, particles, 1, 0., 22);
    twoPartDecayUniform(pdgs, statuses, mothers, particles, 2, 0., 22);
  }

  // uniform ALP decay into two photons
  if (procID == 51) {
    twoPartDecayUniform(pdgs, statuses, mothers, particles, 1, 0., 22);
  }

  if (isPairProductionVM) {
    twoPartDecayVM(pdgs, statuses, mothers, particles, 1);
  }

  // fill genParticles
  std::vector<int> ds{-1, -1};
  genParticles.clear();
  for (int ii = 0; ii < particles.size(); ++ii) {
    ds = std::vector<int>{-1, -1};
    if (daughters.size() > ii) {
      ds[0] = daughters[ii][0];
      ds[1] = daughters[ii][1];
    }

    TParticle tpart = TParticle(
      pdgs[ii], statuses[ii], mothers[ii], mothers[ii], ds[0], ds[1],
      particles[ii].Px(), particles[ii].Py(),
      particles[ii].Pz(), particles[ii].E(), 0.0, 0.0, 0.0, 0.0);
    genParticles.push_back(tpart);
  }

  return 1;
}

void UpcGenerator::generateEvents()
{
  std::vector<int> pdgs;
  std::vector<int> statuses;
  std::vector<int> mothers;
  std::vector<TLorentzVector> particles;

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

    // generate an event
    if (generateEvent(pdgs, statuses, mothers, particles) == 1) {
      writeEvent(evt, pdgs, statuses, mothers, particles);
      evt++;
    } else {
      rejected++;
    }

    pdgs.clear();
    statuses.clear();
    mothers.clear();
    particles.clear();
  }
  // fiducial cross section
  fidCS = totCS * (double)nEvents / (double)(nEvents + rejected);

  PLOG_INFO << "Event generation is finished!";

  bool doAnyCuts = doPtCut || doEtaCut;
  if (doAnyCuts && nEvents > 0) {
    PLOG_INFO << "Kinematic cuts were used";
    PLOG_INFO << "Number of rejected events = " << rejected;
    PLOG_INFO << std::fixed << std::setprecision(6) << "Cross section with cuts = " << fidCS << " mb";
  }

  if (useROOTOut) {
    if (debug > 0) { // write out cross sections
      auto* hNucCSYM = new TH2D("hNucCSYM", "",
                                binEdgesY.size() - 1, binEdgesY.data(),
                                binEdgesM.size() - 1, binEdgesM.data());
      for (int iy = 0; iy < binEdgesY.size() - 1; ++iy)
        for (int im = 0; im < binEdgesM.size() - 1; ++im)
          hNucCSYM->SetBinContent(iy + 1, im + 1, nucCSYM[iy][im]);
      auto* hNucCSM = hNucCSYM->ProjectionY();
      auto* hNucCSY = hNucCSYM->ProjectionX();
      hNucCSM->Write();
      hNucCSY->Write();
      hNucCSYM->Write();
      delete hNucCSM;
      delete hNucCSY;
      delete hNucCSYM;
    }
    mOutFile->Write();
    mOutFile->Close();

    delete mOutFile;
  }

  if (useHepMCOut) {
    delete writerHepMC;
  }
}
