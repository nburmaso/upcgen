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

// A simple interface to Pythia8

#include "UpcPythia8Helper.h"
#include "plog/Appenders/ColorConsoleAppender.h"
#include "plog/Formatters/TxtFormatter.h"
#include "plog/Init.h"
#include "plog/Initializers/RollingFileInitializer.h"
#include "plog/Log.h"

// can be used only if built with Pythia8
#ifdef USE_PYTHIA8

#include "TLorentzVector.h"
#include "TPythia8.h"

UpcPythia8Helper::UpcPythia8Helper() : fPythia8(new TPythia8()),
                                       fDebug(0)
{
}

void UpcPythia8Helper::init()
{
  if (fDoFSR) {
    fPythia8->Pythia8()->readString("ProcessLevel:all = off");
    fPythia8->Pythia8()->readString("SpaceShower:QEDshowerByL = on");
    fPythia8->Pythia8()->readString("TimeShower:QEDshowerByL = on");
  }
  fPythia8->Pythia8()->readString("PartonLevel:FSR = on");
  fPythia8->Pythia8()->readString("SoftQCD:elastic = on");
  fPythia8->Pythia8()->init();
}

void UpcPythia8Helper::decay(std::vector<int>& pdgs,
                             std::vector<TLorentzVector>& particles)
{
  clearEvent();
  for (int i = 0; i < particles.size(); i++) {
    appendParticle(pdgs[i], &particles[i]);
    int idPart = fPythia8->Pythia8()->event[i].id();
    // fPythia8->Pythia8()->particleData.mayDecay(idPart, true);
  }
  if (fDoFSR) {
    fPythia8->Pythia8()->forceTimeShower(0, 0, particles[0].Pt(), 1);
  }
  // fPythia8->Pythia8()->moreDecays();
  if (fDebug > 0)
    fPythia8->EventListing();
}

int UpcPythia8Helper::import(TClonesArray* particles)
{
  return (fPythia8->ImportParticles(particles, "All"));
}

void UpcPythia8Helper::appendParticle(int pdg, TLorentzVector* p)
{
  fPythia8->Pythia8()->event.append(pdg, 11, 0, 0, p->Px(), p->Py(), p->Pz(), p->E(), p->M(), p->Pt());
}

void UpcPythia8Helper::clearEvent()
{
  fPythia8->Pythia8()->event.clear();
}

#endif