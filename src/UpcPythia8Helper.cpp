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

// A simple interface to Pythia8

#include "UpcPythia8Helper.h"
#include "plog/Appenders/ColorConsoleAppender.h"
#include "plog/Formatters/TxtFormatter.h"
#include "plog/Init.h"
#include "plog/Initializers/RollingFileInitializer.h"
#include "plog/Log.h"

// can be used only if built with Pythia8
#ifdef USE_PYTHIA8

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TParticle.h"

UpcPythia8Helper::UpcPythia8Helper() : mPythia8(new Pythia8::Pythia())
{
}

void UpcPythia8Helper::init()
{
  if (mDoFSR) {
    mPythia8->readString("ProcessLevel:all = off");
    mPythia8->readString("PartonShowers:Model = 1");
  }
  if (mSeed == 0) {
    mPythia8->readString("Random:setSeed = on");
    mPythia8->readString("Random:seed = 0");
  } else {
    mPythia8->readString("Random:setSeed = off");
    mPythia8->readString(Form("Random:seed = %ld", mSeed));
  }
  mPythia8->readString("SoftQCD:elastic = on");
  mPythia8->init();
}

void UpcPythia8Helper::process(std::vector<int>& pdgs,
                               std::vector<int>& statuses,
                               std::vector<TLorentzVector>& particles)
{
  clearEvent();
  for (int i = 0; i < particles.size(); i++) {
    appendParticle(pdgs[i], statuses[i], &particles[i]);
    int idPart = mPythia8->event[i].id();
    if (statuses[i] == 23) {
      mPythia8->particleData.mayDecay(idPart, true);
    }
  }
  if (mDoFSR) {
    mPythia8->getShowerModelPtr()->getTimeShower()->showerQED(0, 1, mPythia8->event, particles[0].Pt());
  }
  if (mDoDecays) {
    mPythia8->moreDecays();
  }
}

int UpcPythia8Helper::import(TClonesArray* particles)
{
  if (particles == nullptr) {
    return 0;
  }
  TClonesArray& clonesParticles = *particles;
  clonesParticles.Clear();
  int nparts = 0;
  int ioff = 0;
  if (mPythia8->event[0].id() == 90) {
    ioff = -1;
  }

  for (const auto& part : mPythia8->event) {
    if (part.id() == 90) {
      continue;
    }
    new (clonesParticles[nparts]) TParticle(
      part.id(),
      part.status(),
      part.mother1() + ioff,
      part.mother2() + ioff,
      part.daughter1() + ioff,
      part.daughter2() + ioff,
      part.px(),
      part.py(),
      part.pz(),
      part.e(),
      part.xProd(),
      part.yProd(),
      part.zProd(),
      part.tProd());
    nparts++;
  }
  return nparts;
}

void UpcPythia8Helper::appendParticle(int pdg, int status, TLorentzVector* p)
{
  mPythia8->event.append(pdg, status, 0, 0, p->Px(), p->Py(), p->Pz(), p->E(), p->M(), p->Pt());
}

void UpcPythia8Helper::clearEvent()
{
  mPythia8->event.clear();
}

#endif