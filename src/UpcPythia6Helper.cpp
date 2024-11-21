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

// A simple interface to Pythia6

#include "UpcPythia6Helper.h"

#ifdef USE_PYTHIA6

#include "TParticle.h"
#include "TPythia6.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

void UpcPythia6Helper::init()
{
}

void UpcPythia6Helper::process(std::vector<int>& pdgs, std::vector<int>& statuses, std::vector<TLorentzVector>& particles)
{
  for (int i = 0; i < particles.size(); i++) {
    TPythia6::Instance()->Py1ent(0, pdgs[i], particles[i].Energy(), particles[i].Theta(), particles[i].Phi());
    TPythia6::Instance()->GetPrimaries();
    TClonesArray parts("TParticle");
    TPythia6::Instance()->ImportParticles(&parts, "All");
    mPartHolder.emplace_back(parts);
  }
}

// status codes:
// https://www.pythia.org/download/pythia6/pythia6200.pdf
// https://pythia.org/download/pdf/worksheet8153.pdf
int UpcPythia6Helper::import(TClonesArray* particles)
{
  int nParts = 0;
  int shift = 0;
  TClonesArray& clonesParticles = *particles;
  for (int i = 0; i < mPartHolder.size(); i++) {
    for (int j = 0; j < mPartHolder[i].GetEntriesFast(); j++) {
      auto* part = (TParticle*)mPartHolder[i].At(j);
      bool isPrimary = part->GetFirstMother() == 0;
      bool isFinal = part->GetFirstDaughter() == 0;
      int statusP6 = part->GetStatusCode();
      int statusP8 = -1;
      if (statusP6 == 11 && isPrimary) { // primary, decayed
        statusP8 = -23;
      }
      if (statusP6 == 11 && !isFinal && !isPrimary) { // secondary, decayed
        statusP8 = -91;
      }
      if (statusP6 == 1 && isFinal && !isPrimary) { // secondary, "final state"
        statusP8 = 91;
      }
      part->SetStatusCode(statusP8);
      int mother = part->GetFirstMother() == 0 ? -1 : part->GetFirstMother() + shift - 1;
      part->SetFirstMother(mother);
      new (clonesParticles[nParts]) TParticle(*part);
      nParts++;
    }
    shift = nParts;
  }
  mPartHolder.clear();
  return nParts;
}

#endif
