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

#ifndef UPCGENERATOR_SRC_UPCPYTHIA6HELPER_H_
#define UPCGENERATOR_SRC_UPCPYTHIA6HELPER_H_

#include "UpcPythiaBase.h"

#ifdef USE_PYTHIA6

#include "TString.h"
#include "TArrayF.h"

class UpcPythia6Helper : public UpcPythiaBase
{
 public:
  UpcPythia6Helper() = default;
  ~UpcPythia6Helper() override = default;

  void init() override;

  void process(std::vector<int>& pdgs, std::vector<int>& statuses, std::vector<TLorentzVector>& particles) override;

  int import(TClonesArray* particles) override;

  std::vector<TClonesArray> mPartHolder;
};

#endif

#endif // UPCGENERATOR_SRC_UPCPYTHIA6HELPER_H_
