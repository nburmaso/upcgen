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

// Base class for pythia interfaces

#ifndef UPCGENERATOR_INCLUDE_UPCPYTHIABASE_H_
#define UPCGENERATOR_INCLUDE_UPCPYTHIABASE_H_

#include "TObject.h"
class TClonesArray;
class TLorentzVector;

class UpcPythiaBase
{
 public:
  UpcPythiaBase() { ; }
  virtual ~UpcPythiaBase() { ; }

  virtual void init() = 0;

  virtual void process(std::vector<int>& pdgs, std::vector<int>& statuses, std::vector<TLorentzVector>& particles) = 0;

  virtual int import(TClonesArray* particles) = 0;

  // auxiliary methods
  void setFSR(bool doFSR) { mDoFSR = doFSR; }
  void setSeed(long seed) { mSeed = seed; }
  void setDecays(bool doDecays) { mDoDecays = doDecays; }

 protected:
  long mSeed{0};
  bool mDoFSR{false};
  bool mDoDecays{false};
};

#endif // UPCGENERATOR_INCLUDE_UPCPYTHIABASE_H_
