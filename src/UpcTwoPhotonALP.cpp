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

#include "UpcTwoPhotonALP.h"
#include "UpcPhysConstants.h"

#include <cmath>

// cross section in narrow resonance approximation
double UpcTwoPhotonALP::calcCrossSectionM(double m)
{
  double cs = 4. * M_PI * M_PI * width / (mPart * mPart);
  cs *= phys_consts::hc * phys_consts::hc* 1e7 * phys_consts::alpha * phys_consts::alpha;
  return cs;
}
