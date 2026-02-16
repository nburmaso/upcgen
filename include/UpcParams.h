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

#pragma once

#include "UpcPhysConstants.h"


// todo: refactor generator parameters?

class UpcParams
{
 public:
  UpcParams() = default;
  virtual ~UpcParams() = 0;

  // Woods-Saxon parameters
  inline static double rho0{0.}; // fm^-3
  inline static double R{6.68};  // fm
  inline static double a{0.447}; // fm

  // parameters of the nucleus
  inline static double Z{82};
  inline static double A{208};
  inline static double mNucl{(Z * phc::mProt + (A - Z) * phc::mNeut) / A};

  // beam parameters
  inline static double sqrts{5020.};
  inline static double g1{sqrts / (2. * phc::mProt)};
  inline static double g2{sqrts / (2. * phc::mProt)};

  // Gaussian integration n = 10
  // since cos is symmetric around 0 we only need 5
  // of the points in the gaussian integration.
  static const int ngi10 = 10;
  double weights10[ngi10]{0.0666713443086881,
                          0.1494513491505806,
                          0.2190863625159820,
                          0.2692667193099963,
                          0.2955242247147529,
                          0.2955242247147529,
                          0.2692667193099963,
                          0.2190863625159820,
                          0.1494513491505806,
                          0.0666713443086881};

  double abscissas10[ngi10]{-0.9739065285171717,
                            -0.8650633666889845,
                            -0.6794095682990244,
                            -0.4333953941292472,
                            -0.1488743389816312,
                            0.1488743389816312,
                            0.4333953941292472,
                            0.6794095682990244,
                            0.8650633666889845,
                            0.9739065285171717};
};