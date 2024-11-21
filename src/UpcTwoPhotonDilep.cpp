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

#include "UpcTwoPhotonDilep.h"
#include "UpcPhysConstants.h"

#include <cmath>

using namespace std;

UpcTwoPhotonDilep::UpcTwoPhotonDilep(int partPDG)
{
  this->partPDG = partPDG;
  isCharged = true;
  // populating lepton mass map
  // data from PDG
  if (partPDG == 11) { // electron
    mPart = 0.000510998946;
  }
  if (partPDG == 13) { // muon
    mPart = 0.1056583745;
  }
  if (partPDG == 15) { // tau
    mPart = 1.77686;
  }
}

double UpcTwoPhotonDilep::calcCrossSectionM(double m)
{
  double s = m * m;                 // cms invariant mass squared
  double x = 4 * mPart * mPart / s; // inverse lepton gamma-factor squared = 1/g^2 in cms
  double b = sqrt(1 - x);           // lepton relativistic velocity in cms
  double y = atanh(b);              // lepton rapidity in cms
  double cs = 0;
  cs += (2 + 2 * x - x * x) * y - b * (1 + x);
  cs += 4 * y * aLep;
  cs += (4 * b / x + y) * aLep * aLep;
  cs += (4 * b / x - 2 * y) * aLep * aLep * aLep;
  cs += ((7. / 12.) * b / x + (1. / 6.) * b / x / x - 0.5 * y) * aLep * aLep * aLep * aLep;
  cs *= 4 * phys_consts::hc * phys_consts::hc * 1e7 * phys_consts::alpha * phys_consts::alpha * M_PI / s;
  return cs;
  // [nb] //
}

double UpcTwoPhotonDilep::calcCrossSectionZM(double z, double m)
{
  double s = m * m;
  double k = sqrt(s) / 2.;                // photon/lepton energy in cm system in GeV
  double p = sqrt(k * k - mPart * mPart); // outgoing lepton momentum in GeV
  double norm = 2 * M_PI * phys_consts::alpha * phys_consts::alpha / s * p / k;
  double kt = -2 * k * (k - z * p) / mPart / mPart;
  double ku = -2 * k * (k + z * p) / mPart / mPart;
  double ks = kt + ku;
  double kp = kt * ku;
  double kq = 1. / kt + 1. / ku;
  double kr = ku / kt + kt / ku;
  double cs = 0;
  cs += -8. * (4. * kq * kq + 4. * kq - kr);
  cs += 16. * (2. + kr) * aLep;
  cs += 4. * (2. - 4. * ks + kr) * aLep * aLep;
  cs += -8. * (2. + 2. * ks + kr) * aLep * aLep * aLep;
  cs += -4. * (4. + 2. * ks + 2. * kr - kp) * aLep * aLep * aLep * aLep;
  cs *= norm;
  return cs; // [GeV^-2]
}

double UpcTwoPhotonDilep::calcCrossSectionMPolS(double m)
{
  double r = 2 * mPart / m;
  double r2 = r * r;
  if (r > 1) {
    return 0;
  }
  double cs_s = 4 * M_PI * phys_consts::alpha * phys_consts::alpha * phys_consts::hc * phys_consts::hc / m / m * ((1 + r * r - 3. / 4. * r * r * r * r) * 2 * log(1 / r + sqrt(1 / r / r - 1)) - (1 + 3. / 2. * r * r) * sqrt(1 - r * r));
  return cs_s;
  // fm^2 //
}

double UpcTwoPhotonDilep::calcCrossSectionZMPolS(double z, double m)
{
  double mLep2 = mPart * mPart;
  double m2 = m * m;
  double z2 = z * z;
  double cs_s = 2 * M_PI * phys_consts::alpha * phys_consts::alpha;
  cs_s *= m2 - 4 * mLep2;
  cs_s *= sqrt(m2 - 4 * mLep2);
  cs_s *= (4 * mLep2 * (3 - 2 * z2 + z2 * z2)) + m2 * (1 - z2 * z2);
  cs_s /= m2 * m * (m2 * (1 - z2) + 4 * mLep2 * z2) * (m2 * (1 - z2) + 4 * mLep2 * z2);
  return cs_s;
  // GeV^-2 //
}

double UpcTwoPhotonDilep::calcCrossSectionMPolPS(double m)
{
  double r = 2 * mPart / m;
  double r2 = r * r;
  double r4 = r2 * r2;
  if (r > 1) {
    return 0;
  }
  return 4 * M_PI * phys_consts::alpha * phys_consts::alpha * phys_consts::hc * phys_consts::hc / m / m * ((1 + r * r - 1. / 4. * r * r * r * r) * 2 * log(1 / r + sqrt(1 / r / r - 1)) - (1 + 1. / 2. * r * r) * sqrt(1 - r * r));
  // fm^2 //
}

double UpcTwoPhotonDilep::calcCrossSectionZMPolPS(double z, double m)
{
  double mLep2 = mPart * mPart;
  double m2 = m * m;
  double z2 = z * z;
  double cs_s = 2 * M_PI * phys_consts::alpha * phys_consts::alpha;
  cs_s *= sqrt(m2 - 4 * mLep2);
  cs_s *= m2 * m2 * (1 - z2 * z2) + 8 * m2 * mLep2 * (1 - z2 + z2 * z2) - 16 * mLep2 * mLep2 * (1 - z2) * (1 - z2);
  cs_s /= m2 * m * (m2 * (1 - z2) + 4 * mLep2 * z2) * (m2 * (1 - z2) + 4 * mLep2 * z2);
  return cs_s;
  // GeV^-2 //
}
