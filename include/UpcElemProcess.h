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

/// base (abstract) class for elementary processes

#ifndef UPCGENERATOR_INCLUDE_UPCELEMPROCESS_H_
#define UPCGENERATOR_INCLUDE_UPCELEMPROCESS_H_

class UpcElemProcess
{
 public:
  UpcElemProcess() = default;
  virtual ~UpcElemProcess() {}

  // mass of a particle in final state
  double mPart{};

  // pdg code of a particle
  int partPDG{};

  // final-state is charged or not
  bool isCharged{};

  // `standard` unpolarized cross sections
  virtual double calcCrossSectionM(double m) = 0;

  virtual double calcCrossSectionZM(double z, double m) = 0;

  // polarized cross sections
  // scalar part
  virtual double calcCrossSectionMPolS(double m) = 0;

  virtual double calcCrossSectionZMPolS(double z, double m) = 0;

  // pseudoscalar part
  virtual double calcCrossSectionMPolPS(double m) = 0;

  virtual double calcCrossSectionZMPolPS(double z, double m) = 0;
};

#endif // UPCGENERATOR_INCLUDE_UPCELEMPROCESS_H_
