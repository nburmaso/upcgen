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

/// base (abstract) class for elementary processes

#pragma once

class UpcElemProcess
{
 public:
  UpcElemProcess() = default;
  virtual ~UpcElemProcess() {}

  // mass of the produced particle
  double mPart{};

  // mass of decay particle
  double mDght{};
  int dghtPDG{};

  // pdg code of a particle
  int partPDG{};

  // final-state is charged or not
  bool isCharged{};

  // `standard` unpolarized cross sections
  virtual double calcCrossSectionY(double y) { return 0.; }

  virtual double calcCrossSectionM(double m) { return 0.; }

  virtual double calcCrossSectionZM(double z, double m) = 0;

  // polarized cross sections
  // scalar part
  virtual double calcCrossSectionMPolS(double m) = 0;

  virtual double calcCrossSectionZMPolS(double z, double m) = 0;

  // pseudoscalar part
  virtual double calcCrossSectionMPolPS(double m) = 0;

  virtual double calcCrossSectionZMPolPS(double z, double m) = 0;
};
