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

#include "UpcGenerator.h"
#include <plog/Log.h>
#include <plog/Init.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Appenders/ColorConsoleAppender.h>

int main(int argc, char** argv)
{
  // logging levels:
  //  0 -> no debug logs
  //  1 -> basic debug logs
  //  2 -> verbose logs with intermediate calculation results

  // initialize logger
  static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
  plog::init(plog::verbose, &consoleAppender);

  PLOG_FATAL_IF(argc > 3) << "Wrong number of parameters! Usage: generate debug_level [optional]number_of_threads";
  if (argc > 3) {
    return -1;
  }

  // input
  int debugLevel = std::stoi(argv[1]);
  PLOG_INFO << "Debug level: " << debugLevel;
  int numThreads = 1;
  if (argc == 3) {
    numThreads = std::stoi(argv[2]);
  }

  PLOG_INFO << "Initializing the generator...";
  auto* upcGenerator = new UpcGenerator();
  upcGenerator->setDebugLevel(debugLevel);
  upcGenerator->setNumThreads(numThreads);

  PLOG_WARNING << "Check inputs:";
  upcGenerator->printParameters();

  PLOG_INFO << "Starting generation process...";
  upcGenerator->generateEvents();

  PLOG_INFO << "Event generation is finished!";

  return 0;
}