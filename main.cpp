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

#include "UpcGenerator.h"

#include <plog/Log.h>
#include <plog/Init.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Appenders/ColorConsoleAppender.h>

#include <algorithm>

class InputParser
{
 public:
  InputParser(int& argc, char** argv)
  {
    for (int i = 1; i < argc; ++i)
      this->tokens.emplace_back(argv[i]);
  }

  const std::string& getCmdOption(const std::string& option) const
  {
    std::vector<std::string>::const_iterator itr;
    itr = std::find(this->tokens.begin(), this->tokens.end(), option);
    if (itr != this->tokens.end() && ++itr != this->tokens.end()) {
      return *itr;
    }
    static const std::string empty_string("");
    return empty_string;
  }

  bool cmdOptionExists(const std::string& option) const
  {
    return std::find(this->tokens.begin(), this->tokens.end(), option) != this->tokens.end();
  }

 private:
  std::vector<std::string> tokens;
};

int main(int argc, char** argv)
{
  InputParser input(argc, argv);

  // define default values of command line parameters
  // logging levels:
  //  0 -> no debug logs
  //  1 -> basic debug logs
  //  2 -> verbose logs with intermediate calculation results
  int debugLevel = 0;
  int numThreads = 1;
  std::string parFileName("parameters.in");

  // initialize logger
  static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
  plog::init(plog::verbose, &consoleAppender);

  if (input.cmdOptionExists("-h")) {
    printf("Available options:\n"
           "-h        -- help message\n"
           "-debug    -- set debug level: 0=no debug messages (default),\n"
           "                              1=save calculated cross sections into 'events.root' (for ROOT output only),\n"
           "                              2=print out intermediate calculation results and events info (warning, lots of messages)\n"
           "-nthreads -- set number of threads for cross section and luminosity calculation (default is 1)\n"
           "-parfile -- name of input parameter file\n");
    std::_Exit(0);
  }

  const auto debugLevelOpt = input.getCmdOption("-debug");
  if (!debugLevelOpt.empty()) {
    debugLevel = std::stoi(debugLevelOpt);
  }
  const auto numThreadsOpt = input.getCmdOption("-nthreads");
  if (!numThreadsOpt.empty()) {
    numThreads = std::stoi(numThreadsOpt);
  }
  const auto parFileNameOpt = input.getCmdOption("-parfile");
  if (!parFileNameOpt.empty()) {
    parFileName = std::string(parFileNameOpt);
  }

  PLOG_INFO << "Configuring the generator...";
  auto* upcGenerator = new UpcGenerator();
  upcGenerator->setDebugLevel(debugLevel);
  upcGenerator->setNumThreads(numThreads);
  upcGenerator->setParFile(parFileName);
  upcGenerator->configGeneratorFromFile();

  PLOG_INFO << "Initializing the generator...";
  upcGenerator->init();

  PLOG_INFO << "Starting generation process...";
  upcGenerator->generateEvents();

  delete upcGenerator;

  return 0;
}
