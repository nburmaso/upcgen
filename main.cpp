//
// created by Nazar Burmasov on 6/25/21.
//

#include "LepGenerator.h"
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

  PLOG_FATAL_IF(argc >= 3) << "Wrong number of parameters! Usage: generate debug_level [optional]number_of_threads";

  // input
  int debugLevel = std::stoi(argv[1]);
  PLOG_INFO << "Debug level: " << debugLevel;
  int numThreads = 1;
  if (argc == 3) {
    numThreads = std::stoi(argv[2]);
  }

  PLOG_INFO << "Initializing the generator...";
  auto* lepGenerator = new LepGenerator();
  lepGenerator->setDebugLevel(debugLevel);
  lepGenerator->setNumThreads(numThreads);

  PLOG_WARNING << "Check inputs:";
  lepGenerator->printParameters();

  PLOG_INFO << "Starting generation process...";
  lepGenerator->generateEvents();

  PLOG_INFO << "Event generation is finished!";

  return 0;
}