/*-------------------------------------------------------------------------------
This file is part of ranger.

Copyright (c) [2014-2018] [Marvin N. Wright]

This software may be modified and distributed under the terms of the MIT license.

Please note that the C++ core of ranger is distributed under MIT license and the
R package "ranger" under GPL3 license.
#-------------------------------------------------------------------------------*/

#ifndef ARGUMENTHANDLER_H_
#define ARGUMENTHANDLER_H_

#include <getopt.h>
#include <string>
#include <vector>

#include "globals.h"

/*
 * Encapsulate getopt.
 * To add an option:
 *    Add to short_options
 *    Add to long_options
 *    Add member variable
 *    Add default value to constructor
 *    Add case in processArguments() function, use try-catch
 *    Add to checkArguments() function?
 *    Add to help function
 *    Add in R version?
 * Access via public members
 */
class ArgumentHandler {
public:
  ArgumentHandler(int argc, char **argv);
  virtual ~ArgumentHandler();

  // Get arguments and catch conversion exceptions
  int processArguments();

  // Check required arguments, ranges, files, ..
  void checkArguments();

  // All command line arguments as member: Capital letters
  std::vector<std::string> alwayssplitvars;
  std::string caseweights;
  std::string depvarname;
  double fraction;
  bool holdout;
  MemoryMode memmode;
  bool savemem;
  std::string predict;
  PredictionType predictiontype;
  uint randomsplits;
  std::string splitweights;
  uint nthreads;
  bool predall;

  // All command line arguments as member: Small letters
  double alpha;
  double minprop;
  std::vector<std::string> catvars;
  std::string file;
  ImportanceMode impmeasure;
  uint targetpartitionsize;
  uint mtry;
  std::string outprefix;
  bool probability;
  SplitRule splitrule;
  std::string statusvarname;
  uint ntree;
  bool replace;
  bool verbose;
  bool write;
  TreeType treetype;
  uint seed;

private:
  // Display messages
  void displayHelp();
  void displayVersion();

  int argc;
  char** argv;

  DISALLOW_COPY_AND_ASSIGN(ArgumentHandler);
};

#endif /* ARGUMENTHANDLER_H_ */
