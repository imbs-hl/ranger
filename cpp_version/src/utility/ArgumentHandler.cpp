/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#include <fstream>
#include <iostream>
#include <stdexcept>

#include "ArgumentHandler.h"
#include "version.h"
#include "utility.h"

namespace ranger {

ArgumentHandler::ArgumentHandler(int argc, char **argv) :
    caseweights(""), depvarname(""), fraction(0), holdout(false), memmode(MEM_DOUBLE), savemem(false), skipoob(false), predict(
        ""), predictiontype(DEFAULT_PREDICTIONTYPE), randomsplits(DEFAULT_NUM_RANDOM_SPLITS), splitweights(""), nthreads(
        DEFAULT_NUM_THREADS), predall(false), alpha(DEFAULT_ALPHA), minprop(DEFAULT_MINPROP), maxdepth(
        DEFAULT_MAXDEPTH), file(""), impmeasure(DEFAULT_IMPORTANCE_MODE), targetpartitionsize(0), minbucket(0), mtry(0), outprefix(
        "ranger_out"), probability(false), splitrule(DEFAULT_SPLITRULE), statusvarname(""), ntree(DEFAULT_NUM_TREE), replace(
        true), verbose(false), write(false), treetype(TREE_CLASSIFICATION), seed(0), usedepth(false) {
  this->argc = argc;
  this->argv = argv;
}

int ArgumentHandler::processArguments() {

  // short options
  char const *short_options = "A:C:D:F:HM:NOP:Q:R:S:U:XZa:b:c:d:f:hi:j:kl:m:n:o:pr:s:t:uvwy:z:";

// long options: longname, no/optional/required argument?, flag(not used!), shortname
    const struct option long_options[] = {

      { "alwayssplitvars",      required_argument,  0, 'A'},
      { "caseweights",          required_argument,  0, 'C'},
      { "depvarname",           required_argument,  0, 'D'},
      { "fraction",             required_argument,  0, 'F'},
      { "holdout",              no_argument,        0, 'H'},
      { "memmode",              required_argument,  0, 'M'},
      { "savemem",              no_argument,        0, 'N'},
      { "skipoob",              no_argument,        0, 'O'},
      { "predict",              required_argument,  0, 'P'},
      { "predictiontype",       required_argument,  0, 'Q'},
      { "randomsplits",         required_argument,  0, 'R'},
      { "splitweights",         required_argument,  0, 'S'},
      { "nthreads",             required_argument,  0, 'U'},
      { "predall",              no_argument,        0, 'X'},
      { "version",              no_argument,        0, 'Z'},
      { "alpha",                required_argument,  0, 'a'},
      { "minprop",              required_argument,  0, 'b'},
      { "catvars",              required_argument,  0, 'c'},
      { "maxdepth",             required_argument,  0, 'd'},
      { "file",                 required_argument,  0, 'f'},
      { "help",                 no_argument,        0, 'h'},
      { "impmeasure",           required_argument,  0, 'i'},
      { "regcoef",              required_argument,  0, 'j'},
      { "usedepth",             no_argument,        0, 'k'},
      { "targetpartitionsize",  required_argument,  0, 'l'},
      { "mtry",                 required_argument,  0, 'm'},
      { "minbucket",            required_argument,  0, 'n'},
      { "outprefix",            required_argument,  0, 'o'},
      { "probability",          no_argument,        0, 'p'},
      { "splitrule",            required_argument,  0, 'r'},
      { "statusvarname",        required_argument,  0, 's'},
      { "ntree",                required_argument,  0, 't'},
      { "noreplace",            no_argument,        0, 'u'},
      { "verbose",              no_argument,        0, 'v'},
      { "write",                no_argument,        0, 'w'},
      { "treetype",             required_argument,  0, 'y'},
      { "seed",                 required_argument,  0, 'z'},

      { 0, 0, 0, 0}
    };

  while (1) {
    int option_index = 0;
    int c = getopt_long(argc, argv, short_options, long_options, &option_index);

    // stop if no options left
    if (c == -1) {
      break;
    }

    switch (c) {

    // upper case options
    case 'A':
      splitString(alwayssplitvars, optarg, ',');
      break;

    case 'C':
      caseweights = optarg;
      break;

    case 'D':
      depvarname = optarg;
      break;

    case 'F':
      try {
        fraction = std::stod(optarg);
        if (fraction > 1 || fraction <= 0) {
          throw std::runtime_error("");
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'fraction'. Please give a value in (0,1]. See '--help' for details.");
      }
      break;

    case 'H':
      holdout = true;
      break;

    case 'M':
      try {
        memmode = (MemoryMode) std::stoi(optarg);
        if (memmode > MAX_MEM_MODE) {
          throw std::runtime_error("");
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'memmode'. Please give a positive integer. See '--help' for details.");
      }
      break;

    case 'N':
      savemem = true;
      break;

    case 'O':
      skipoob = true;
      break;

    case 'P':
      predict = optarg;
      break;

    case 'Q':
      try {
        switch (std::stoi(optarg)) {
        case 1:
          predictiontype = RESPONSE;
          break;
        case 2:
          predictiontype = TERMINALNODES;
          break;
        default:
          throw std::runtime_error("");
          break;
        }
      } catch (...) {
        throw std::runtime_error("Illegal prediction type selected. See '--help' for details.");
      }
      break;

    case 'R':
      try {
        int temp = std::stoi(optarg);
        if (temp < 1) {
          throw std::runtime_error("");
        } else {
          randomsplits = temp;
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'randomsplits'. Please give a positive integer. See '--help' for details.");
      }
      break;

    case 'S':
      splitweights = optarg;
      break;

    case 'U':
      try {
        int temp = std::stoi(optarg);
        if (temp < 1) {
          throw std::runtime_error("");
        } else {
          nthreads = temp;
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'nthreads'. Please give a positive integer. See '--help' for details.");
      }
      break;

    case 'X':
      predall = true;
      break;

    case 'Z':
      displayVersion();
      return -1;
      break;

      // lower case options
    case 'a':
      try {
        double temp = std::stod(optarg);
        if (temp < 0 || temp > 1) {
          throw std::runtime_error("");
        } else {
          alpha = temp;
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'alpha'. Please give a value between 0 and 1. See '--help' for details.");
      }
      break;

    case 'b':
      try {
        double temp = std::stod(optarg);
        if (temp < 0 || temp > 0.5) {
          throw std::runtime_error("");
        } else {
          minprop = temp;
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'minprop'. Please give a value between 0 and 0.5. See '--help' for details.");
      }
      break;

    case 'c':
      splitString(catvars, optarg, ',');
      break;

    case 'd':
      try {
        int temp = std::stoi(optarg);
        if (temp < 0) {
          throw std::runtime_error("");
        } else {
          maxdepth = temp;
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'maxdepth'. Please give a positive integer. See '--help' for details.");
      }
      break;

    case 'f':
      file = optarg;
      break;

    case 'h':
      displayHelp();
      return -1;
      break;

    case 'i':
      try {
        impmeasure = (ImportanceMode) std::stoi(optarg);
        if (impmeasure > MAX_IMP_MODE) {
          throw std::runtime_error("");
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'impmeasure'. Please give a positive integer. See '--help' for details.");
      }
      break;

    case 'j':
      splitString(regcoef, optarg, ',');
      break;

    case 'k':
      usedepth = true;
      break;

    case 'l':
      try {
        int temp = std::stoi(optarg);
        if (temp < 1) {
          throw std::runtime_error("");
        } else {
          targetpartitionsize = temp;
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'targetpartitionsize'. Please give a positive integer. See '--help' for details.");
      }
      break;

    case 'm':
      try {
        int temp = std::stoi(optarg);
        if (temp < 1) {
          throw std::runtime_error("");
        } else {
          mtry = temp;
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'mtry'. Please give a positive integer. See '--help' for details.");
      }
      break;

    case 'n':
      try {
        int temp = std::stoi(optarg);
        if (temp < 1) {
          throw std::runtime_error("");
        } else {
          minbucket = temp;
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'minbucket'. Please give a positive integer. See '--help' for details.");
      }
      break;
      
    case 'o':
      outprefix = optarg;
      break;

    case 'p':
      probability = true;
      break;

    case 'r':
      try {
        switch (std::stoi(optarg)) {
        case 1:
          splitrule = LOGRANK;
          break;
        case 2:
          splitrule = AUC;
          break;
        case 3:
          splitrule = AUC_IGNORE_TIES;
          break;
        case 4:
          splitrule = MAXSTAT;
          break;
        case 5:
          splitrule = EXTRATREES;
          break;
        case 6:
          splitrule = BETA;
          break;
        case 7:
          splitrule = HELLINGER;
          break;
        default:
          throw std::runtime_error("");
          break;
        }
      } catch (...) {
        throw std::runtime_error("Illegal splitrule selected. See '--help' for details.");
      }
      break;

    case 's':
      statusvarname = optarg;
      break;

    case 't':
      try {
        int temp = std::stoi(optarg);
        if (temp < 1) {
          throw std::runtime_error("");
        } else {
          ntree = temp;
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'ntree'. Please give a positive integer. See '--help' for details.");
      }
      break;

    case 'u':
      replace = false;
      break;

    case 'v':
      verbose = true;
      break;

    case 'w':
      write = true;
      break;

    case 'y':
      try {
        switch (std::stoi(optarg)) {
        case 1:
          treetype = TREE_CLASSIFICATION;
          break;
        case 3:
          treetype = TREE_REGRESSION;
          break;
        case 5:
          treetype = TREE_SURVIVAL;
          break;
        default:
          throw std::runtime_error("");
          break;
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'treetype'. Please give a positive integer. See '--help' for details.");
      }
      break;

    case 'z':
      try {
        int temp = std::stoi(optarg);
        if (temp < 0) {
          throw std::runtime_error("");
        } else {
          seed = temp;
        }
      } catch (...) {
        throw std::runtime_error(
            "Illegal argument for option 'seed'. Please give a positive integer. See '--help' for details.");
      }
      break;

    default:
      break;

    }
  }

  // print all other parameters
  while (optind < argc) {
    std::cout << "Other parameter, not processed: " << argv[optind++] << std::endl;
  }

  return 0;
}

void ArgumentHandler::checkArguments() {

  // required arguments
  if (file.empty()) {
    throw std::runtime_error("Please specify an input filename with '--file'. See '--help' for details.");
  }
  if (predict.empty() && depvarname.empty()) {
    throw std::runtime_error("Please specify a dependent variable name with '--depvarname'. See '--help' for details.");
  }

  if (predict.empty() && treetype == TREE_SURVIVAL && statusvarname.empty()) {
    throw std::runtime_error("Please specify a status variable name with '--statusvarname'. See '--help' for details.");
  }
  if (treetype != TREE_SURVIVAL && !statusvarname.empty()) {
    throw std::runtime_error("Option '--statusvarname' only applicable for survival forest. See '--help' for details.");
  }

  if (treetype == TREE_SURVIVAL && splitrule == MAXSTAT && impmeasure == IMP_GINI) {
    throw std::runtime_error(
        "Node impurity variable importance not supported for survival forests with MAXSTAT splitrule. See '--help' for details.");
  }

  if (treetype != TREE_CLASSIFICATION && probability) {
    throw std::runtime_error("Probability estimation is only applicable to classification forests.");
  }

  // Get treetype for prediction
  if (!predict.empty()) {
    std::ifstream infile;
    infile.open(predict, std::ios::binary);
    if (!infile.good()) {
      throw std::runtime_error("Could not read from input file: " + predict + ".");
    }

    // Do not read dependent variable names
    uint num_dependent_variables;
    infile.read((char*) &num_dependent_variables, sizeof(num_dependent_variables));
    for (size_t i = 0; i < num_dependent_variables; ++i) {
      size_t length;
      infile.read((char*) &length, sizeof(size_t));
      infile.ignore(length);
    }

    // Do not read num_trees
    infile.ignore(sizeof(size_t));

    // Do not read is_ordered_variable
    size_t length;
    infile.read((char*) &length, sizeof(length));
    infile.ignore(length * sizeof(bool));

    // Do not read number of variables
    infile.ignore(sizeof(size_t));

    // Get treetype
    infile.read((char*) &treetype, sizeof(treetype));
    infile.close();
  }

  if (predict.empty() && predall) {
    throw std::runtime_error("Option '--predall' only available in prediction mode.");
  }

  if (!alwayssplitvars.empty() && !splitweights.empty()) {
    throw std::runtime_error("Please use only one option of splitweights and alwayssplitvars.");
  }

  // Check splitrule
  if (((splitrule == AUC || splitrule == AUC_IGNORE_TIES) && treetype != TREE_SURVIVAL)
      || (splitrule == MAXSTAT && (treetype != TREE_SURVIVAL && treetype != TREE_REGRESSION))
      || (splitrule == BETA && treetype != TREE_REGRESSION)
      || (splitrule == HELLINGER && treetype != TREE_CLASSIFICATION && treetype != TREE_PROBABILITY)) {
    throw std::runtime_error("Illegal splitrule selected. See '--help' for details.");
  }

  // Check holdout mode
  if (holdout && caseweights.empty()) {
    throw std::runtime_error("Case weights required to use holdout mode.");
  }

  // Unordered survival splitting only available for logrank or extratrees splitrule
  if (treetype == TREE_SURVIVAL && !catvars.empty() && (splitrule != LOGRANK && splitrule != EXTRATREES)) {
    throw std::runtime_error("Unordered splitting in survival trees only available for LOGRANK splitrule.");
  }

  // Memory save option not allowed in unordered extratrees mode
  if (splitrule == EXTRATREES && !catvars.empty() && savemem) {
    throw std::runtime_error("savemem option not possible in extraTrees mode with unordered predictors.");
  }

  // Corrected impurity importance not allowed if split weights used
  if (!splitweights.empty() && impmeasure == IMP_GINI_CORRECTED) {
    throw std::runtime_error("Corrected impurity importance not supported in combination with splitweights.");
  }

  // Regularization coefficient
  if (regcoef.size() > 0) {
    for (auto& r : regcoef) {
      if (r > 1) {
        throw std::runtime_error("The regularization coefficients cannot be greater than 1.");
      }
      if (r <= 0) {
        throw std::runtime_error("The regularization coefficients must be positive.");
      }
    }

    if (nthreads != 1) {
      std::cout << "Warning: Parallelization deactivated (regularization used)." << std::endl;
      nthreads = 1;
    }
  }
}

void ArgumentHandler::displayHelp() {
  std::cout << "Usage: " << std::endl;
  std::cout << "    " << argv[0] << " [options]" << std::endl;
  std::cout << std::endl;

  std::cout << "Options:" << std::endl;
  std::cout << "    " << "--help                        Print this help." << std::endl;
  std::cout << "    " << "--version                     Print version and citation information." << std::endl;
  std::cout << "    " << "--verbose                     Turn on verbose mode." << std::endl;
  std::cout << "    " << "--file FILE                   Filename of input data. Only numerical values are supported."
      << std::endl;
  std::cout << "    " << "--treetype TYPE               Set tree type to:" << std::endl;
  std::cout << "    " << "                              TYPE = 1: Classification." << std::endl;
  std::cout << "    " << "                              TYPE = 3: Regression." << std::endl;
  std::cout << "    " << "                              TYPE = 5: Survival." << std::endl;
  std::cout << "    " << "                              (Default: 1)" << std::endl;
  std::cout << "    "
      << "--probability                 Grow a Classification forest with probability estimation for the classes."
      << std::endl;
  std::cout << "    " << "                              Use in combination with --treetype 1." << std::endl;
  std::cout << "    "
      << "--depvarname NAME             Name of dependent variable. For survival trees this is the time variable."
      << std::endl;
  std::cout << "    " << "--statusvarname NAME          Name of status variable, only applicable for survival trees."
      << std::endl;
  std::cout << "    " << "                              Coding is 1 for event and 0 for censored." << std::endl;
  std::cout << "    " << "--ntree N                     Set number of trees to N." << std::endl;
  std::cout << "    " << "                              (Default: 500)" << std::endl;
  std::cout << "    " << "--mtry N                      Number of variables to possibly split at in each node."
      << std::endl;
  std::cout << "    " << "                              (Default: sqrt(p) with p = number of independent variables)"
      << std::endl;
  std::cout << "    " << "--targetpartitionsize N       Set minimal node size to N." << std::endl;
  std::cout << "    " << "--minbucket N                 Set min bucket size to N." << std::endl;
  std::cout << "    "
      << "                              For Classification and Regression growing is stopped if a node reaches a size smaller than N."
      << std::endl;
  std::cout << "    "
      << "                              For Survival growing is stopped if one child would reach a size smaller than N."
      << std::endl;
  std::cout << "    "
      << "                              This means nodes with size smaller N can occur for Classification and Regression."
      << std::endl;
  std::cout << "    "
      << "                              (Default: 1 for Classification, 5 for Regression, and 3 for Survival)"
      << std::endl;
  std::cout << "    " << "--maxdepth N                  Set maximal tree depth to N." << std::endl;
  std::cout << "    "
      << "                              Set to 0 for unlimited depth. A value of 1 corresponds to tree stumps (1 split)."
      << std::endl;
  std::cout << "    "
      << "--catvars V1,V2,..            Comma separated list of names of (unordered) categorical variables. "
      << std::endl;
  std::cout << "    "
      << "                              Categorical variables must contain only positive integer values." << std::endl;
  std::cout << "    " << "--write                       Save forest to file <outprefix>.forest." << std::endl;
  std::cout << "    "
      << "--predict FILE                Load forest from FILE and predict with new data. The new data is expected in the exact same "
      << std::endl;
  std::cout << "    "
      << "                              shape as the training data. If the outcome of your new dataset is unknown, add a dummy column."
      << std::endl;
  std::cout << "    "
      << "--predall                     Return a matrix with individual predictions for each tree instead of aggregated "
      << std::endl;
  std::cout << "    " << "                              predictions for all trees (classification and regression only)."
      << std::endl;
  std::cout << "    " << "--predictiontype TYPE         Set type of prediction to:" << std::endl;
  std::cout << "    " << "                              TYPE = 1: Return predicted classes or values." << std::endl;
  std::cout << "    "
      << "                              TYPE = 2: Return terminal node IDs per tree for new observations." << std::endl;
  std::cout << "    " << "                              (Default: 1)" << std::endl;
  std::cout << "    " << "--impmeasure TYPE             Set importance mode to:" << std::endl;
  std::cout << "    " << "                              TYPE = 0: none." << std::endl;
  std::cout << "    "
      << "                              TYPE = 1: Node impurity: Gini for Classification, variance for Regression, sum of test statistic for Survival."
      << std::endl;
  std::cout << "    " << "                              TYPE = 2: Permutation importance, scaled by standard errors."
      << std::endl;
  std::cout << "    " << "                              TYPE = 3: Permutation importance, no scaling." << std::endl;
  std::cout << "    "
      << "                              TYPE = 5: Corrected node impurity: Bias-corrected version of node impurity importance."
      << std::endl;
  std::cout << "    " << "                              TYPE = 6: Local (casewise) permutation importance." << std::endl;
  std::cout << "    " << "                              (Default: 0)" << std::endl;
  std::cout << "    " << "--noreplace                   Sample without replacement." << std::endl;
  std::cout << "    "
      << "--fraction X                  Fraction of observations to sample. Default is 1 for sampling with replacement "
      << std::endl;
  std::cout << "    " << "                              and 0.632 for sampling without replacement." << std::endl;
  std::cout << "    " << "--splitrule RULE              Splitting rule:" << std::endl;
  std::cout << "    "
      << "                              RULE = 1: Gini for Classification, variance for Regression, logrank for Survival."
      << std::endl;
  std::cout << "    "
      << "                              RULE = 2: AUC for Survival, not available for Classification and Regression."
      << std::endl;
  std::cout << "    "
      << "                              RULE = 3: AUC (ignore ties) for Survival, not available for Classification and Regression."
      << std::endl;
  std::cout << "    "
      << "                              RULE = 4: MAXSTAT for Survival and Regression, not available for Classification."
      << std::endl;
  std::cout << "    " << "                              RULE = 5: ExtraTrees for all tree types." << std::endl;
  std::cout << "    " << "                              RULE = 6: BETA for regression, only for (0,1) bounded outcomes." << std::endl;
  std::cout << "    " << "                              RULE = 7: Hellinger for Classification, not available for Regression and Survival." << std::endl;
  std::cout << "    " << "                              (Default: 1)" << std::endl;
  std::cout << "    "
      << "--randomsplits N              Number of random splits to consider for each splitting variable (ExtraTrees splitrule only)."
      << std::endl;
  std::cout << "    "
      << "--alpha VAL                   Significance threshold to allow splitting (MAXSTAT splitrule only)."
      << std::endl;
  std::cout << "    "
      << "--minprop VAL                 Lower quantile of covariate distribtuion to be considered for splitting (MAXSTAT splitrule only)."
      << std::endl;
  std::cout << "    " << "--caseweights FILE            Filename of case weights file." << std::endl;
  std::cout << "    "
      << "--holdout                     Hold-out mode. Hold-out all samples with case weight 0 and use these for variable "
      << std::endl;
  std::cout << "    " << "                              importance and prediction error." << std::endl;
  std::cout << "    " << "--splitweights FILE           Filename of split select weights file." << std::endl;
  std::cout << "    "
      << "--alwayssplitvars V1,V2,..    Comma separated list of variable names to be always considered for splitting."
      << std::endl;
  std::cout << "    "
      << "--regcoef r1,r2,..            Comma separated list of regularization coefficients (one for all variables or one for each variable)."
      << std::endl;
  std::cout << "    " << "--usedepth                    Use node depth for regularization." << std::endl;
  std::cout << "    " << "--skipoob                     Skip computation of OOB error." << std::endl;
  std::cout << "    " << "--nthreads N                  Set number of parallel threads to N." << std::endl;
  std::cout << "    " << "                              (Default: Number of CPUs available)" << std::endl;
  std::cout << "    " << "--seed SEED                   Set random seed to SEED." << std::endl;
  std::cout << "    " << "                              (Default: No seed)" << std::endl;
  std::cout << "    " << "--outprefix PREFIX            Prefix for output files." << std::endl;
  std::cout << "    " << "--memmode MODE                Set memory mode to:" << std::endl;
  std::cout << "    " << "                              MODE = 0: double." << std::endl;
  std::cout << "    " << "                              MODE = 1: float." << std::endl;
  std::cout << "    " << "                              MODE = 2: char." << std::endl;
  std::cout << "    " << "                              (Default: 0)" << std::endl;
  std::cout << "    " << "--savemem                     Use memory saving (but slower) splitting mode." << std::endl;
  std::cout << std::endl;

  std::cout << "See README file for details and examples." << std::endl;
}

void ArgumentHandler::displayVersion() {
  std::cout << "Ranger version: " << RANGER_VERSION << std::endl;
  std::cout << std::endl;
  std::cout << "Please cite Ranger: " << std::endl;
  std::cout
      << "Wright, M. N. & Ziegler, A. (2017). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. Journal of Statistical Software 77:1-17."
      << std::endl;
  std::cout << std::endl;
  std::cout << "BibTeX:" << std::endl;
  std::cout << "@Article{," << std::endl;
  std::cout
      << "    title = {{ranger}: A Fast Implementation of Random Forests for High Dimensional Data in {C++} and {R},"
      << std::endl;
  std::cout << "    author = {Wright, Marvin N. and Ziegler, Andreas}," << std::endl;
  std::cout << "    journal = {Journal of Statistical Software}," << std::endl;
  std::cout << "    year = {2017}," << std::endl;
  std::cout << "    number = {1}," << std::endl;
  std::cout << "    pages = {1--17}," << std::endl;
  std::cout << "    doi = {10.18637/jss.v077.i01}," << std::endl;
  std::cout << "}" << std::endl;
}

} // namespace ranger
