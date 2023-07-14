/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#include <math.h>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <ctime>
#include <functional>
#include <thread>
#include <chrono>

#include "utility.h"
#include "Forest.h"
#include "DataChar.h"
#include "DataDouble.h"
#include "DataFloat.h"

namespace ranger {

Forest::Forest() :
    verbose_out(0), num_trees(DEFAULT_NUM_TREE), mtry(0), min_node_size(0), min_bucket(0), num_independent_variables(0), seed(0), num_samples(
        0), prediction_mode(false), memory_mode(MEM_DOUBLE), sample_with_replacement(true), memory_saving_splitting(
        false), splitrule(DEFAULT_SPLITRULE), predict_all(false), keep_inbag(false), sample_fraction( { 1 }), holdout(
        false), prediction_type(DEFAULT_PREDICTIONTYPE), num_random_splits(DEFAULT_NUM_RANDOM_SPLITS), max_depth(
        DEFAULT_MAXDEPTH), alpha(DEFAULT_ALPHA), minprop(DEFAULT_MINPROP), num_threads(DEFAULT_NUM_THREADS), data { }, overall_prediction_error(
    NAN), importance_mode(DEFAULT_IMPORTANCE_MODE), regularization_usedepth(false), progress(0) {
}

// #nocov start
void Forest::initCpp(std::string dependent_variable_name, MemoryMode memory_mode, std::string input_file, uint mtry,
    std::string output_prefix, uint num_trees, std::ostream* verbose_out, uint seed, uint num_threads,
    std::string load_forest_filename, ImportanceMode importance_mode, uint min_node_size, uint min_bucket,
    std::string split_select_weights_file, const std::vector<std::string>& always_split_variable_names,
    std::string status_variable_name, bool sample_with_replacement,
    const std::vector<std::string>& unordered_variable_names, bool memory_saving_splitting, SplitRule splitrule,
    std::string case_weights_file, bool predict_all, double sample_fraction, double alpha, double minprop, bool holdout,
    PredictionType prediction_type, uint num_random_splits, uint max_depth,
    const std::vector<double>& regularization_factor, bool regularization_usedepth) {

  this->memory_mode = memory_mode;
  this->verbose_out = verbose_out;

  if (!dependent_variable_name.empty()) {
    if (status_variable_name.empty()) {
      this->dependent_variable_names = {dependent_variable_name};
    } else {
      this->dependent_variable_names = {dependent_variable_name, status_variable_name};
    }
  }

  // Set prediction mode
  bool prediction_mode = false;
  if (!load_forest_filename.empty()) {
    prediction_mode = true;
  }

  // Sample fraction default and convert to vector
  if (sample_fraction == 0) {
    if (sample_with_replacement) {
      sample_fraction = DEFAULT_SAMPLE_FRACTION_REPLACE;
    } else {
      sample_fraction = DEFAULT_SAMPLE_FRACTION_NOREPLACE;
    }
  }
  std::vector<double> sample_fraction_vector = { sample_fraction };

  if (prediction_mode) {
    loadDependentVariableNamesFromFile(load_forest_filename);
  }

  // Call other init function
  init(loadDataFromFile(input_file), mtry, output_prefix, num_trees, seed, num_threads, importance_mode,
      min_node_size, min_bucket, prediction_mode, sample_with_replacement, unordered_variable_names, memory_saving_splitting,
      splitrule, predict_all, sample_fraction_vector, alpha, minprop, holdout, prediction_type, num_random_splits,
      false, max_depth, regularization_factor, regularization_usedepth);

  if (prediction_mode) {
    loadFromFile(load_forest_filename);
  }
  // Set variables to be always considered for splitting
  if (!always_split_variable_names.empty()) {
    setAlwaysSplitVariables(always_split_variable_names);
  }

  // TODO: Read 2d weights for tree-wise split select weights
  // Load split select weights from file
  if (!split_select_weights_file.empty()) {
    std::vector<std::vector<double>> split_select_weights;
    split_select_weights.resize(1);
    loadDoubleVectorFromFile(split_select_weights[0], split_select_weights_file);
    if (split_select_weights[0].size() != num_independent_variables) {
      throw std::runtime_error("Number of split select weights is not equal to number of independent variables.");
    }
    setSplitWeightVector(split_select_weights);
  }

  // Load case weights from file
  if (!case_weights_file.empty()) {
    loadDoubleVectorFromFile(case_weights, case_weights_file);
    if (case_weights.size() != num_samples) {
      throw std::runtime_error("Number of case weights is not equal to number of samples.");
    }
  }

  // Sample from non-zero weights in holdout mode
  if (holdout && !case_weights.empty()) {
    size_t nonzero_weights = 0;
    for (auto& weight : case_weights) {
      if (weight > 0) {
        ++nonzero_weights;
      }
    }
    this->sample_fraction[0] = this->sample_fraction[0] * ((double) nonzero_weights / (double) num_samples);
  }

  // Check if all catvars are coded in integers starting at 1
  if (!unordered_variable_names.empty()) {
    std::string error_message = checkUnorderedVariables(*data, unordered_variable_names);
    if (!error_message.empty()) {
      throw std::runtime_error(error_message);
    }
  }
}
// #nocov end

void Forest::initR(std::unique_ptr<Data> input_data, uint mtry, uint num_trees, std::ostream* verbose_out, uint seed,
    uint num_threads, ImportanceMode importance_mode, uint min_node_size, uint min_bucket,
    std::vector<std::vector<double>>& split_select_weights, const std::vector<std::string>& always_split_variable_names,
    bool prediction_mode, bool sample_with_replacement, const std::vector<std::string>& unordered_variable_names,
    bool memory_saving_splitting, SplitRule splitrule, std::vector<double>& case_weights,
    std::vector<std::vector<size_t>>& manual_inbag, bool predict_all, bool keep_inbag,
    std::vector<double>& sample_fraction, double alpha, double minprop, bool holdout, PredictionType prediction_type,
    uint num_random_splits, bool order_snps, uint max_depth, const std::vector<double>& regularization_factor,
    bool regularization_usedepth) {

  this->verbose_out = verbose_out;

  // Call other init function
  init(std::move(input_data), mtry, "", num_trees, seed, num_threads, importance_mode, min_node_size, min_bucket,
      prediction_mode, sample_with_replacement, unordered_variable_names, memory_saving_splitting, splitrule,
      predict_all, sample_fraction, alpha, minprop, holdout, prediction_type, num_random_splits, order_snps, max_depth,
      regularization_factor, regularization_usedepth);

  // Set variables to be always considered for splitting
  if (!always_split_variable_names.empty()) {
    setAlwaysSplitVariables(always_split_variable_names);
  }

  // Set split select weights
  if (!split_select_weights.empty()) {
    setSplitWeightVector(split_select_weights);
  }

  // Set case weights
  if (!case_weights.empty()) {
    if (case_weights.size() != num_samples) {
      throw std::runtime_error("Number of case weights not equal to number of samples.");
    }
    this->case_weights = case_weights;
  }

  // Set manual inbag
  if (!manual_inbag.empty()) {
    this->manual_inbag = manual_inbag;
  }

  // Keep inbag counts
  this->keep_inbag = keep_inbag;
}

void Forest::init(std::unique_ptr<Data> input_data, uint mtry, std::string output_prefix,
    uint num_trees, uint seed, uint num_threads, ImportanceMode importance_mode, uint min_node_size, uint min_bucket,
    bool prediction_mode, bool sample_with_replacement, const std::vector<std::string>& unordered_variable_names,
    bool memory_saving_splitting, SplitRule splitrule, bool predict_all, std::vector<double>& sample_fraction,
    double alpha, double minprop, bool holdout, PredictionType prediction_type, uint num_random_splits, bool order_snps,
    uint max_depth, const std::vector<double>& regularization_factor, bool regularization_usedepth) {

  // Initialize data with memmode
  this->data = std::move(input_data);

  // Initialize random number generator and set seed
  if (seed == 0) {
    std::random_device random_device;
    random_number_generator.seed(random_device());
  } else {
    random_number_generator.seed(seed);
  }

  // Set number of threads
  if (num_threads == DEFAULT_NUM_THREADS) {
    this->num_threads = std::thread::hardware_concurrency();
  } else {
    this->num_threads = num_threads;
  }

  // Set member variables
  this->num_trees = num_trees;
  this->mtry = mtry;
  this->seed = seed;
  this->output_prefix = output_prefix;
  this->importance_mode = importance_mode;
  this->min_node_size = min_node_size;
  this->min_bucket = min_bucket;
  this->prediction_mode = prediction_mode;
  this->sample_with_replacement = sample_with_replacement;
  this->memory_saving_splitting = memory_saving_splitting;
  this->splitrule = splitrule;
  this->predict_all = predict_all;
  this->sample_fraction = sample_fraction;
  this->holdout = holdout;
  this->alpha = alpha;
  this->minprop = minprop;
  this->prediction_type = prediction_type;
  this->num_random_splits = num_random_splits;
  this->max_depth = max_depth;
  this->regularization_factor = regularization_factor;
  this->regularization_usedepth = regularization_usedepth;

  // Set number of samples and variables
  num_samples = data->getNumRows();
  num_independent_variables = data->getNumCols();

  // Set unordered factor variables
  if (!prediction_mode) {
    data->setIsOrderedVariable(unordered_variable_names);
  }

  initInternal();

  // Init split select weights
  split_select_weights.push_back(std::vector<double>());

  // Init manual inbag
  manual_inbag.push_back(std::vector<size_t>());

  // Check if mtry is in valid range
  if (this->mtry > num_independent_variables) {
    throw std::runtime_error("mtry can not be larger than number of variables in data.");
  }

  // Check if any observations samples
  if ((size_t) num_samples * sample_fraction[0] < 1) {
    throw std::runtime_error("sample_fraction too small, no observations sampled.");
  }

  // Permute samples for corrected Gini importance
  if (importance_mode == IMP_GINI_CORRECTED) {
    data->permuteSampleIDs(random_number_generator);
  }

  // Order SNP levels if in "order" splitting
  if (!prediction_mode && order_snps) {
    data->orderSnpLevels((importance_mode == IMP_GINI_CORRECTED));
  }

  // Regularization
  if (regularization_factor.size() > 0) {
    if (regularization_factor.size() == 1 && num_independent_variables > 1) {
      double single_regularization_factor = regularization_factor[0];
      this->regularization_factor.resize(num_independent_variables, single_regularization_factor);
    } else if (regularization_factor.size() != num_independent_variables) {
      throw std::runtime_error("Use 1 or p (the number of predictor variables) regularization factors.");
    }

    // Set all variables to not used
    split_varIDs_used.resize(num_independent_variables, false);
  }
}

void Forest::run(bool verbose, bool compute_oob_error) {

  if (prediction_mode) {
    if (verbose && verbose_out) {
      *verbose_out << "Predicting .." << std::endl;
    }
    predict();
  } else {
    if (verbose && verbose_out) {
      *verbose_out << "Growing trees .." << std::endl;
    }

    grow();

    if (verbose && verbose_out) {
      *verbose_out << "Computing prediction error .." << std::endl;
    }

    if (compute_oob_error) {
      computePredictionError();
    }

    if (importance_mode == IMP_PERM_BREIMAN || importance_mode == IMP_PERM_LIAW || importance_mode == IMP_PERM_RAW
        || importance_mode == IMP_PERM_CASEWISE) {
      if (verbose && verbose_out) {
        *verbose_out << "Computing permutation variable importance .." << std::endl;
      }
      computePermutationImportance();
    }
  }
}

// #nocov start
void Forest::writeOutput() {

  if (verbose_out)
    *verbose_out << std::endl;
  writeOutputInternal();
  if (verbose_out) {
    if (dependent_variable_names.size() >= 1) {
      *verbose_out << "Dependent variable name:           " << dependent_variable_names[0] << std::endl;
    }
    *verbose_out << "Number of trees:                   " << num_trees << std::endl;
    *verbose_out << "Sample size:                       " << num_samples << std::endl;
    *verbose_out << "Number of independent variables:   " << num_independent_variables << std::endl;
    *verbose_out << "Mtry:                              " << mtry << std::endl;
    *verbose_out << "Target node size:                  " << min_node_size << std::endl;
    *verbose_out << "Variable importance mode:          " << importance_mode << std::endl;
    *verbose_out << "Memory mode:                       " << memory_mode << std::endl;
    *verbose_out << "Seed:                              " << seed << std::endl;
    *verbose_out << "Number of threads:                 " << num_threads << std::endl;
    *verbose_out << std::endl;
  }

  if (prediction_mode) {
    writePredictionFile();
  } else {
    if (verbose_out) {
      *verbose_out << "Overall OOB prediction error:      " << overall_prediction_error << std::endl;
      *verbose_out << std::endl;
    }

    if (!split_select_weights.empty() && !split_select_weights[0].empty()) {
      if (verbose_out) {
        *verbose_out
            << "Warning: Split select weights used. Variable importance measures are only comparable for variables with equal weights."
            << std::endl;
      }
    }

    if (importance_mode != IMP_NONE) {
      writeImportanceFile();
    }

    writeConfusionFile();
  }
}

void Forest::writeImportanceFile() {

  // Open importance file for writing
  std::string filename = output_prefix + ".importance";
  std::ofstream importance_file;
  importance_file.open(filename, std::ios::out);
  if (!importance_file.good()) {
    throw std::runtime_error("Could not write to importance file: " + filename + ".");
  }

  if (importance_mode == IMP_PERM_CASEWISE) {
    // Write variable names
    for (auto& variable_name : data->getVariableNames()) {
      importance_file << variable_name << " ";
    }
    importance_file << std::endl;

    // Write importance values
    for (size_t i = 0; i < num_samples; ++i) {
      for (size_t j = 0; j < num_independent_variables; ++j) {
        if (variable_importance_casewise.size() <= (j * num_samples + i)) {
          throw std::runtime_error("Memory error in local variable importance.");
        }
        importance_file << variable_importance_casewise[j * num_samples + i] << " ";
      }
      importance_file << std::endl;
    }
  } else {
    // Write importance to file
    for (size_t i = 0; i < variable_importance.size(); ++i) {
      std::string variable_name = data->getVariableNames()[i];
      importance_file << variable_name << ": " << variable_importance[i] << std::endl;
    }
  }

  importance_file.close();
  if (verbose_out)
    *verbose_out << "Saved variable importance to file " << filename << "." << std::endl;
}

void Forest::saveToFile() {

  // Open file for writing
  std::string filename = output_prefix + ".forest";
  std::ofstream outfile;
  outfile.open(filename, std::ios::binary);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to output file: " + filename + ".");
  }

  // Write dependent variable names
  uint num_dependent_variables = dependent_variable_names.size();
  if (num_dependent_variables >= 1) {
    outfile.write((char*) &num_dependent_variables, sizeof(num_dependent_variables));
    for (auto& var_name : dependent_variable_names) {
      size_t length = var_name.size();
      outfile.write((char*) &length, sizeof(length));
      outfile.write((char*) var_name.c_str(), length * sizeof(char));
    }
  } else {
    throw std::runtime_error("Missing dependent variable name.");
  }

  // Write num_trees
  outfile.write((char*) &num_trees, sizeof(num_trees));

  // Write is_ordered_variable
  saveVector1D(data->getIsOrderedVariable(), outfile);

  saveToFileInternal(outfile);

  // Write tree data for each tree
  for (auto& tree : trees) {
    tree->appendToFile(outfile);
  }

  // Close file
  outfile.close();
  if (verbose_out)
    *verbose_out << "Saved forest to file " << filename << "." << std::endl;
}
// #nocov end

void Forest::grow() {

  // Create thread ranges
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);

  // Call special grow functions of subclasses. There trees must be created.
  growInternal();

  // Init trees, create a seed for each tree, based on main seed
  std::uniform_int_distribution<uint> udist;
  for (size_t i = 0; i < num_trees; ++i) {
    uint tree_seed;
    if (seed == 0) {
      tree_seed = udist(random_number_generator);
    } else {
      tree_seed = (i + 1) * seed;
    }

    // Get split select weights for tree
    std::vector<double>* tree_split_select_weights;
    if (split_select_weights.size() > 1) {
      tree_split_select_weights = &split_select_weights[i];
    } else {
      tree_split_select_weights = &split_select_weights[0];
    }

    // Get inbag counts for tree
    std::vector<size_t>* tree_manual_inbag;
    if (manual_inbag.size() > 1) {
      tree_manual_inbag = &manual_inbag[i];
    } else {
      tree_manual_inbag = &manual_inbag[0];
    }

    trees[i]->init(data.get(), mtry, num_samples, tree_seed, &deterministic_varIDs, tree_split_select_weights,
        importance_mode, min_node_size, min_bucket, sample_with_replacement, memory_saving_splitting, splitrule, &case_weights,
        tree_manual_inbag, keep_inbag, &sample_fraction, alpha, minprop, holdout, num_random_splits, max_depth,
        &regularization_factor, regularization_usedepth, &split_varIDs_used);
  }

  // Init variable importance
  variable_importance.resize(num_independent_variables, 0);

  // Grow trees in multiple threads
  progress = 0;
#ifdef R_BUILD
  aborted = false;
  aborted_threads = 0;
#endif

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  // Initialize importance per thread
  std::vector<std::vector<double>> variable_importance_threads(num_threads);

  for (uint i = 0; i < num_threads; ++i) {
    if (importance_mode == IMP_GINI || importance_mode == IMP_GINI_CORRECTED) {
      variable_importance_threads[i].resize(num_independent_variables, 0);
    }
    threads.emplace_back(&Forest::growTreesInThread, this, i, &(variable_importance_threads[i]));
  }
  showProgress("Growing trees..", num_trees);
  for (auto &thread : threads) {
    thread.join();
  }

#ifdef R_BUILD
  if (aborted_threads > 0) {
    throw std::runtime_error("User interrupt.");
  }
#endif

  // Sum thread importances
  if (importance_mode == IMP_GINI || importance_mode == IMP_GINI_CORRECTED) {
    variable_importance.resize(num_independent_variables, 0);
    for (size_t i = 0; i < num_independent_variables; ++i) {
      for (uint j = 0; j < num_threads; ++j) {
        variable_importance[i] += variable_importance_threads[j][i];
      }
    }
    variable_importance_threads.clear();
  }

  // Divide importance by number of trees
  if (importance_mode == IMP_GINI || importance_mode == IMP_GINI_CORRECTED) {
    for (auto& v : variable_importance) {
      v /= num_trees;
    }
  }
}

void Forest::predict() {

  // Predict trees in multiple threads and join the threads with the main thread
  progress = 0;
#ifdef R_BUILD
  aborted = false;
  aborted_threads = 0;
#endif

  // Predict
  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  for (uint i = 0; i < num_threads; ++i) {
    threads.emplace_back(&Forest::predictTreesInThread, this, i, data.get(), false);
  }
  showProgress("Predicting..", num_trees);
  for (auto &thread : threads) {
    thread.join();
  }

  // Aggregate predictions
  allocatePredictMemory();
  threads.clear();
  threads.reserve(num_threads);
  progress = 0;
  for (uint i = 0; i < num_threads; ++i) {
    threads.emplace_back(&Forest::predictInternalInThread, this, i);
  }
  showProgress("Aggregating predictions..", num_samples);
  for (auto &thread : threads) {
    thread.join();
  }

#ifdef R_BUILD
  if (aborted_threads > 0) {
    throw std::runtime_error("User interrupt.");
  }
#endif
}

void Forest::computePredictionError() {

  // Predict trees in multiple threads
  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  progress = 0;
  for (uint i = 0; i < num_threads; ++i) {
    threads.emplace_back(&Forest::predictTreesInThread, this, i, data.get(), true);
  }
  showProgress("Computing prediction error..", num_trees);
  for (auto &thread : threads) {
    thread.join();
  }

#ifdef R_BUILD
  if (aborted_threads > 0) {
    throw std::runtime_error("User interrupt.");
  }
#endif

  // Call special function for subclasses
  computePredictionErrorInternal();
}

void Forest::computePermutationImportance() {

  // Compute tree permutation importance in multiple threads
  progress = 0;
#ifdef R_BUILD
  aborted = false;
  aborted_threads = 0;
#endif

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  // Initialize importance and variance
  std::vector<std::vector<double>> variable_importance_threads(num_threads);
  std::vector<std::vector<double>> variance_threads(num_threads);
  std::vector<std::vector<double>> variable_importance_casewise_threads(num_threads);

  // Compute importance
  for (uint i = 0; i < num_threads; ++i) {
    variable_importance_threads[i].resize(num_independent_variables, 0);
    if (importance_mode == IMP_PERM_BREIMAN || importance_mode == IMP_PERM_LIAW) {
      variance_threads[i].resize(num_independent_variables, 0);
    }
    if (importance_mode == IMP_PERM_CASEWISE) {
      variable_importance_casewise_threads[i].resize(num_independent_variables * num_samples, 0);
    }
    threads.emplace_back(&Forest::computeTreePermutationImportanceInThread, this, i,
        std::ref(variable_importance_threads[i]), std::ref(variance_threads[i]),
        std::ref(variable_importance_casewise_threads[i]));
  }
  showProgress("Computing permutation importance..", num_trees);
  for (auto &thread : threads) {
    thread.join();
  }

#ifdef R_BUILD
  if (aborted_threads > 0) {
    throw std::runtime_error("User interrupt.");
  }
#endif

  // Sum thread importances
  variable_importance.resize(num_independent_variables, 0);
  for (size_t i = 0; i < num_independent_variables; ++i) {
    for (uint j = 0; j < num_threads; ++j) {
      variable_importance[i] += variable_importance_threads[j][i];
    }
  }
  variable_importance_threads.clear();

  // Sum thread variances
  std::vector<double> variance(num_independent_variables, 0);
  if (importance_mode == IMP_PERM_BREIMAN || importance_mode == IMP_PERM_LIAW) {
    for (size_t i = 0; i < num_independent_variables; ++i) {
      for (uint j = 0; j < num_threads; ++j) {
        variance[i] += variance_threads[j][i];
      }
    }
    variance_threads.clear();
  }

  // Sum thread casewise importances
  if (importance_mode == IMP_PERM_CASEWISE) {
    variable_importance_casewise.resize(num_independent_variables * num_samples, 0);
    for (size_t i = 0; i < variable_importance_casewise.size(); ++i) {
      for (uint j = 0; j < num_threads; ++j) {
        variable_importance_casewise[i] += variable_importance_casewise_threads[j][i];
      }
    }
    variable_importance_casewise_threads.clear();
  }

  for (size_t i = 0; i < variable_importance.size(); ++i) {
    variable_importance[i] /= num_trees;

    // Normalize by variance for scaled permutation importance
    if (importance_mode == IMP_PERM_BREIMAN || importance_mode == IMP_PERM_LIAW) {
      if (variance[i] != 0) {
        variance[i] = variance[i] / num_trees - variable_importance[i] * variable_importance[i];
        variable_importance[i] /= sqrt(variance[i] / num_trees);
      }
    }
  }

  if (importance_mode == IMP_PERM_CASEWISE) {
    for (size_t i = 0; i < variable_importance_casewise.size(); ++i) {
      variable_importance_casewise[i] /= num_trees;
    }
  }
}

void Forest::growTreesInThread(uint thread_idx, std::vector<double>* variable_importance) {
  if (thread_ranges.size() > thread_idx + 1) {
    for (size_t i = thread_ranges[thread_idx]; i < thread_ranges[thread_idx + 1]; ++i) {
      trees[i]->grow(variable_importance);

      // Check for user interrupt
#ifdef R_BUILD
      if (aborted) {
        std::unique_lock<std::mutex> lock(mutex);
        ++aborted_threads;
        condition_variable.notify_one();
        return;
      }
#endif

      // Increase progress by 1 tree
      std::unique_lock<std::mutex> lock(mutex);
      ++progress;
      condition_variable.notify_one();
    }
  }
}

void Forest::predictTreesInThread(uint thread_idx, const Data* prediction_data, bool oob_prediction) {
  if (thread_ranges.size() > thread_idx + 1) {
    for (size_t i = thread_ranges[thread_idx]; i < thread_ranges[thread_idx + 1]; ++i) {
      trees[i]->predict(prediction_data, oob_prediction);

      // Check for user interrupt
#ifdef R_BUILD
      if (aborted) {
        std::unique_lock<std::mutex> lock(mutex);
        ++aborted_threads;
        condition_variable.notify_one();
        return;
      }
#endif

      // Increase progress by 1 tree
      std::unique_lock<std::mutex> lock(mutex);
      ++progress;
      condition_variable.notify_one();
    }
  }
}

void Forest::predictInternalInThread(uint thread_idx) {
  // Create thread ranges
  std::vector<uint> predict_ranges;
  equalSplit(predict_ranges, 0, num_samples - 1, num_threads);

  if (predict_ranges.size() > thread_idx + 1) {
    for (size_t i = predict_ranges[thread_idx]; i < predict_ranges[thread_idx + 1]; ++i) {
      predictInternal(i);

      // Check for user interrupt
#ifdef R_BUILD
      if (aborted) {
        std::unique_lock<std::mutex> lock(mutex);
        ++aborted_threads;
        condition_variable.notify_one();
        return;
      }
#endif

      // Increase progress by 1 tree
      std::unique_lock<std::mutex> lock(mutex);
      ++progress;
      condition_variable.notify_one();
    }
  }
}

void Forest::computeTreePermutationImportanceInThread(uint thread_idx, std::vector<double>& importance,
    std::vector<double>& variance, std::vector<double>& importance_casewise) {
  if (thread_ranges.size() > thread_idx + 1) {
    for (size_t i = thread_ranges[thread_idx]; i < thread_ranges[thread_idx + 1]; ++i) {
      trees[i]->computePermutationImportance(importance, variance, importance_casewise);

      // Check for user interrupt
#ifdef R_BUILD
      if (aborted) {
        std::unique_lock<std::mutex> lock(mutex);
        ++aborted_threads;
        condition_variable.notify_one();
        return;
      }
#endif

      // Increase progress by 1 tree
      std::unique_lock<std::mutex> lock(mutex);
      ++progress;
      condition_variable.notify_one();
    }
  }
}

// #nocov start
void Forest::loadFromFile(std::string filename) {
  if (verbose_out)
    *verbose_out << "Loading forest from file " << filename << "." << std::endl;

  // Open file for reading
  std::ifstream infile;
  infile.open(filename, std::ios::binary);
  if (!infile.good()) {
    throw std::runtime_error("Could not read from input file: " + filename + ".");
  }

  // Skip dependent variable names (already read)
  uint num_dependent_variables;
  infile.read((char*) &num_dependent_variables, sizeof(num_dependent_variables));
  for (size_t i = 0; i < num_dependent_variables; ++i) {
    size_t length;
    infile.read((char*) &length, sizeof(size_t));
    infile.ignore(length);
  }

  // Read num_trees
  infile.read((char*) &num_trees, sizeof(num_trees));

  // Read is_ordered_variable
  readVector1D(data->getIsOrderedVariable(), infile);

  // Read tree data. This is different for tree types -> virtual function
  loadFromFileInternal(infile);

  infile.close();

  // Create thread ranges
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);
}

void Forest::loadDependentVariableNamesFromFile(std::string filename) {

  // Open file for reading
  std::ifstream infile;
  infile.open(filename, std::ios::binary);
  if (!infile.good()) {
    throw std::runtime_error("Could not read from input file: " + filename + ".");
  }

  // Read dependent variable names
  dependent_variable_names.clear();
  uint num_dependent_variables = 0;
  infile.read((char*) &num_dependent_variables, sizeof(num_dependent_variables));
  for (size_t i = 0; i < num_dependent_variables; ++i) {
    size_t length;
    infile.read((char*) &length, sizeof(size_t));
    char* temp = new char[length + 1];
    infile.read((char*) temp, length * sizeof(char));
    temp[length] = '\0';
    dependent_variable_names.push_back(temp);
    delete[] temp;
  }

  infile.close();
}

std::unique_ptr<Data> Forest::loadDataFromFile(const std::string& data_path) {
  std::unique_ptr<Data> result { };
  switch (memory_mode) {
  case MEM_DOUBLE:
    result = std::make_unique<DataDouble>();
    break;
  case MEM_FLOAT:
    result = std::make_unique<DataFloat>();
    break;
  case MEM_CHAR:
    result = std::make_unique<DataChar>();
    break;
  }

  if (verbose_out)
    *verbose_out << "Loading input file: " << data_path << "." << std::endl;
  bool found_rounding_error = result->loadFromFile(data_path, dependent_variable_names);
  if (found_rounding_error && verbose_out) {
    *verbose_out << "Warning: Rounding or Integer overflow occurred. Use FLOAT or DOUBLE precision to avoid this."
        << std::endl;
  }
  return result;
}
// #nocov end

void Forest::setSplitWeightVector(std::vector<std::vector<double>>& split_select_weights) {

  // Size should be 1 x num_independent_variables or num_trees x num_independent_variables
  if (split_select_weights.size() != 1 && split_select_weights.size() != num_trees) {
    throw std::runtime_error("Size of split select weights not equal to 1 or number of trees.");
  }

  // Reserve space
  size_t num_weights = num_independent_variables;
  if (importance_mode == IMP_GINI_CORRECTED) {
    num_weights = 2 * num_independent_variables;
  }
  if (split_select_weights.size() == 1) {
    this->split_select_weights[0].resize(num_weights);
  } else {
    this->split_select_weights.clear();
    this->split_select_weights.resize(num_trees, std::vector<double>(num_weights));
  }
  
  // Deterministic varIDs 
  std::vector<bool> is_deterministic(num_weights, false);
  for (auto it = deterministic_varIDs.cbegin(); it != deterministic_varIDs.cend(); ++it) {
    is_deterministic[*it] = true;
  }

  // Split up in deterministic and weighted variables, ignore zero weights
  for (size_t i = 0; i < split_select_weights.size(); ++i) {
    size_t num_zero_weights = 0;

    // Size should be 1 x num_independent_variables or num_trees x num_independent_variables
    if (split_select_weights[i].size() != num_independent_variables) {
      throw std::runtime_error("Number of split select weights not equal to number of independent variables.");
    }

    for (size_t j = 0; j < split_select_weights[i].size(); ++j) {
      double weight = split_select_weights[i][j];

      if (weight == 0 || is_deterministic[j]) {
        ++num_zero_weights;
      } else if (weight < 0 || weight > 1) {
        throw std::runtime_error("One or more split select weights not in range [0,1].");
      } else {
        this->split_select_weights[i][j] = weight;
      }
    }

    // Copy weights for corrected impurity importance
    if (importance_mode == IMP_GINI_CORRECTED) {
      std::vector<double>* sw = &(this->split_select_weights[i]);
      std::copy_n(sw->begin(), num_independent_variables, sw->begin() + num_independent_variables);
    }

    if (num_weights - num_zero_weights < mtry) {
      throw std::runtime_error("Too many zeros in split select weights. Need at least mtry variables to split at.");
    }
  }
}

void Forest::setAlwaysSplitVariables(const std::vector<std::string>& always_split_variable_names) {

  deterministic_varIDs.reserve(num_independent_variables);

  for (auto& variable_name : always_split_variable_names) {
    size_t varID = data->getVariableID(variable_name);
    deterministic_varIDs.push_back(varID);
  }

  if (deterministic_varIDs.size() + this->mtry > num_independent_variables) {
    throw std::runtime_error(
        "Number of variables to be always considered for splitting plus mtry cannot be larger than number of independent variables.");
  }

  // Also add variables for corrected impurity importance
  if (importance_mode == IMP_GINI_CORRECTED) {
    size_t num_deterministic_varIDs = deterministic_varIDs.size();
    for (size_t k = 0; k < num_deterministic_varIDs; ++k) {
      deterministic_varIDs.push_back(k + num_independent_variables);
    }
  }
  
  // Sort in reverse order for removing with erase later
  std::sort(deterministic_varIDs.rbegin(), deterministic_varIDs.rend());
}

void Forest::showProgress(std::string operation, size_t max_progress) {
  using std::chrono::steady_clock;
  using std::chrono::duration_cast;
  using std::chrono::seconds;

  steady_clock::time_point start_time = steady_clock::now();
  steady_clock::time_point last_time = steady_clock::now();
  std::unique_lock<std::mutex> lock(mutex);

  // Wait for message from threads and show output if enough time elapsed
  while (progress < max_progress) {
    condition_variable.wait(lock);
    seconds elapsed_time = duration_cast<seconds>(steady_clock::now() - last_time);

    // Check for user interrupt
#ifdef R_BUILD
    if (!aborted && checkInterrupt()) {
      aborted = true;
    }
    if (aborted && aborted_threads >= num_threads) {
      return;
    }
#endif

    if (progress > 0 && elapsed_time.count() > STATUS_INTERVAL) {
      double relative_progress = (double) progress / (double) max_progress;
      seconds time_from_start = duration_cast<seconds>(steady_clock::now() - start_time);
      uint remaining_time = (1 / relative_progress - 1) * time_from_start.count();
      if (verbose_out) {
        *verbose_out << operation << " Progress: " << round(100 * relative_progress) << "%. Estimated remaining time: "
            << beautifyTime(remaining_time) << "." << std::endl;
      }
      last_time = steady_clock::now();
    }
  }
}

} // namespace ranger
