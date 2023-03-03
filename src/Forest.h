/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#ifndef FOREST_H_
#define FOREST_H_

#include <vector>
#include <iostream>
#include <random>
#include <ctime>
#include <memory>
#include <thread>
#include <chrono>
#include <mutex>
#include <condition_variable>

#include "globals.h"
#include "Tree.h"
#include "Data.h"

namespace ranger {

class Forest {
public:
  Forest();

  Forest(const Forest&) = delete;
  Forest& operator=(const Forest&) = delete;

  virtual ~Forest() = default;

  // Init from c++ main or Rcpp from R
  void initCpp(std::string dependent_variable_name, MemoryMode memory_mode, std::string input_file, uint mtry,
      std::string output_prefix, uint num_trees, std::ostream* verbose_out, uint seed, uint num_threads,
      std::string load_forest_filename, ImportanceMode importance_mode, uint min_node_size, uint min_bucket,
      std::string split_select_weights_file, const std::vector<std::string>& always_split_variable_names,
      std::string status_variable_name, bool sample_with_replacement,
      const std::vector<std::string>& unordered_variable_names, bool memory_saving_splitting, SplitRule splitrule,
      std::string case_weights_file, bool predict_all, double sample_fraction, double alpha, double minprop,
      bool holdout, PredictionType prediction_type, uint num_random_splits, uint max_depth,
      const std::vector<double>& regularization_factor, bool regularization_usedepth);
  void initR(std::unique_ptr<Data> input_data, uint mtry, uint num_trees, std::ostream* verbose_out, uint seed,
      uint num_threads, ImportanceMode importance_mode, uint min_node_size, uint min_bucket,
      std::vector<std::vector<double>>& split_select_weights,
      const std::vector<std::string>& always_split_variable_names, bool prediction_mode, bool sample_with_replacement,
      const std::vector<std::string>& unordered_variable_names, bool memory_saving_splitting, SplitRule splitrule,
      std::vector<double>& case_weights, std::vector<std::vector<size_t>>& manual_inbag, bool predict_all,
      bool keep_inbag, std::vector<double>& sample_fraction, double alpha, double minprop, bool holdout,
      PredictionType prediction_type, uint num_random_splits, bool order_snps, uint max_depth,
      const std::vector<double>& regularization_factor, bool regularization_usedepth);
  void init(std::unique_ptr<Data> input_data, uint mtry, std::string output_prefix,
      uint num_trees, uint seed, uint num_threads, ImportanceMode importance_mode, uint min_node_size, uint min_bucket,
      bool prediction_mode, bool sample_with_replacement, const std::vector<std::string>& unordered_variable_names,
      bool memory_saving_splitting, SplitRule splitrule, bool predict_all, std::vector<double>& sample_fraction,
      double alpha, double minprop, bool holdout, PredictionType prediction_type, uint num_random_splits,
      bool order_snps, uint max_depth, const std::vector<double>& regularization_factor, bool regularization_usedepth);
  virtual void initInternal() = 0;

  // Grow or predict
  void run(bool verbose, bool compute_oob_error);

  // Write results to output files
  void writeOutput();
  virtual void writeOutputInternal() = 0;
  virtual void writeConfusionFile() = 0;
  virtual void writePredictionFile() = 0;
  void writeImportanceFile();

  // Save forest to file
  void saveToFile();
  virtual void saveToFileInternal(std::ofstream& outfile) = 0;

  std::vector<std::vector<std::vector<size_t>>> getChildNodeIDs() {
    std::vector<std::vector<std::vector<size_t>>> result;
    for (auto& tree : trees) {
      result.push_back(tree->getChildNodeIDs());
    }
    return result;
  }
  std::vector<std::vector<size_t>> getSplitVarIDs() {
    std::vector<std::vector<size_t>> result;
    for (auto& tree : trees) {
      result.push_back(tree->getSplitVarIDs());
    }
    return result;
  }
  std::vector<std::vector<double>> getSplitValues() {
    std::vector<std::vector<double>> result;
    for (auto& tree : trees) {
      result.push_back(tree->getSplitValues());
    }
    return result;
  }
  const std::vector<double>& getVariableImportance() const {
    return variable_importance;
  }
  const std::vector<double>& getVariableImportanceCasewise() const {
    return variable_importance_casewise;
  }
  double getOverallPredictionError() const {
    return overall_prediction_error;
  }
  const std::vector<std::vector<std::vector<double>>>& getPredictions() const {
    return predictions;
  }
  size_t getNumTrees() const {
    return num_trees;
  }
  uint getMtry() const {
    return mtry;
  }
  uint getMinNodeSize() const {
    return min_node_size;
  }
  uint getMinBucket() const {
    return min_bucket;
  }
  size_t getNumIndependentVariables() const {
    return num_independent_variables;
  }

  const std::vector<bool>& getIsOrderedVariable() const {
    return data->getIsOrderedVariable();
  }

  std::vector<std::vector<size_t>> getInbagCounts() const {
    std::vector<std::vector<size_t>> result;
    for (auto& tree : trees) {
      result.push_back(tree->getInbagCounts());
    }
    return result;
  }

  const std::vector<std::vector<size_t>>& getSnpOrder() const {
    return data->getSnpOrder();
  }

protected:
  void grow();
  virtual void growInternal() = 0;

  // Predict using existing tree from file and data as prediction data
  void predict();
  virtual void allocatePredictMemory() = 0;
  virtual void predictInternal(size_t sample_idx) = 0;

  void computePredictionError();
  virtual void computePredictionErrorInternal() = 0;

  void computePermutationImportance();

  // Multithreading methods for growing/prediction/importance, called by each thread
  void growTreesInThread(uint thread_idx, std::vector<double>* variable_importance);
  void predictTreesInThread(uint thread_idx, const Data* prediction_data, bool oob_prediction);
  void predictInternalInThread(uint thread_idx);
  void computeTreePermutationImportanceInThread(uint thread_idx, std::vector<double>& importance,
      std::vector<double>& variance, std::vector<double>& importance_casewise);

  // Load forest from file
  void loadFromFile(std::string filename);
  virtual void loadFromFileInternal(std::ifstream& infile) = 0;
  void loadDependentVariableNamesFromFile(std::string filename);

  // Load data from file
  std::unique_ptr<Data> loadDataFromFile(const std::string& data_path);

  // Set split select weights and variables to be always considered for splitting
  void setSplitWeightVector(std::vector<std::vector<double>>& split_select_weights);
  void setAlwaysSplitVariables(const std::vector<std::string>& always_split_variable_names);

  // Show progress every few seconds
  void showProgress(std::string operation, size_t max_progress);

  // Verbose output stream, cout if verbose==true, logfile if not
  std::ostream* verbose_out;

  std::vector<std::string> dependent_variable_names; // time,status for survival
  size_t num_trees;
  uint mtry;
  uint min_node_size;
  uint min_bucket;
  size_t num_independent_variables;
  uint seed;
  size_t num_samples;
  bool prediction_mode;
  MemoryMode memory_mode;
  bool sample_with_replacement;
  bool memory_saving_splitting;
  SplitRule splitrule;
  bool predict_all;
  bool keep_inbag;
  std::vector<double> sample_fraction;
  bool holdout;
  PredictionType prediction_type;
  uint num_random_splits;
  uint max_depth;

  // MAXSTAT splitrule
  double alpha;
  double minprop;

  // Multithreading
  uint num_threads;
  std::vector<uint> thread_ranges;
  std::mutex mutex;
  std::condition_variable condition_variable;

  std::vector<std::unique_ptr<Tree>> trees;
  std::unique_ptr<Data> data;

  std::vector<std::vector<std::vector<double>>> predictions;
  double overall_prediction_error;

  // Weight vector for selecting possible split variables, one weight between 0 (never select) and 1 (always select) for each variable
  // Deterministic variables are always selected
  std::vector<size_t> deterministic_varIDs;
  std::vector<std::vector<double>> split_select_weights;

  // Bootstrap weights
  std::vector<double> case_weights;

  // Pre-selected bootstrap samples (per tree)
  std::vector<std::vector<size_t>> manual_inbag;

  // Random number generator
  std::mt19937_64 random_number_generator;

  std::string output_prefix;
  ImportanceMode importance_mode;

  // Regularization
  std::vector<double> regularization_factor;
  bool regularization_usedepth;
  std::vector<bool> split_varIDs_used;
  
  // Variable importance for all variables in forest
  std::vector<double> variable_importance;

  // Casewise variable importance for all variables in forest
  std::vector<double> variable_importance_casewise;

  // Computation progress (finished trees)
  size_t progress;
#ifdef R_BUILD
  size_t aborted_threads;
  bool aborted;
#endif
};

} // namespace ranger

#endif /* FOREST_H_ */
