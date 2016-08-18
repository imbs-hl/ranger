/*-------------------------------------------------------------------------------
 This file is part of Ranger.

 Ranger is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Ranger is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with Ranger. If not, see <http://www.gnu.org/licenses/>.

 Written by:

 Marvin N. Wright
 Institut für Medizinische Biometrie und Statistik
 Universität zu Lübeck
 Ratzeburger Allee 160
 23562 Lübeck

 http://www.imbs-luebeck.de
 wright@imbs.uni-luebeck.de
 #-------------------------------------------------------------------------------*/

#ifndef TREE_H_
#define TREE_H_

#include <vector>
#include <random>
#include <iostream>

#include "globals.h"
#include "Data.h"

class Tree {
public:
  Tree();

  // Create from loaded forest
  Tree(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
      std::vector<double>& split_values, std::vector<bool>* is_ordered_variable);

  virtual ~Tree();

  void init(Data* data, uint mtry, size_t dependent_varID, size_t num_samples, uint seed,
      std::vector<size_t>* deterministic_varIDs, std::vector<size_t>* split_select_varIDs,
      std::vector<double>* split_select_weights, ImportanceMode importance_mode, uint min_node_size,
      std::vector<size_t>* no_split_variables, bool sample_with_replacement, std::vector<bool>* is_unordered,
      bool memory_saving_splitting, SplitRule splitrule, std::vector<double>* case_weights, bool keep_inbag,
      double sample_fraction, double alpha, double minprop, bool holdout);

  virtual void initInternal() = 0;

  void grow(std::vector<double>* variable_importance);

  void predict(const Data* prediction_data, bool oob_prediction);

  void computePermutationImportance(std::vector<double>* forest_importance, std::vector<double>* forest_variance);

  void appendToFile(std::ofstream& file);
  virtual void appendToFileInternal(std::ofstream& file) = 0;

  const std::vector<std::vector<size_t> >& getChildNodeIDs() const {
    return child_nodeIDs;
  }
  const std::vector<double>& getSplitValues() const {
    return split_values;
  }
  const std::vector<size_t>& getSplitVarIDs() const {
    return split_varIDs;
  }

  const std::vector<size_t>& getOobSampleIDs() const {
    return oob_sampleIDs;
  }
  size_t getNumSamplesOob() const {
    return num_samples_oob;
  }

  const std::vector<size_t>& getInbagCounts() const {
    return inbag_counts;
  }

protected:
  void createPossibleSplitVarSubset(std::vector<size_t>& result);

  bool splitNode(size_t nodeID);
  virtual bool splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) = 0;

  void createEmptyNode();
  virtual void createEmptyNodeInternal() = 0;

  size_t dropDownSamplePermuted(size_t permuted_varID, size_t sampleID, size_t permuted_sampleID);
  void permuteAndPredictOobSamples(size_t permuted_varID, std::vector<size_t>& permutations);

  virtual double computePredictionAccuracyInternal() = 0;

  void bootstrap();
  void bootstrapWithoutReplacement();

  void bootstrapWeighted();
  void bootstrapWithoutReplacementWeighted();

  virtual void cleanUpInternal() = 0;

  size_t dependent_varID;
  uint mtry;

  // Number of samples (all samples, not only inbag for this tree)
  size_t num_samples;

  // Number of OOB samples
  size_t num_samples_oob;

  // For each varID true if unordered
  std::vector<bool>* is_ordered_variable;

  // Variable to not split at (only dependent_varID for non-survival trees)
  std::vector<size_t>* no_split_variables;

  // Minimum node size to split, like in original RF nodes of smaller size can be produced
  uint min_node_size;

  // Weight vector for selecting possible split variables, one weight between 0 (never select) and 1 (always select) for each variable
  // Deterministic variables are always selected
  std::vector<size_t>* deterministic_varIDs;
  std::vector<size_t>* split_select_varIDs;
  std::vector<double>* split_select_weights;

  // Bootstrap weights
  std::vector<double>* case_weights;

  // Splitting variable for each node
  std::vector<size_t> split_varIDs;

  // Value to split at for each node, for now only binary split
  // For terminal nodes the prediction value is saved here
  std::vector<double> split_values;

  // Vector of left and right child node IDs, 0 for no child
  std::vector<std::vector<size_t>> child_nodeIDs;

  // For each node a vector with IDs of samples in node
  std::vector<std::vector<size_t>> sampleIDs;

  // IDs of OOB individuals, sorted
  std::vector<size_t> oob_sampleIDs;

  // Holdout mode
  bool holdout;

  // Inbag counts
  bool keep_inbag;
  std::vector<size_t> inbag_counts;

  // Random number generator
  std::mt19937_64 random_number_generator;

  // Pointer to original data
  Data* data;

  // Variable importance for all variables
  std::vector<double>* variable_importance;
  ImportanceMode importance_mode;

  // When growing here the OOB set is used
  // Terminal nodeIDs for prediction samples
  std::vector<size_t> prediction_terminal_nodeIDs;

  bool sample_with_replacement;
  double sample_fraction;

  bool memory_saving_splitting;
  SplitRule splitrule;
  double alpha;
  double minprop;

private:
  DISALLOW_COPY_AND_ASSIGN(Tree);
};

#endif /* TREE_H_ */
