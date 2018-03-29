/*-------------------------------------------------------------------------------
This file is part of ranger.

Copyright (c) [2014-2018] [Marvin N. Wright]

This software may be modified and distributed under the terms of the MIT license.

Please note that the C++ core of ranger is distributed under MIT license and the
R package "ranger" under GPL3 license.
#-------------------------------------------------------------------------------*/

#ifndef TREEPROBABILITY_H_
#define TREEPROBABILITY_H_

#include <map>

#include "globals.h"
#include "Tree.h"

class TreeProbability: public Tree {
public:
  TreeProbability(std::vector<double>* class_values, std::vector<uint>* response_classIDs,
      std::vector<std::vector<size_t>>* sampleIDs_per_class, std::vector<double>* class_weights);

  // Create from loaded forest
  TreeProbability(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
      std::vector<double>& split_values, std::vector<double>* class_values, std::vector<uint>* response_classIDs,
      std::vector<std::vector<double>>& terminal_class_counts);

  virtual ~TreeProbability();

  void allocateMemory();

  void addToTerminalNodes(size_t nodeID);
  void computePermutationImportanceInternal(std::vector<std::vector<size_t>>* permutations);
  void appendToFileInternal(std::ofstream& file);

  const std::vector<double>& getPrediction(size_t sampleID) const {
    size_t terminal_nodeID = prediction_terminal_nodeIDs[sampleID];
    return terminal_class_counts[terminal_nodeID];
  }

  size_t getPredictionTerminalNodeID(size_t sampleID) const {
    return prediction_terminal_nodeIDs[sampleID];
  }

  const std::vector<std::vector<double>>& getTerminalClassCounts() const {
    return terminal_class_counts;
  }

private:
  bool splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  void createEmptyNodeInternal();

  double computePredictionAccuracyInternal();

  // Called by splitNodeInternal(). Sets split_varIDs and split_values.
  bool findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  void findBestSplitValueSmallQ(size_t nodeID, size_t varID, size_t num_classes, size_t* class_counts,
      size_t num_samples_node, double& best_value, size_t& best_varID, double& best_decrease);
  void findBestSplitValueLargeQ(size_t nodeID, size_t varID, size_t num_classes, size_t* class_counts,
      size_t num_samples_node, double& best_value, size_t& best_varID, double& best_decrease);
  void findBestSplitValueUnordered(size_t nodeID, size_t varID, size_t num_classes, size_t* class_counts,
      size_t num_samples_node, double& best_value, size_t& best_varID, double& best_decrease);

  bool findBestSplitExtraTrees(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  void findBestSplitValueExtraTrees(size_t nodeID, size_t varID, size_t num_classes, size_t* class_counts,
      size_t num_samples_node, double& best_value, size_t& best_varID, double& best_decrease);
  void findBestSplitValueExtraTreesUnordered(size_t nodeID, size_t varID, size_t num_classes, size_t* class_counts,
      size_t num_samples_node, double& best_value, size_t& best_varID, double& best_decrease);

  void addImpurityImportance(size_t nodeID, size_t varID, double decrease);

  void bootstrapClassWise();
  void bootstrapWithoutReplacementClassWise();

  void cleanUpInternal() {
    if (counter != 0) {
      delete[] counter;
    }
    if (counter_per_class != 0) {
      delete[] counter_per_class;
    }
  }

  // Classes of the dependent variable and classIDs for responses
  std::vector<double>* class_values;
  std::vector<uint>* response_classIDs;
  std::vector<std::vector<size_t>>* sampleIDs_per_class;

  // Class counts in terminal nodes. Empty for non-terminal nodes.
  std::vector<std::vector<double>> terminal_class_counts;

  // Splitting weights
  std::vector<double>* class_weights;

  size_t* counter;
  size_t* counter_per_class;

  DISALLOW_COPY_AND_ASSIGN(TreeProbability);
};

#endif /* TREEPROBABILITY_H_ */
