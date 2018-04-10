/*-------------------------------------------------------------------------------
This file is part of ranger.

Copyright (c) [2014-2018] [Marvin N. Wright]

This software may be modified and distributed under the terms of the MIT license.

Please note that the C++ core of ranger is distributed under MIT license and the
R package "ranger" under GPL3 license.
#-------------------------------------------------------------------------------*/

#ifndef TREEREGRESSION_H_
#define TREEREGRESSION_H_

#include "globals.h"
#include "Tree.h"

namespace ranger {

class TreeRegression: public Tree {
public:
  TreeRegression();

  // Create from loaded forest
  TreeRegression(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
      std::vector<double>& split_values);

  virtual ~TreeRegression();

  void allocateMemory();

  double estimate(size_t nodeID);
  void computePermutationImportanceInternal(std::vector<std::vector<size_t>>* permutations);
  void appendToFileInternal(std::ofstream& file);

  double getPrediction(size_t sampleID) const {
    size_t terminal_nodeID = prediction_terminal_nodeIDs[sampleID];
    return (split_values[terminal_nodeID]);
  }

  size_t getPredictionTerminalNodeID(size_t sampleID) const {
    return prediction_terminal_nodeIDs[sampleID];
  }

private:
  bool splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  void createEmptyNodeInternal();

  double computePredictionAccuracyInternal();

  // Called by splitNodeInternal(). Sets split_varIDs and split_values.
  bool findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  void findBestSplitValueSmallQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
      double& best_value, size_t& best_varID, double& best_decrease);
  void findBestSplitValueLargeQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
      double& best_value, size_t& best_varID, double& best_decrease);
  void findBestSplitValueUnordered(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
      double& best_value, size_t& best_varID, double& best_decrease);

  bool findBestSplitMaxstat(size_t nodeID, std::vector<size_t>& possible_split_varIDs);

  bool findBestSplitExtraTrees(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  void findBestSplitValueExtraTrees(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
      double& best_value, size_t& best_varID, double& best_decrease);
  void findBestSplitValueExtraTreesUnordered(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
      double& best_value, size_t& best_varID, double& best_decrease);

  void addImpurityImportance(size_t nodeID, size_t varID, double decrease);

  double computePredictionMSE();

  void cleanUpInternal() {
    if (counter != 0) {
      delete[] counter;
    }
    if (sums != 0) {
      delete[] sums;
    }
  }

  size_t* counter;
  double* sums;

  DISALLOW_COPY_AND_ASSIGN(TreeRegression);
};

} // namespace ranger

#endif /* TREEREGRESSION_H_ */
