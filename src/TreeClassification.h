/*-------------------------------------------------------------------------------
This file is part of ranger.

Copyright (c) [2014-2018] [Marvin N. Wright]

This software may be modified and distributed under the terms of the MIT license.

Please note that the C++ core of ranger is distributed under MIT license and the
R package "ranger" under GPL3 license.
#-------------------------------------------------------------------------------*/

#ifndef TREECLASSIFICATION_H_
#define TREECLASSIFICATION_H_

#include "globals.h"
#include "Tree.h"

namespace ranger {

class TreeClassification: public Tree {
public:
  TreeClassification(std::vector<double>* class_values, std::vector<uint>* response_classIDs,
      std::vector<std::vector<size_t>>* sampleIDs_per_class, std::vector<double>* class_weights);

  // Create from loaded forest
  TreeClassification(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
      std::vector<double>& split_values, std::vector<double>* class_values, std::vector<uint>* response_classIDs);

  virtual ~TreeClassification();

  TreeClassification(const TreeClassification&)            = delete;
  TreeClassification& operator=(const TreeClassification&) = delete;

  void allocateMemory();

  double estimate(size_t nodeID);
  void computePermutationImportanceInternal(std::vector<std::vector<size_t>>* permutations);
  void appendToFileInternal(std::ofstream& file);

  double getPrediction(size_t sampleID) const {
    size_t terminal_nodeID = prediction_terminal_nodeIDs[sampleID];
    return split_values[terminal_nodeID];
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

  void addGiniImportance(size_t nodeID, size_t varID, double decrease);

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

  // Splitting weights
  std::vector<double>* class_weights;

  size_t* counter;
  size_t* counter_per_class;
};

} // namespace ranger

#endif /* TREECLASSIFICATION_H_ */
