/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#ifndef TREECLASSIFICATION_H_
#define TREECLASSIFICATION_H_

#include <vector>

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

  TreeClassification(const TreeClassification&) = delete;
  TreeClassification& operator=(const TreeClassification&) = delete;

  virtual ~TreeClassification() override = default;

  void allocateMemory() override;

  double estimate(size_t nodeID);
  void computePermutationImportanceInternal(std::vector<std::vector<size_t>>* permutations);
  void appendToFileInternal(std::ofstream& file) override;

  double getPrediction(size_t sampleID) const {
    size_t terminal_nodeID = prediction_terminal_nodeIDs[sampleID];
    return split_values[terminal_nodeID];
  }

  size_t getPredictionTerminalNodeID(size_t sampleID) const {
    return prediction_terminal_nodeIDs[sampleID];
  }

private:
  bool splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) override;
  void createEmptyNodeInternal() override;

  double computePredictionAccuracyInternal(std::vector<double>* prediction_error_casewise) override;

  // Called by splitNodeInternal(). Sets split_varIDs and split_values.
  bool findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  void findBestSplitValueSmallQ(size_t nodeID, size_t varID, size_t num_classes,
      const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
      double& best_decrease);
  void findBestSplitValueSmallQ(size_t nodeID, size_t varID, size_t num_classes,
      const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
      double& best_decrease, const std::vector<double>& possible_split_values, std::vector<size_t>& counter_per_class,
      std::vector<size_t>& counter);
  void findBestSplitValueLargeQ(size_t nodeID, size_t varID, size_t num_classes,
      const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
      double& best_decrease);
  void findBestSplitValueUnordered(size_t nodeID, size_t varID, size_t num_classes,
      const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
      double& best_decrease);

  bool findBestSplitExtraTrees(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  void findBestSplitValueExtraTrees(size_t nodeID, size_t varID, size_t num_classes,
      const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
      double& best_decrease);
  void findBestSplitValueExtraTrees(size_t nodeID, size_t varID, size_t num_classes,
      const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
      double& best_decrease, const std::vector<double>& possible_split_values, std::vector<size_t>& class_counts_right,
      std::vector<size_t>& n_right);
  void findBestSplitValueExtraTreesUnordered(size_t nodeID, size_t varID, size_t num_classes,
      const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
      double& best_decrease);

  void addGiniImportance(size_t nodeID, size_t varID, double decrease);

  void bootstrapClassWise() override;
  void bootstrapWithoutReplacementClassWise() override;

  void cleanUpInternal() override {
    counter.clear();
    counter.shrink_to_fit();
    counter_per_class.clear();
    counter_per_class.shrink_to_fit();
  }

  // Classes of the dependent variable and classIDs for responses
  const std::vector<double>* class_values;
  const std::vector<uint>* response_classIDs;
  const std::vector<std::vector<size_t>>* sampleIDs_per_class;

  // Splitting weights
  const std::vector<double>* class_weights;

  std::vector<size_t> counter;
  std::vector<size_t> counter_per_class;
};

} // namespace ranger

#endif /* TREECLASSIFICATION_H_ */
