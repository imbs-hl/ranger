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
#include <vector>

#include "globals.h"
#include "Tree.h"

namespace ranger {

class TreeProbability: public Tree {
public:
  TreeProbability(std::vector<double>* class_values, std::vector<uint>* response_classIDs,
      std::vector<std::vector<size_t>>* sampleIDs_per_class, std::vector<double>* class_weights);

  // Create from loaded forest
  TreeProbability(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
      std::vector<double>& split_values, std::vector<double>* class_values, std::vector<uint>* response_classIDs,
      std::vector<std::vector<double>>& terminal_class_counts);

  TreeProbability(const TreeProbability&) = delete;
  TreeProbability& operator=(const TreeProbability&) = delete;

  virtual ~TreeProbability() override = default;

  void allocateMemory() override;

  void addToTerminalNodes(size_t nodeID);
  void computePermutationImportanceInternal(std::vector<std::vector<size_t>>* permutations);
  void appendToFileInternal(std::ofstream& file) override;

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
  bool splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) override;
  void createEmptyNodeInternal() override;

  double computePredictionAccuracyInternal() override;

  // Called by splitNodeInternal(). Sets split_varIDs and split_values.
  bool findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  void findBestSplitValueSmallQ(size_t nodeID, size_t varID, size_t num_classes,
      const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
      double& best_decrease);
  void findBestSplitValueSmallQ(size_t nodeID, size_t varID, size_t num_classes,
      const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
      double& best_decrease, const std::vector<double>& possible_split_values, std::vector<size_t>& class_counts_right,
      std::vector<size_t>& n_right);
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

  void addImpurityImportance(size_t nodeID, size_t varID, double decrease);

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

  // Class counts in terminal nodes. Empty for non-terminal nodes.
  std::vector<std::vector<double>> terminal_class_counts;

  // Splitting weights
  const std::vector<double>* class_weights;

  std::vector<size_t> counter;
  std::vector<size_t> counter_per_class;
};

} // namespace ranger

#endif /* TREEPROBABILITY_H_ */
