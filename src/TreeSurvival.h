/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#ifndef TREESURVIVAL_H_
#define TREESURVIVAL_H_

#include <vector>

#include "globals.h"
#include "Tree.h"

namespace ranger {

class TreeSurvival: public Tree {
public:
  TreeSurvival(std::vector<double>* unique_timepoints, std::vector<size_t>* response_timepointIDs);

  // Create from loaded forest
  TreeSurvival(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
      std::vector<double>& split_values, std::vector<std::vector<double>> chf, std::vector<double>* unique_timepoints,
      std::vector<size_t>* response_timepointIDs);

  TreeSurvival(const TreeSurvival&) = delete;
  TreeSurvival& operator=(const TreeSurvival&) = delete;

  virtual ~TreeSurvival() override = default;

  void allocateMemory() override;

  void appendToFileInternal(std::ofstream& file) override;
  void computePermutationImportanceInternal(std::vector<std::vector<size_t>>* permutations);

  const std::vector<std::vector<double>>& getChf() const {
    return chf;
  }

  const std::vector<double>& getPrediction(size_t sampleID) const {
    size_t terminal_nodeID = prediction_terminal_nodeIDs[sampleID];
    return chf[terminal_nodeID];
  }

  size_t getPredictionTerminalNodeID(size_t sampleID) const {
    return prediction_terminal_nodeIDs[sampleID];
  }

private:

  void createEmptyNodeInternal() override;
  void computeSurvival(size_t nodeID);
  double computePredictionAccuracyInternal(std::vector<double>* prediction_error_casewise) override;
  
  bool splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) override;

  bool findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  bool findBestSplitMaxstat(size_t nodeID, std::vector<size_t>& possible_split_varIDs);

  void findBestSplitValueLogRank(size_t nodeID, size_t varID, std::vector<double>& possible_split_values,
      double& best_value, size_t& best_varID, double& best_logrank);
  void findBestSplitValueLogRankUnordered(size_t nodeID, size_t varID, std::vector<double>& factor_levels,
      double& best_value, size_t& best_varID, double& best_logrank);
  void findBestSplitValueAUC(size_t nodeID, size_t varID, double& best_value, size_t& best_varID, double& best_auc);

  void computeDeathCounts(size_t nodeID);
  void computeChildDeathCounts(size_t nodeID, size_t varID, std::vector<double>& possible_split_values,
      std::vector<size_t>& num_samples_right_child, std::vector<size_t>& num_samples_at_risk_right_child,
      std::vector<size_t>& num_deaths_right_child, size_t num_splits);

  void computeAucSplit(double time_k, double time_l, double status_k, double status_l, double value_k, double value_l,
      size_t num_splits, std::vector<double>& possible_split_values, std::vector<double>& num_count,
      std::vector<double>& num_total);

  void findBestSplitValueLogRank(size_t nodeID, size_t varID, double& best_value, size_t& best_varID,
      double& best_logrank);
  void findBestSplitValueLogRankUnordered(size_t nodeID, size_t varID, double& best_value, size_t& best_varID,
      double& best_logrank);

  bool findBestSplitExtraTrees(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  void findBestSplitValueExtraTrees(size_t nodeID, size_t varID, double& best_value, size_t& best_varID,
      double& best_logrank);
  void findBestSplitValueExtraTreesUnordered(size_t nodeID, size_t varID, double& best_value, size_t& best_varID,
      double& best_logrank);

  void addImpurityImportance(size_t nodeID, size_t varID, double decrease);

  void cleanUpInternal() override {
    num_deaths.clear();
    num_deaths.shrink_to_fit();
    num_samples_at_risk.clear();
    num_samples_at_risk.shrink_to_fit();
  }

  // Unique time points for all individuals (not only this bootstrap), sorted
  const std::vector<double>* unique_timepoints;
  size_t num_timepoints;
  const std::vector<size_t>* response_timepointIDs;

  // For all terminal nodes CHF for all unique timepoints. For other nodes empty vector.
  std::vector<std::vector<double>> chf;

  // Fields to save to while tree growing
  std::vector<size_t> num_deaths;
  std::vector<size_t> num_samples_at_risk;
};

} // namespace ranger

#endif /* TREESURVIVAL_H_ */
