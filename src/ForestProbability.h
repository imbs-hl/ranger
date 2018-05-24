/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#ifndef FORESTPROBABILITY_H_
#define FORESTPROBABILITY_H_

#include <map>
#include <utility>
#include <vector>

#include "globals.h"
#include "Forest.h"
#include "TreeProbability.h"

namespace ranger {

class ForestProbability: public Forest {
public:
  ForestProbability() = default;

  ForestProbability(const ForestProbability&) = delete;
  ForestProbability& operator=(const ForestProbability&) = delete;

  virtual ~ForestProbability() override = default;

  void loadForest(size_t dependent_varID, size_t num_trees,
      std::vector<std::vector<std::vector<size_t>> >& forest_child_nodeIDs,
      std::vector<std::vector<size_t>>& forest_split_varIDs, std::vector<std::vector<double>>& forest_split_values,
      std::vector<double>& class_values, std::vector<std::vector<std::vector<double>>>& forest_terminal_class_counts,
      std::vector<bool>& is_ordered_variable);

  std::vector<std::vector<std::vector<double>>> getTerminalClassCounts() const;

  const std::vector<double>& getClassValues() const {
    return class_values;
  }

  void setClassWeights(std::vector<double>& class_weights) {
    this->class_weights = class_weights;
  }

protected:
  void initInternal(std::string status_variable_name) override;
  void growInternal() override;
  void allocatePredictMemory() override;
  void predictInternal(size_t sample_idx) override;
  void computePredictionErrorInternal() override;
  void writeOutputInternal() override;
  void writeConfusionFile() override;
  void writePredictionFile() override;
  void saveToFileInternal(std::ofstream& outfile) override;
  void loadFromFileInternal(std::ifstream& infile) override;

  // Classes of the dependent variable and classIDs for responses
  std::vector<double> class_values;
  std::vector<uint> response_classIDs;
  std::vector<std::vector<size_t>> sampleIDs_per_class;

  // Splitting weights
  std::vector<double> class_weights;

private:
  const std::vector<double>& getTreePrediction(size_t tree_idx, size_t sample_idx) const;
  size_t getTreePredictionTerminalNodeID(size_t tree_idx, size_t sample_idx) const;
};

} // namespace ranger

#endif /* FORESTPROBABILITY_H_ */
