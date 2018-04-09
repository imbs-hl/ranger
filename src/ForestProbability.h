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

class ForestProbability: public Forest {
public:
  ForestProbability();
  virtual ~ForestProbability();

  ForestProbability(const ForestProbability&)            = delete;
  ForestProbability& operator=(const ForestProbability&) = delete;

  void loadForest(size_t dependent_varID, size_t num_trees,
      std::vector<std::vector<std::vector<size_t>> >& forest_child_nodeIDs,
      std::vector<std::vector<size_t>>& forest_split_varIDs, std::vector<std::vector<double>>& forest_split_values,
      std::vector<double>& class_values, std::vector<std::vector<std::vector<double>>>& forest_terminal_class_counts,
      std::vector<bool>& is_ordered_variable);

  std::vector<std::vector<std::vector<double>>> getTerminalClassCounts() {
    std::vector<std::vector<std::vector<double>>> result;
    result.reserve(num_trees);
    for (Tree* tree : trees) {
      TreeProbability* temp = (TreeProbability*) tree;
      result.push_back(temp->getTerminalClassCounts());
    }
    return result;
  }

  const std::vector<double>& getClassValues() const {
    return class_values;
  }

  void setClassWeights(std::vector<double>& class_weights) {
    this->class_weights = class_weights;
  }

protected:
  void initInternal(std::string status_variable_name);
  void growInternal();
  void allocatePredictMemory();
  void predictInternal(size_t sample_idx);
  void computePredictionErrorInternal();
  void writeOutputInternal();
  void writeConfusionFile();
  void writePredictionFile();
  void saveToFileInternal(std::ofstream& outfile);
  void loadFromFileInternal(std::ifstream& infile);

  // Classes of the dependent variable and classIDs for responses
  std::vector<double> class_values;
  std::vector<uint> response_classIDs;
  std::vector<std::vector<size_t>> sampleIDs_per_class;

  // Splitting weights
  std::vector<double> class_weights;

  // Table with classifications and true classes
  std::map<std::pair<double, double>, size_t> classification_table;
};

#endif /* FORESTPROBABILITY_H_ */
