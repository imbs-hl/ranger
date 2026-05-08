/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#ifndef FORESTREGRESSION_H_
#define FORESTREGRESSION_H_

#include <iostream>
#include <vector>

#include "globals.h"
#include "Forest.h"

namespace ranger {

class ForestRegression: public Forest {
public:
  ForestRegression() = default;

  ForestRegression(const ForestRegression&) = delete;
  ForestRegression& operator=(const ForestRegression&) = delete;

  virtual ~ForestRegression() override = default;

  void loadForest(size_t num_trees, std::vector<std::vector<std::vector<size_t>> >& forest_child_nodeIDs,
      std::vector<std::vector<size_t>>& forest_split_varIDs, std::vector<std::vector<double>>& forest_split_values,
      std::vector<bool>& is_ordered_variable);

private:
  void initInternal() override;
  void growInternal() override;
  void allocatePredictMemory() override;
  void predictInternal(size_t sample_idx) override;
  void computePredictionErrorInternal() override;
  void writeOutputInternal() override;
  void writeConfusionFile() override;
  void writePredictionFile() override;
  void saveToFileInternal(std::ofstream& outfile) override;
  void loadFromFileInternal(std::ifstream& infile) override;

private:
  double getTreePrediction(size_t tree_idx, size_t sample_idx) const;
  size_t getTreePredictionTerminalNodeID(size_t tree_idx, size_t sample_idx) const;
};

} // namespace ranger

#endif /* FORESTREGRESSION_H_ */
