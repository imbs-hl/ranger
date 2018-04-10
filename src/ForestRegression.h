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
  ForestRegression();
  virtual ~ForestRegression();

  void loadForest(size_t dependent_varID, size_t num_trees,
      std::vector<std::vector<std::vector<size_t>> >& forest_child_nodeIDs,
      std::vector<std::vector<size_t>>& forest_split_varIDs, std::vector<std::vector<double>>& forest_split_values,
      std::vector<bool>& is_ordered_variable);

private:
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

  DISALLOW_COPY_AND_ASSIGN(ForestRegression);
};

} // namespace ranger

#endif /* FORESTREGRESSION_H_ */
