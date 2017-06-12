/*-------------------------------------------------------------------------------
 This file is part of Ranger.

 Ranger is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Ranger is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with Ranger. If not, see <http://www.gnu.org/licenses/>.

 Written by:

 Marvin N. Wright
 Institut für Medizinische Biometrie und Statistik
 Universität zu Lübeck
 Ratzeburger Allee 160
 23562 Lübeck
 Germany

 http://www.imbs-luebeck.de
 #-------------------------------------------------------------------------------*/

#include "TreeProbability.h"
#include "utility.h"
#include "Data.h"

TreeProbability::TreeProbability(std::vector<double>* class_values, std::vector<uint>* response_classIDs) :
    class_values(class_values), response_classIDs(response_classIDs), counter(0), counter_per_class(0) {
}

TreeProbability::TreeProbability(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
    std::vector<double>& split_values, std::vector<double>* class_values, std::vector<uint>* response_classIDs,
    std::vector<std::vector<double>>& terminal_class_counts, std::vector<bool>* is_ordered_variable) :
    Tree(child_nodeIDs, split_varIDs, split_values, is_ordered_variable), class_values(class_values), response_classIDs(
        response_classIDs), terminal_class_counts(terminal_class_counts), counter(0), counter_per_class(0) {
}

TreeProbability::~TreeProbability() {
  // Empty on purpose
}

void TreeProbability::initInternal() {
  // Init counters if not in memory efficient mode
  if (!memory_saving_splitting) {
    size_t num_classes = class_values->size();
    size_t max_num_splits = data->getMaxNumUniqueValues();

    // Use number of random splits for extratrees
    if (splitrule == EXTRATREES && num_random_splits > max_num_splits) {
      max_num_splits = num_random_splits;
    }

    counter = new size_t[max_num_splits];
    counter_per_class = new size_t[num_classes * max_num_splits];
  }
}

void TreeProbability::addToTerminalNodes(size_t nodeID) {

  size_t num_samples_in_node = sampleIDs[nodeID].size();
  terminal_class_counts[nodeID].resize(class_values->size(), 0);

  // Compute counts
  for (size_t i = 0; i < num_samples_in_node; ++i) {
    size_t node_sampleID = sampleIDs[nodeID][i];
    size_t classID = (*response_classIDs)[node_sampleID];
    ++terminal_class_counts[nodeID][classID];
  }

  // Compute fractions
  for (size_t i = 0; i < terminal_class_counts[nodeID].size(); ++i) {
    terminal_class_counts[nodeID][i] /= num_samples_in_node;
  }
}

void TreeProbability::appendToFileInternal(std::ofstream& file) { // #nocov start

  // Add Terminal node class counts
  // Convert to vector without empty elements and save
  std::vector<size_t> terminal_nodes;
  std::vector<std::vector<double>> terminal_class_counts_vector;
  for (size_t i = 0; i < terminal_class_counts.size(); ++i) {
    if (!terminal_class_counts[i].empty()) {
      terminal_nodes.push_back(i);
      terminal_class_counts_vector.push_back(terminal_class_counts[i]);
    }
  }
  saveVector1D(terminal_nodes, file);
  saveVector2D(terminal_class_counts_vector, file);
} // #nocov end

bool TreeProbability::splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  // Check node size, stop if maximum reached
  if (sampleIDs[nodeID].size() <= min_node_size) {
    addToTerminalNodes(nodeID);
    return true;
  }

  // Check if node is pure and set split_value to estimate and stop if pure
  bool pure = true;
  double pure_value = 0;
  for (size_t i = 0; i < sampleIDs[nodeID].size(); ++i) {
    double value = data->get(sampleIDs[nodeID][i], dependent_varID);
    if (i != 0 && value != pure_value) {
      pure = false;
      break;
    }
    pure_value = value;
  }
  if (pure) {
    addToTerminalNodes(nodeID);
    return true;
  }

  // Find best split, stop if no decrease of impurity
  bool stop;
  if (splitrule == EXTRATREES) {
    stop = findBestSplitExtraTrees(nodeID, possible_split_varIDs);
  } else {
    stop = findBestSplit(nodeID, possible_split_varIDs);
  }

  if (stop) {
    addToTerminalNodes(nodeID);
    return true;
  }

  return false;
}

void TreeProbability::createEmptyNodeInternal() {
  terminal_class_counts.push_back(std::vector<double>());
}

double TreeProbability::computePredictionAccuracyInternal() {

  size_t num_predictions = prediction_terminal_nodeIDs.size();
  double sum_of_squares = 0;
  for (size_t i = 0; i < num_predictions; ++i) {
    size_t sampleID = oob_sampleIDs[i];
    size_t real_classID = (*response_classIDs)[sampleID];
    size_t terminal_nodeID = prediction_terminal_nodeIDs[i];
    double predicted_value = terminal_class_counts[terminal_nodeID][real_classID];
    sum_of_squares += (1 - predicted_value) * (1 - predicted_value);
  }
  return (1.0 - sum_of_squares / (double) num_predictions);
}

bool TreeProbability::findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  size_t num_samples_node = sampleIDs[nodeID].size();
  size_t num_classes = class_values->size();
  double best_decrease = -1;
  size_t best_varID = 0;
  double best_value = 0;

  size_t* class_counts = new size_t[num_classes]();
  // Compute overall class counts
  for (size_t i = 0; i < num_samples_node; ++i) {
    size_t sampleID = sampleIDs[nodeID][i];
    uint sample_classID = (*response_classIDs)[sampleID];
    ++class_counts[sample_classID];
  }

  // For all possible split variables
  for (auto& varID : possible_split_varIDs) {
    // Find best split value, if ordered consider all values as split values, else all 2-partitions
    if ((*is_ordered_variable)[varID]) {

      // Use memory saving method if option set
      if (memory_saving_splitting) {
        findBestSplitValueSmallQ(nodeID, varID, num_classes, class_counts, num_samples_node, best_value, best_varID,
            best_decrease);
      } else {
        // Use faster method for both cases
        double q = (double) num_samples_node / (double) data->getNumUniqueDataValues(varID);
        if (q < Q_THRESHOLD) {
          findBestSplitValueSmallQ(nodeID, varID, num_classes, class_counts, num_samples_node, best_value, best_varID,
              best_decrease);
        } else {
          findBestSplitValueLargeQ(nodeID, varID, num_classes, class_counts, num_samples_node, best_value, best_varID,
              best_decrease);
        }
      }
    } else {
      findBestSplitValueUnordered(nodeID, varID, num_classes, class_counts, num_samples_node, best_value, best_varID,
          best_decrease);
    }
  }

  delete[] class_counts;

  // Stop if no good split found
  if (best_decrease < 0) {
    return true;
  }

  // Save best values
  split_varIDs[nodeID] = best_varID;
  split_values[nodeID] = best_value;

  // Compute decrease of impurity for this node and add to variable importance if needed
  if (importance_mode == IMP_GINI) {
    addImpurityImportance(nodeID, best_varID, best_decrease);
  }
  return false;
}

void TreeProbability::findBestSplitValueSmallQ(size_t nodeID, size_t varID, size_t num_classes, size_t* class_counts,
    size_t num_samples_node, double& best_value, size_t& best_varID, double& best_decrease) {

  // Create possible split values
  std::vector<double> possible_split_values;
  data->getAllValues(possible_split_values, sampleIDs[nodeID], varID);

  // Try next variable if all equal for this
  if (possible_split_values.size() < 2) {
    return;
  }

  // Initialize with 0, if not in memory efficient mode, use pre-allocated space
  // -1 because no split possible at largest value
  size_t num_splits = possible_split_values.size() - 1;
  size_t* class_counts_right;
  size_t* n_right;
  if (memory_saving_splitting) {
    class_counts_right = new size_t[num_splits * num_classes]();
    n_right = new size_t[num_splits]();
  } else {
    class_counts_right = counter_per_class;
    n_right = counter;
    std::fill(class_counts_right, class_counts_right + num_splits * num_classes, 0);
    std::fill(n_right, n_right + num_splits, 0);
  }

  // Count samples in right child per class and possbile split
  for (auto& sampleID : sampleIDs[nodeID]) {
    double value = data->get(sampleID, varID);
    uint sample_classID = (*response_classIDs)[sampleID];

    // Count samples until split_value reached
    for (size_t i = 0; i < num_splits; ++i) {
      if (value > possible_split_values[i]) {
        ++n_right[i];
        ++class_counts_right[i * num_classes + sample_classID];
      } else {
        break;
      }
    }
  }

  // Compute decrease of impurity for each possible split
  for (size_t i = 0; i < num_splits; ++i) {

    // Stop if one child empty
    size_t n_left = num_samples_node - n_right[i];
    if (n_left == 0 || n_right[i] == 0) {
      continue;
    }

    // Sum of squares
    double sum_left = 0;
    double sum_right = 0;
    for (size_t j = 0; j < num_classes; ++j) {
      size_t class_count_right = class_counts_right[i * num_classes + j];
      size_t class_count_left = class_counts[j] - class_count_right;

      sum_right += class_count_right * class_count_right;
      sum_left += class_count_left * class_count_left;
    }

    // Decrease of impurity
    double decrease = sum_left / (double) n_left + sum_right / (double) n_right[i];

    // If better than before, use this
    if (decrease > best_decrease) {
      best_value = (possible_split_values[i] + possible_split_values[i + 1]) / 2;
      best_varID = varID;
      best_decrease = decrease;

      // Use smaller value if average is numerically the same as the larger value
      if (best_value == possible_split_values[i + 1]) {
        best_value = possible_split_values[i];
      }
    }
  }

  if (memory_saving_splitting) {
    delete[] class_counts_right;
    delete[] n_right;
  }
}

void TreeProbability::findBestSplitValueLargeQ(size_t nodeID, size_t varID, size_t num_classes, size_t* class_counts,
    size_t num_samples_node, double& best_value, size_t& best_varID, double& best_decrease) {

  // Set counters to 0
  size_t num_unique = data->getNumUniqueDataValues(varID);
  std::fill(counter_per_class, counter_per_class + num_unique * num_classes, 0);
  std::fill(counter, counter + num_unique, 0);

  // Count values
  for (auto& sampleID : sampleIDs[nodeID]) {
    size_t index = data->getIndex(sampleID, varID);
    size_t classID = (*response_classIDs)[sampleID];

    ++counter[index];
    ++counter_per_class[index * num_classes + classID];
  }

  size_t n_left = 0;
  size_t* class_counts_left = new size_t[num_classes]();

  // Compute decrease of impurity for each split
  for (size_t i = 0; i < num_unique - 1; ++i) {

    // Stop if nothing here
    if (counter[i] == 0) {
      continue;
    }

    n_left += counter[i];

    // Stop if right child empty
    size_t n_right = num_samples_node - n_left;
    if (n_right == 0) {
      break;
    }

    // Sum of squares
    double sum_left = 0;
    double sum_right = 0;
    for (size_t j = 0; j < num_classes; ++j) {
      class_counts_left[j] += counter_per_class[i * num_classes + j];
      size_t class_count_right = class_counts[j] - class_counts_left[j];

      sum_left += class_counts_left[j] * class_counts_left[j];
      sum_right += class_count_right * class_count_right;
    }

    // Decrease of impurity
    double decrease = sum_right / (double) n_right + sum_left / (double) n_left;

    // If better than before, use this
    if (decrease > best_decrease) {
      // Find next value in this node
      size_t j = i + 1;
      while(j < num_unique && counter[j] == 0) {
        ++j;
      }

      // Use mid-point split
      best_value = (data->getUniqueDataValue(varID, i) + data->getUniqueDataValue(varID, j)) / 2;
      best_varID = varID;
      best_decrease = decrease;

      // Use smaller value if average is numerically the same as the larger value
      if (best_value == data->getUniqueDataValue(varID, j)) {
        best_value = data->getUniqueDataValue(varID, i);
      }
    }
  }

  delete[] class_counts_left;
}

void TreeProbability::findBestSplitValueUnordered(size_t nodeID, size_t varID, size_t num_classes, size_t* class_counts,
    size_t num_samples_node, double& best_value, size_t& best_varID, double& best_decrease) {

  // Create possible split values
  std::vector<double> factor_levels;
  data->getAllValues(factor_levels, sampleIDs[nodeID], varID);

  // Try next variable if all equal for this
  if (factor_levels.size() < 2) {
    return;
  }

  // Number of possible splits is 2^num_levels
  size_t num_splits = (1 << factor_levels.size());

  // Compute decrease of impurity for each possible split
  // Split where all left (0) or all right (1) are excluded
  // The second half of numbers is just left/right switched the first half -> Exclude second half
  for (size_t local_splitID = 1; local_splitID < num_splits / 2; ++local_splitID) {

    // Compute overall splitID by shifting local factorIDs to global positions
    size_t splitID = 0;
    for (size_t j = 0; j < factor_levels.size(); ++j) {
      if ((local_splitID & (1 << j))) {
        double level = factor_levels[j];
        size_t factorID = floor(level) - 1;
        splitID = splitID | (1 << factorID);
      }
    }

    // Initialize
    size_t* class_counts_right = new size_t[num_classes]();
    size_t n_right = 0;

    // Count classes in left and right child
    for (auto& sampleID : sampleIDs[nodeID]) {
      uint sample_classID = (*response_classIDs)[sampleID];
      double value = data->get(sampleID, varID);
      size_t factorID = floor(value) - 1;

      // If in right child, count
      // In right child, if bitwise splitID at position factorID is 1
      if ((splitID & (1 << factorID))) {
        ++n_right;
        ++class_counts_right[sample_classID];
      }
    }
    size_t n_left = num_samples_node - n_right;

    // Sum of squares
    double sum_left = 0;
    double sum_right = 0;
    for (size_t j = 0; j < num_classes; ++j) {
      size_t class_count_right = class_counts_right[j];
      size_t class_count_left = class_counts[j] - class_count_right;

      sum_right += class_count_right * class_count_right;
      sum_left += class_count_left * class_count_left;
    }

    // Decrease of impurity
    double decrease = sum_left / (double) n_left + sum_right / (double) n_right;

    // If better than before, use this
    if (decrease > best_decrease) {
      best_value = splitID;
      best_varID = varID;
      best_decrease = decrease;
    }
  }
}

bool TreeProbability::findBestSplitExtraTrees(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  size_t num_samples_node = sampleIDs[nodeID].size();
  size_t num_classes = class_values->size();
  double best_decrease = -1;
  size_t best_varID = 0;
  double best_value = 0;

  size_t* class_counts = new size_t[num_classes]();
  // Compute overall class counts
  for (size_t i = 0; i < num_samples_node; ++i) {
    size_t sampleID = sampleIDs[nodeID][i];
    uint sample_classID = (*response_classIDs)[sampleID];
    ++class_counts[sample_classID];
  }

  // For all possible split variables
  for (auto& varID : possible_split_varIDs) {
    // Find best split value, if ordered consider all values as split values, else all 2-partitions
    if ((*is_ordered_variable)[varID]) {
      findBestSplitValueExtraTrees(nodeID, varID, num_classes, class_counts, num_samples_node, best_value, best_varID,
          best_decrease);
    } else {
      findBestSplitValueExtraTreesUnordered(nodeID, varID, num_classes, class_counts, num_samples_node, best_value,
          best_varID, best_decrease);
    }
  }

  delete[] class_counts;

  // Stop if no good split found
  if (best_decrease < 0) {
    return true;
  }

  // Save best values
  split_varIDs[nodeID] = best_varID;
  split_values[nodeID] = best_value;

  // Compute decrease of impurity for this node and add to variable importance if needed
  if (importance_mode == IMP_GINI) {
    addImpurityImportance(nodeID, best_varID, best_decrease);
  }
  return false;
}

void TreeProbability::findBestSplitValueExtraTrees(size_t nodeID, size_t varID, size_t num_classes,
    size_t* class_counts, size_t num_samples_node, double& best_value, size_t& best_varID, double& best_decrease) {

  // Get min/max values of covariate in node
  double min;
  double max;
  data->getMinMaxValues(min, max, sampleIDs[nodeID], varID);

  // Try next variable if all equal for this
  if (min == max) {
    return;
  }

  // Create possible split values: Draw randomly between min and max
  std::vector<double> possible_split_values;
  std::uniform_real_distribution<double> udist(min, max);
  possible_split_values.reserve(num_random_splits);
  for (size_t i = 0; i < num_random_splits; ++i) {
    possible_split_values.push_back(udist(random_number_generator));
  }

  // Initialize with 0, if not in memory efficient mode, use pre-allocated space
  size_t num_splits = possible_split_values.size();
  size_t* class_counts_right;
  size_t* n_right;
  if (memory_saving_splitting) {
    class_counts_right = new size_t[num_splits * num_classes]();
    n_right = new size_t[num_splits]();
  } else {
    class_counts_right = counter_per_class;
    n_right = counter;
    std::fill(class_counts_right, class_counts_right + num_splits * num_classes, 0);
    std::fill(n_right, n_right + num_splits, 0);
  }

  // Count samples in right child per class and possbile split
  for (auto& sampleID : sampleIDs[nodeID]) {
    double value = data->get(sampleID, varID);
    uint sample_classID = (*response_classIDs)[sampleID];

    // Count samples until split_value reached
    for (size_t i = 0; i < num_splits; ++i) {
      if (value > possible_split_values[i]) {
        ++n_right[i];
        ++class_counts_right[i * num_classes + sample_classID];
      } else {
        break;
      }
    }
  }

  // Compute decrease of impurity for each possible split
  for (size_t i = 0; i < num_splits; ++i) {

    // Stop if one child empty
    size_t n_left = num_samples_node - n_right[i];
    if (n_left == 0 || n_right[i] == 0) {
      continue;
    }

    // Sum of squares
    double sum_left = 0;
    double sum_right = 0;
    for (size_t j = 0; j < num_classes; ++j) {
      size_t class_count_right = class_counts_right[i * num_classes + j];
      size_t class_count_left = class_counts[j] - class_count_right;

      sum_right += class_count_right * class_count_right;
      sum_left += class_count_left * class_count_left;
    }

    // Decrease of impurity
    double decrease = sum_left / (double) n_left + sum_right / (double) n_right[i];

    // If better than before, use this
    if (decrease > best_decrease) {
      best_value = possible_split_values[i];
      best_varID = varID;
      best_decrease = decrease;
    }
  }

  if (memory_saving_splitting) {
    delete[] class_counts_right;
    delete[] n_right;
  }
}

void TreeProbability::findBestSplitValueExtraTreesUnordered(size_t nodeID, size_t varID, size_t num_classes,
    size_t* class_counts, size_t num_samples_node, double& best_value, size_t& best_varID, double& best_decrease) {

  size_t num_unique_values = data->getNumUniqueDataValues(varID);

  // Get all factor indices in node
  std::vector<bool> factor_in_node(num_unique_values, false);
  for (auto& sampleID : sampleIDs[nodeID]) {
    size_t index = data->getIndex(sampleID, varID);
    factor_in_node[index] = true;
  }

  // Vector of indices in and out of node
  std::vector<size_t> indices_in_node;
  std::vector<size_t> indices_out_node;
  indices_in_node.reserve(num_unique_values);
  indices_out_node.reserve(num_unique_values);
  for (size_t i = 0; i < num_unique_values; ++i) {
    if (factor_in_node[i]) {
      indices_in_node.push_back(i);
    } else {
      indices_out_node.push_back(i);
    }
  }

  // Generate num_random_splits splits
  for (size_t i = 0; i < num_random_splits; ++i) {
    std::vector<size_t> split_subset;
    split_subset.reserve(num_unique_values);

    // Draw random subsets, sample all partitions with equal probability
    if (indices_in_node.size() > 1) {
      size_t num_partitions = (2 << (indices_in_node.size() - 1)) - 2; // 2^n-2 (don't allow full or empty)
      std::uniform_int_distribution<size_t> udist(1, num_partitions);
      size_t splitID_in_node = udist(random_number_generator);
      for (size_t j = 0; j < indices_in_node.size(); ++j) {
        if ((splitID_in_node & (1 << j)) > 0) {
          split_subset.push_back(indices_in_node[j]);
        }
      }
    }
    if (indices_out_node.size() > 1) {
      size_t num_partitions = (2 << (indices_out_node.size() - 1)) - 1; // 2^n-1 (allow full or empty)
      std::uniform_int_distribution<size_t> udist(0, num_partitions);
      size_t splitID_out_node = udist(random_number_generator);
      for (size_t j = 0; j < indices_out_node.size(); ++j) {
        if ((splitID_out_node & (1 << j)) > 0) {
          split_subset.push_back(indices_out_node[j]);
        }
      }
    }

    // Assign union of the two subsets to right child
    size_t splitID = 0;
    for (auto& idx : split_subset) {
      splitID |= 1 << idx;
    }

    // Initialize
    size_t* class_counts_right = new size_t[num_classes]();
    size_t n_right = 0;

    // Count classes in left and right child
    for (auto& sampleID : sampleIDs[nodeID]) {
      uint sample_classID = (*response_classIDs)[sampleID];
      double value = data->get(sampleID, varID);
      size_t factorID = floor(value) - 1;

      // If in right child, count
      // In right child, if bitwise splitID at position factorID is 1
      if ((splitID & (1 << factorID))) {
        ++n_right;
        ++class_counts_right[sample_classID];
      }
    }
    size_t n_left = num_samples_node - n_right;

    // Sum of squares
    double sum_left = 0;
    double sum_right = 0;
    for (size_t j = 0; j < num_classes; ++j) {
      size_t class_count_right = class_counts_right[j];
      size_t class_count_left = class_counts[j] - class_count_right;

      sum_right += class_count_right * class_count_right;
      sum_left += class_count_left * class_count_left;
    }

    // Decrease of impurity
    double decrease = sum_left / (double) n_left + sum_right / (double) n_right;

    // If better than before, use this
    if (decrease > best_decrease) {
      best_value = splitID;
      best_varID = varID;
      best_decrease = decrease;
    }
  }
}

void TreeProbability::addImpurityImportance(size_t nodeID, size_t varID, double decrease) {

  std::vector<size_t> class_counts;
  class_counts.resize(class_values->size(), 0);

  for (auto& sampleID : sampleIDs[nodeID]) {
    uint sample_classID = (*response_classIDs)[sampleID];
    class_counts[sample_classID]++;
  }
  double sum_node = 0;
  for (auto& class_count : class_counts) {
    sum_node += class_count * class_count;
  }
  double best_gini = decrease - sum_node / (double) sampleIDs[nodeID].size();

// No variable importance for no split variables
  size_t tempvarID = varID;
  for (auto& skip : *no_split_variables) {
    if (varID >= skip) {
      --tempvarID;
    }
  }
  (*variable_importance)[tempvarID] += best_gini;
}

