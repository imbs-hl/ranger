/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#include "TreeProbability.h"
#include "utility.h"
#include "Data.h"

namespace ranger {

TreeProbability::TreeProbability(std::vector<double>* class_values, std::vector<uint>* response_classIDs,
    std::vector<std::vector<size_t>>* sampleIDs_per_class, std::vector<double>* class_weights) :
    class_values(class_values), response_classIDs(response_classIDs), sampleIDs_per_class(sampleIDs_per_class), class_weights(
        class_weights), counter(0), counter_per_class(0) {
}

TreeProbability::TreeProbability(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
    std::vector<double>& split_values, std::vector<double>* class_values, std::vector<uint>* response_classIDs,
    std::vector<std::vector<double>>& terminal_class_counts) :
    Tree(child_nodeIDs, split_varIDs, split_values), class_values(class_values), response_classIDs(response_classIDs), sampleIDs_per_class(
        0), terminal_class_counts(terminal_class_counts), class_weights(0), counter(0), counter_per_class(0) {
}

void TreeProbability::allocateMemory() {
  // Init counters if not in memory efficient mode
  if (!memory_saving_splitting) {
    size_t num_classes = class_values->size();
    size_t max_num_splits = data->getMaxNumUniqueValues();

    // Use number of random splits for extratrees
    if (splitrule == EXTRATREES && num_random_splits > max_num_splits) {
      max_num_splits = num_random_splits;
    }

    counter.resize(max_num_splits);
    counter_per_class.resize(num_classes * max_num_splits);
  }
}

void TreeProbability::addToTerminalNodes(size_t nodeID) {

  size_t num_samples_in_node = end_pos[nodeID] - start_pos[nodeID];
  terminal_class_counts[nodeID].resize(class_values->size(), 0);

  // Compute counts
  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    size_t classID = (*response_classIDs)[sampleID];
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

  // Stop if maximum node size or depth reached
  size_t num_samples_node = end_pos[nodeID] - start_pos[nodeID];
  if (num_samples_node <= min_node_size || (nodeID >= last_left_nodeID && max_depth > 0 && depth >= max_depth)) {
    addToTerminalNodes(nodeID);
    return true;
  }

  // Check if node is pure and set split_value to estimate and stop if pure
  bool pure = true;
  double pure_value = 0;
  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    double value = data->get_y(sampleID, 0);
    if (pos != start_pos[nodeID] && value != pure_value) {
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

double TreeProbability::computePredictionAccuracyInternal(std::vector<double>* prediction_error_casewise) {

  size_t num_predictions = prediction_terminal_nodeIDs.size();
  double sum_of_squares = 0;
  for (size_t i = 0; i < num_predictions; ++i) {
    size_t sampleID = oob_sampleIDs[i];
    size_t real_classID = (*response_classIDs)[sampleID];
    size_t terminal_nodeID = prediction_terminal_nodeIDs[i];
    double predicted_value = terminal_class_counts[terminal_nodeID][real_classID];
    double err = (1 - predicted_value) * (1 - predicted_value);
    if (prediction_error_casewise) {
      (*prediction_error_casewise)[i] = err;
    }
    sum_of_squares += err;
  }
  return (1.0 - sum_of_squares / (double) num_predictions);
}

bool TreeProbability::findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  size_t num_samples_node = end_pos[nodeID] - start_pos[nodeID];
  size_t num_classes = class_values->size();
  double best_decrease = -1;
  size_t best_varID = 0;
  double best_value = 0;

  std::vector<size_t> class_counts(num_classes);
  // Compute overall class counts
  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    uint sample_classID = (*response_classIDs)[sampleID];
    ++class_counts[sample_classID];
  }

  // Stop early if no split posssible
  if (num_samples_node >= 2 * min_bucket) {

    // For all possible split variables
    for (auto& varID : possible_split_varIDs) {
      // Find best split value, if ordered consider all values as split values, else all 2-partitions
      if (data->isOrderedVariable(varID)) {

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
  }

  // Stop if no good split found
  if (best_decrease < 0) {
    return true;
  }

  // Save best values
  split_varIDs[nodeID] = best_varID;
  split_values[nodeID] = best_value;

  // Compute decrease of impurity for this node and add to variable importance if needed
  if (importance_mode == IMP_GINI || importance_mode == IMP_GINI_CORRECTED) {
    addImpurityImportance(nodeID, best_varID, best_decrease);
  }

  // Regularization
  saveSplitVarID(best_varID);

  return false;
}

void TreeProbability::findBestSplitValueSmallQ(size_t nodeID, size_t varID, size_t num_classes,
    const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
    double& best_decrease) {

  // Create possible split values
  std::vector<double> possible_split_values;
  data->getAllValues(possible_split_values, sampleIDs, varID, start_pos[nodeID], end_pos[nodeID]);

  // Try next variable if all equal for this
  if (possible_split_values.size() < 2) {
    return;
  }

  const size_t num_splits = possible_split_values.size();
  if (memory_saving_splitting) {
    std::vector<size_t> class_counts_right(num_splits * num_classes), n_right(num_splits);
    findBestSplitValueSmallQ(nodeID, varID, num_classes, class_counts, num_samples_node, best_value, best_varID,
        best_decrease, possible_split_values, class_counts_right, n_right);
  } else {
    std::fill_n(counter_per_class.begin(), num_splits * num_classes, 0);
    std::fill_n(counter.begin(), num_splits, 0);
    findBestSplitValueSmallQ(nodeID, varID, num_classes, class_counts, num_samples_node, best_value, best_varID,
        best_decrease, possible_split_values, counter_per_class, counter);
  }
}

void TreeProbability::findBestSplitValueSmallQ(size_t nodeID, size_t varID, size_t num_classes,
    const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
    double& best_decrease, const std::vector<double>& possible_split_values, std::vector<size_t>& counter_per_class,
    std::vector<size_t>& counter) {

  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    uint sample_classID = (*response_classIDs)[sampleID];
    size_t idx = std::lower_bound(possible_split_values.begin(), possible_split_values.end(),
        data->get_x(sampleID, varID)) - possible_split_values.begin();

    ++counter_per_class[idx * num_classes + sample_classID];
    ++counter[idx];
  }

  size_t n_left = 0;
  std::vector<size_t> class_counts_left(num_classes);

  // Compute decrease of impurity for each split
  for (size_t i = 0; i < possible_split_values.size() - 1; ++i) {

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

    // Stop if minimal bucket size reached
    if (n_left < min_bucket || n_right < min_bucket) {
      continue;
    }

    double decrease;
    if (splitrule == HELLINGER) {
      for (size_t j = 0; j < num_classes; ++j) {
        class_counts_left[j] += counter_per_class[i * num_classes + j];
      }

      // TPR is number of outcome 1s in one node / total number of 1s
      // FPR is number of outcome 0s in one node / total number of 0s
      double tpr = (double) (class_counts[1] - class_counts_left[1]) / (double) class_counts[1];
      double fpr = (double) (class_counts[0] - class_counts_left[0]) / (double) class_counts[0];

      // Decrease of impurity
      double a1 = sqrt(tpr) - sqrt(fpr);
      double a2 = sqrt(1 - tpr) - sqrt(1 - fpr);
      decrease = sqrt(a1 * a1 + a2 * a2);
    } else {
      // Sum of squares
      double sum_left = 0;
      double sum_right = 0;
      for (size_t j = 0; j < num_classes; ++j) {
        class_counts_left[j] += counter_per_class[i * num_classes + j];
        size_t class_count_right = class_counts[j] - class_counts_left[j];

        sum_left += (*class_weights)[j] * class_counts_left[j] * class_counts_left[j];
        sum_right += (*class_weights)[j] * class_count_right * class_count_right;
      }

      // Decrease of impurity
      decrease = sum_right / (double) n_right + sum_left / (double) n_left;
    }

    // Regularization
    regularize(decrease, varID);

    // If better than before, use this
    if (decrease > best_decrease) {
      // Use mid-point split
      best_value = (possible_split_values[i] + possible_split_values[i + 1]) / 2;
      best_varID = varID;
      best_decrease = decrease;

      // Use smaller value if average is numerically the same as the larger value
      if (best_value == possible_split_values[i + 1]) {
        best_value = possible_split_values[i];
      }
    }
  }
}

void TreeProbability::findBestSplitValueLargeQ(size_t nodeID, size_t varID, size_t num_classes,
    const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
    double& best_decrease) {

  // Set counters to 0
  size_t num_unique = data->getNumUniqueDataValues(varID);
  std::fill_n(counter_per_class.begin(), num_unique * num_classes, 0);
  std::fill_n(counter.begin(), num_unique, 0);

  // Count values
  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    size_t index = data->getIndex(sampleID, varID);
    size_t classID = (*response_classIDs)[sampleID];

    ++counter[index];
    ++counter_per_class[index * num_classes + classID];
  }

  size_t n_left = 0;
  std::vector<size_t> class_counts_left(num_classes);

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

    // Stop if minimal bucket size reached
    if (n_left < min_bucket || n_right < min_bucket) {
      continue;
    }

    double decrease;
    if (splitrule == HELLINGER) {
      for (size_t j = 0; j < num_classes; ++j) {
        class_counts_left[j] += counter_per_class[i * num_classes + j];
      }

      // TPR is number of outcome 1s in one node / total number of 1s
      // FPR is number of outcome 0s in one node / total number of 0s
      double tpr = (double) (class_counts[1] - class_counts_left[1]) / (double) class_counts[1];
      double fpr = (double) (class_counts[0] - class_counts_left[0]) / (double) class_counts[0];

      // Decrease of impurity
      double a1 = sqrt(tpr) - sqrt(fpr);
      double a2 = sqrt(1 - tpr) - sqrt(1 - fpr);
      decrease = sqrt(a1 * a1 + a2 * a2);
    } else {
      // Sum of squares
      double sum_left = 0;
      double sum_right = 0;
      for (size_t j = 0; j < num_classes; ++j) {
        class_counts_left[j] += counter_per_class[i * num_classes + j];
        size_t class_count_right = class_counts[j] - class_counts_left[j];

        sum_left += (*class_weights)[j] * class_counts_left[j] * class_counts_left[j];
        sum_right += (*class_weights)[j] * class_count_right * class_count_right;
      }

      // Decrease of impurity
      decrease = sum_right / (double) n_right + sum_left / (double) n_left;
    }

    // Regularization
    regularize(decrease, varID);

    // If better than before, use this
    if (decrease > best_decrease) {
      // Find next value in this node
      size_t j = i + 1;
      while (j < num_unique && counter[j] == 0) {
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
}

void TreeProbability::findBestSplitValueUnordered(size_t nodeID, size_t varID, size_t num_classes,
    const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
    double& best_decrease) {

  // Create possible split values
  std::vector<double> factor_levels;
  data->getAllValues(factor_levels, sampleIDs, varID, start_pos[nodeID], end_pos[nodeID]);

  // Try next variable if all equal for this
  if (factor_levels.size() < 2) {
    return;
  }

  // Number of possible splits is 2^num_levels
  size_t num_splits = (1ULL << factor_levels.size());

  // Compute decrease of impurity for each possible split
  // Split where all left (0) or all right (1) are excluded
  // The second half of numbers is just left/right switched the first half -> Exclude second half
  for (size_t local_splitID = 1; local_splitID < num_splits / 2; ++local_splitID) {

    // Compute overall splitID by shifting local factorIDs to global positions
    size_t splitID = 0;
    for (size_t j = 0; j < factor_levels.size(); ++j) {
      if ((local_splitID & (1ULL << j))) {
        double level = factor_levels[j];
        size_t factorID = floor(level) - 1;
        splitID = splitID | (1ULL << factorID);
      }
    }

    // Initialize
    std::vector<size_t> class_counts_right(num_classes);
    size_t n_right = 0;

    // Count classes in left and right child
    for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
      size_t sampleID = sampleIDs[pos];
      uint sample_classID = (*response_classIDs)[sampleID];
      double value = data->get_x(sampleID, varID);
      size_t factorID = floor(value) - 1;

      // If in right child, count
      // In right child, if bitwise splitID at position factorID is 1
      if ((splitID & (1ULL << factorID))) {
        ++n_right;
        ++class_counts_right[sample_classID];
      }
    }
    size_t n_left = num_samples_node - n_right;

    // Stop if minimal bucket size reached
    if (n_left < min_bucket || n_right < min_bucket) {
      continue;
    }

    double decrease;
    if (splitrule == HELLINGER) {
      // TPR is number of outcome 1s in one node / total number of 1s
      // FPR is number of outcome 0s in one node / total number of 0s
      double tpr = (double) class_counts_right[1] / (double) class_counts[1];
      double fpr = (double) class_counts_right[0] / (double) class_counts[0];

      // Decrease of impurity
      double a1 = sqrt(tpr) - sqrt(fpr);
      double a2 = sqrt(1 - tpr) - sqrt(1 - fpr);
      decrease = sqrt(a1 * a1 + a2 * a2);
    } else {
      // Sum of squares
      double sum_left = 0;
      double sum_right = 0;
      for (size_t j = 0; j < num_classes; ++j) {
        size_t class_count_right = class_counts_right[j];
        size_t class_count_left = class_counts[j] - class_count_right;

        sum_right += (*class_weights)[j] * class_count_right * class_count_right;
        sum_left += (*class_weights)[j] * class_count_left * class_count_left;
      }

      // Decrease of impurity
      decrease = sum_left / (double) n_left + sum_right / (double) n_right;
    }

    // Regularization
    regularize(decrease, varID);

    // If better than before, use this
    if (decrease > best_decrease) {
      best_value = splitID;
      best_varID = varID;
      best_decrease = decrease;
    }
  }
}

bool TreeProbability::findBestSplitExtraTrees(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  size_t num_samples_node = end_pos[nodeID] - start_pos[nodeID];
  size_t num_classes = class_values->size();
  double best_decrease = -1;
  size_t best_varID = 0;
  double best_value = 0;

  std::vector<size_t> class_counts(num_classes);
  // Compute overall class counts
  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    uint sample_classID = (*response_classIDs)[sampleID];
    ++class_counts[sample_classID];
  }

  // Stop early if no split posssible
  if (num_samples_node >= 2 * min_bucket) {

    // For all possible split variables
    for (auto& varID : possible_split_varIDs) {
      // Find best split value, if ordered consider all values as split values, else all 2-partitions
      if (data->isOrderedVariable(varID)) {
        findBestSplitValueExtraTrees(nodeID, varID, num_classes, class_counts, num_samples_node, best_value, best_varID,
            best_decrease);
      } else {
        findBestSplitValueExtraTreesUnordered(nodeID, varID, num_classes, class_counts, num_samples_node, best_value,
            best_varID, best_decrease);
      }
    }
  }

  // Stop if no good split found
  if (best_decrease < 0) {
    return true;
  }

  // Save best values
  split_varIDs[nodeID] = best_varID;
  split_values[nodeID] = best_value;

  // Compute decrease of impurity for this node and add to variable importance if needed
  if (importance_mode == IMP_GINI || importance_mode == IMP_GINI_CORRECTED) {
    addImpurityImportance(nodeID, best_varID, best_decrease);
  }

  // Regularization
  saveSplitVarID(best_varID);

  return false;
}

void TreeProbability::findBestSplitValueExtraTrees(size_t nodeID, size_t varID, size_t num_classes,
    const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
    double& best_decrease) {

  // Get min/max values of covariate in node
  double min;
  double max;
  data->getMinMaxValues(min, max, sampleIDs, varID, start_pos[nodeID], end_pos[nodeID]);

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
  if (num_random_splits > 1) {
    std::sort(possible_split_values.begin(), possible_split_values.end());
  }

  const size_t num_splits = possible_split_values.size();
  if (memory_saving_splitting) {
    std::vector<size_t> class_counts_right(num_splits * num_classes), n_right(num_splits);
    findBestSplitValueExtraTrees(nodeID, varID, num_classes, class_counts, num_samples_node, best_value, best_varID,
        best_decrease, possible_split_values, class_counts_right, n_right);
  } else {
    std::fill_n(counter_per_class.begin(), num_splits * num_classes, 0);
    std::fill_n(counter.begin(), num_splits, 0);
    findBestSplitValueExtraTrees(nodeID, varID, num_classes, class_counts, num_samples_node, best_value, best_varID,
        best_decrease, possible_split_values, counter_per_class, counter);
  }
}

void TreeProbability::findBestSplitValueExtraTrees(size_t nodeID, size_t varID, size_t num_classes,
    const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
    double& best_decrease, const std::vector<double>& possible_split_values, std::vector<size_t>& class_counts_right,
    std::vector<size_t>& n_right) {
  const size_t num_splits = possible_split_values.size();

  // Count samples in right child per class and possbile split
  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    double value = data->get_x(sampleID, varID);
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

    // Stop if minimal bucket size reached
    if (n_left < min_bucket || n_right[i] < min_bucket) {
      continue;
    }

    // Sum of squares
    double sum_left = 0;
    double sum_right = 0;
    for (size_t j = 0; j < num_classes; ++j) {
      size_t class_count_right = class_counts_right[i * num_classes + j];
      size_t class_count_left = class_counts[j] - class_count_right;

      sum_right += (*class_weights)[j] * class_count_right * class_count_right;
      sum_left += (*class_weights)[j] * class_count_left * class_count_left;
    }

    // Decrease of impurity
    double decrease = sum_left / (double) n_left + sum_right / (double) n_right[i];

    // Regularization
    regularize(decrease, varID);

    // If better than before, use this
    if (decrease > best_decrease) {
      best_value = possible_split_values[i];
      best_varID = varID;
      best_decrease = decrease;
    }
  }
}

void TreeProbability::findBestSplitValueExtraTreesUnordered(size_t nodeID, size_t varID, size_t num_classes,
    const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
    double& best_decrease) {

  size_t num_unique_values = data->getNumUniqueDataValues(varID);

  // Get all factor indices in node
  std::vector<bool> factor_in_node(num_unique_values, false);
  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
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
      size_t num_partitions = (2ULL << (indices_in_node.size() - 1ULL)) - 2ULL; // 2^n-2 (don't allow full or empty)
      std::uniform_int_distribution<size_t> udist(1, num_partitions);
      size_t splitID_in_node = udist(random_number_generator);
      for (size_t j = 0; j < indices_in_node.size(); ++j) {
        if ((splitID_in_node & (1ULL << j)) > 0) {
          split_subset.push_back(indices_in_node[j]);
        }
      }
    }
    if (indices_out_node.size() > 1) {
      size_t num_partitions = (2ULL << (indices_out_node.size() - 1ULL)) - 1ULL; // 2^n-1 (allow full or empty)
      std::uniform_int_distribution<size_t> udist(0, num_partitions);
      size_t splitID_out_node = udist(random_number_generator);
      for (size_t j = 0; j < indices_out_node.size(); ++j) {
        if ((splitID_out_node & (1ULL << j)) > 0) {
          split_subset.push_back(indices_out_node[j]);
        }
      }
    }

    // Assign union of the two subsets to right child
    size_t splitID = 0;
    for (auto& idx : split_subset) {
      splitID |= 1ULL << idx;
    }

    // Initialize
    std::vector<size_t> class_counts_right(num_classes);
    size_t n_right = 0;

    // Count classes in left and right child
    for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
      size_t sampleID = sampleIDs[pos];
      uint sample_classID = (*response_classIDs)[sampleID];
      double value = data->get_x(sampleID, varID);
      size_t factorID = floor(value) - 1;

      // If in right child, count
      // In right child, if bitwise splitID at position factorID is 1
      if ((splitID & (1ULL << factorID))) {
        ++n_right;
        ++class_counts_right[sample_classID];
      }
    }
    size_t n_left = num_samples_node - n_right;

    // Stop if minimal bucket size reached
    if (n_left < min_bucket || n_right < min_bucket) {
      continue;
    }

    // Sum of squares
    double sum_left = 0;
    double sum_right = 0;
    for (size_t j = 0; j < num_classes; ++j) {
      size_t class_count_right = class_counts_right[j];
      size_t class_count_left = class_counts[j] - class_count_right;

      sum_right += (*class_weights)[j] * class_count_right * class_count_right;
      sum_left += (*class_weights)[j] * class_count_left * class_count_left;
    }

    // Decrease of impurity
    double decrease = sum_left / (double) n_left + sum_right / (double) n_right;

    // Regularization
    regularize(decrease, varID);

    // If better than before, use this
    if (decrease > best_decrease) {
      best_value = splitID;
      best_varID = varID;
      best_decrease = decrease;
    }
  }
}

void TreeProbability::addImpurityImportance(size_t nodeID, size_t varID, double decrease) {

  double best_decrease = decrease;
  if (splitrule != HELLINGER) {
    size_t num_samples_node = end_pos[nodeID] - start_pos[nodeID];
    std::vector<size_t> class_counts;
    class_counts.resize(class_values->size(), 0);

    for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
      size_t sampleID = sampleIDs[pos];
      uint sample_classID = (*response_classIDs)[sampleID];
      class_counts[sample_classID]++;
    }
    double sum_node = 0;
    for (size_t i = 0; i < class_counts.size(); ++i) {
      sum_node += (*class_weights)[i] * class_counts[i] * class_counts[i];
    }
    best_decrease = decrease - sum_node / (double) num_samples_node;
  }

  // No variable importance for no split variables
  size_t tempvarID = data->getUnpermutedVarID(varID);

  // Subtract if corrected importance and permuted variable, else add
  if (importance_mode == IMP_GINI_CORRECTED && varID >= data->getNumCols()) {
    (*variable_importance)[tempvarID] -= best_decrease;
  } else {
    (*variable_importance)[tempvarID] += best_decrease;
  }
}

void TreeProbability::bootstrapClassWise() {
  // Number of samples is sum of sample fraction * number of samples
  size_t num_samples_inbag = 0;
  double sum_sample_fraction = 0;
  for (auto& s : *sample_fraction) {
    num_samples_inbag += (size_t) num_samples * s;
    sum_sample_fraction += s;
  }

  // Reserve space, reserve a little more to be save)
  sampleIDs.reserve(num_samples_inbag);
  oob_sampleIDs.reserve(num_samples * (exp(-sum_sample_fraction) + 0.1));

  // Start with all samples OOB
  inbag_counts.resize(num_samples, 0);

  // Draw samples for each class
  for (size_t i = 0; i < sample_fraction->size(); ++i) {
    // Draw samples of class with replacement as inbag and mark as not OOB
    size_t num_samples_class = (*sampleIDs_per_class)[i].size();
    size_t num_samples_inbag_class = round(num_samples * (*sample_fraction)[i]);
    std::uniform_int_distribution<size_t> unif_dist(0, num_samples_class - 1);
    for (size_t s = 0; s < num_samples_inbag_class; ++s) {
      size_t draw = (*sampleIDs_per_class)[i][unif_dist(random_number_generator)];
      sampleIDs.push_back(draw);
      ++inbag_counts[draw];
    }
  }

  // Save OOB samples
  for (size_t s = 0; s < inbag_counts.size(); ++s) {
    if (inbag_counts[s] == 0) {
      oob_sampleIDs.push_back(s);
    }
  }
  num_samples_oob = oob_sampleIDs.size();

  if (!keep_inbag) {
    inbag_counts.clear();
    inbag_counts.shrink_to_fit();
  }
}

void TreeProbability::bootstrapWithoutReplacementClassWise() {
  // Draw samples for each class
  for (size_t i = 0; i < sample_fraction->size(); ++i) {
    size_t num_samples_class = (*sampleIDs_per_class)[i].size();
    size_t num_samples_inbag_class = round(num_samples * (*sample_fraction)[i]);

    shuffleAndSplitAppend(sampleIDs, oob_sampleIDs, num_samples_class, num_samples_inbag_class,
        (*sampleIDs_per_class)[i], random_number_generator);
  }
  num_samples_oob = oob_sampleIDs.size();

  if (keep_inbag) {
    // All observation are 0 or 1 times inbag
    inbag_counts.resize(num_samples, 1);
    for (size_t i = 0; i < oob_sampleIDs.size(); i++) {
      inbag_counts[oob_sampleIDs[i]] = 0;
    }
  }
}

} // namespace ranger
