/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#include <algorithm>
#include <iostream>
#include <iterator>

#include <ctime>

#include "utility.h"
#include "TreeRegression.h"
#include "Data.h"

namespace ranger {

TreeRegression::TreeRegression(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
    std::vector<double>& split_values) :
    Tree(child_nodeIDs, split_varIDs, split_values), counter(0), sums(0) {
}

void TreeRegression::allocateMemory() {
  // Init counters if not in memory efficient mode
  if (!memory_saving_splitting) {
    size_t max_num_splits = data->getMaxNumUniqueValues();

    // Use number of random splits for extratrees
    if (splitrule == EXTRATREES && num_random_splits > max_num_splits) {
      max_num_splits = num_random_splits;
    }

    counter.resize(max_num_splits);
    sums.resize(max_num_splits);
  }
}

double TreeRegression::estimate(size_t nodeID) {

  // Mean of responses of samples in node
  double sum_responses_in_node = sumNodeResponse(nodeID);

  size_t num_samples_in_node = end_pos[nodeID] - start_pos[nodeID];
  if (splitrule == POISSON && sum_responses_in_node == 0.) {
    // Poisson is not allowed to predict 0.
    // We use a weighted average of parent and child mean values,
    // see vignette "Introduction to Rpart" Chapter 8.2 and
    // https://ssrn.com/abstract=2870308 Chapter 6.1.3
    
    // Search for parent's nodeID: loop over all nodeIDs
    size_t parent_nodeID = 0;
    bool found = false;
    // Loop over left child nodes
    for(std::size_t i = 0; i < child_nodeIDs[0].size(); ++i) {
      // Break if parent node found
      if (child_nodeIDs[0][i] == nodeID) {
        parent_nodeID = i;
        found = true;
        break;
      }
    }
    if (!found) {
      // Loop over right child nodes
      for(std::size_t i = 0; i < child_nodeIDs[1].size(); ++i) {
        // Break if parent node found
        if (child_nodeIDs[1][i] == nodeID) {
          parent_nodeID = i;
          found = true;
          break;
        }
      }
    }
    
    double sum_responses_in_parent = sumNodeResponse(parent_nodeID);
    size_t num_samples_in_parent = end_pos[parent_nodeID] - start_pos[parent_nodeID];
    double mean_node = (sum_responses_in_node / (double) num_samples_in_node);
    double mean_parent = (sum_responses_in_parent / (double) num_samples_in_parent);
    double alpha = num_samples_in_node * mean_parent/(num_samples_in_node * mean_parent + poisson_tau);
    return alpha * mean_node + (1 - alpha) * mean_parent;
  } else {
    return (sum_responses_in_node / (double) num_samples_in_node);
  }
}

void TreeRegression::appendToFileInternal(std::ofstream& file) { // #nocov start
  // Empty on purpose
} // #nocov end

bool TreeRegression::splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  size_t num_samples_node = end_pos[nodeID] - start_pos[nodeID];
  
  // Save node statistics
  if (save_node_stats) {
    num_samples_nodes[nodeID] = num_samples_node;
    node_predictions[nodeID] = estimate(nodeID);
  }

  // Stop if maximum node size or depth reached
  if (num_samples_node <= (*min_node_size)[0] || (nodeID >= last_left_nodeID && max_depth > 0 && depth >= max_depth)) {
    split_values[nodeID] = estimate(nodeID);
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
    if (splitrule == POISSON && pure_value == 0.) {
      split_values[nodeID] = estimate(nodeID);
    } else {
      split_values[nodeID] = pure_value;
    }
    return true;
  }

  // Find best split, stop if no decrease of impurity
  bool stop;
  if (splitrule == MAXSTAT) {
    stop = findBestSplitMaxstat(nodeID, possible_split_varIDs);
  } else if (splitrule == EXTRATREES) {
    stop = findBestSplitExtraTrees(nodeID, possible_split_varIDs);
  } else if (splitrule == BETA) {
    stop = findBestSplitBeta(nodeID, possible_split_varIDs);
  } else if (splitrule == POISSON) {
    stop = findBestSplitPoisson(nodeID, possible_split_varIDs);
  } else {
    stop = findBestSplit(nodeID, possible_split_varIDs);
  }

  if (stop) {
    split_values[nodeID] = estimate(nodeID);
    return true;
  }

  return false;
}

void TreeRegression::createEmptyNodeInternal() {
  if (save_node_stats) {
    node_predictions.push_back(0);
  }
}

double TreeRegression::computePredictionAccuracyInternal(std::vector<double>* prediction_error_casewise) {

  size_t num_predictions = prediction_terminal_nodeIDs.size();
  double sum_of_squares = 0;
  for (size_t i = 0; i < num_predictions; ++i) {
    size_t terminal_nodeID = prediction_terminal_nodeIDs[i];
    double predicted_value = split_values[terminal_nodeID];
    double real_value = data->get_y(oob_sampleIDs[i], 0);
    if (predicted_value != real_value) {
      double diff = (predicted_value - real_value) * (predicted_value - real_value);
      if (prediction_error_casewise) {
        (*prediction_error_casewise)[i] = diff;
      }
      sum_of_squares += diff;
    }
  }
  return (1.0 - sum_of_squares / (double) num_predictions);
}

bool TreeRegression::findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  size_t num_samples_node = end_pos[nodeID] - start_pos[nodeID];
  double best_decrease = -1;
  size_t best_varID = 0;
  double best_value = 0;

  // Compute sum of responses in node
  double sum_node = sumNodeResponse(nodeID);

  // Stop early if no split posssible
  if (num_samples_node >= 2 * (*min_bucket)[0]) {

    // For all possible split variables
    for (auto& varID : possible_split_varIDs) {

      // Find best split value, if ordered consider all values as split values, else all 2-partitions
      if (data->isOrderedVariable(varID)) {

        // Use memory saving method if option set
        if (memory_saving_splitting) {
          findBestSplitValueSmallQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease);
        } else {
          // Use faster method for both cases
          double q = (double) num_samples_node / (double) data->getNumUniqueDataValues(varID);
          if (q < Q_THRESHOLD) {
            if (data->hasNA()) {
              findBestSplitValueNanSmallQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease);
            } else {
              findBestSplitValueSmallQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease);
            }
          } else {
            if (data->hasNA()) {
              findBestSplitValueNanLargeQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease); 
            } else {
              findBestSplitValueLargeQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease); 
            }
          }
        }
      } else {
        findBestSplitValueUnordered(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease);
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
  
  // Save split statistics
  if (save_node_stats) {
    split_stats[nodeID] = best_decrease;
  }

  // Compute decrease of impurity for this node and add to variable importance if needed
  if (importance_mode == IMP_GINI || importance_mode == IMP_GINI_CORRECTED) {
    addImpurityImportance(nodeID, best_varID, best_decrease);
  }

  // Regularization
  saveSplitVarID(best_varID);

  return false;
}

void TreeRegression::findBestSplitValueSmallQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
    double& best_value, size_t& best_varID, double& best_decrease) {

  // Create possible split values
  std::vector<double> possible_split_values;
  data->getAllValues(possible_split_values, sampleIDs, varID, start_pos[nodeID], end_pos[nodeID]);

  // Try next variable if all equal for this
  if (possible_split_values.size() < 2) {
    return;
  }

  const size_t num_splits = possible_split_values.size();
  if (memory_saving_splitting) {
    std::vector<double> sums_right(num_splits);
    std::vector<size_t> n_right(num_splits);
    findBestSplitValueSmallQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease,
        possible_split_values, sums_right, n_right);
  } else {
    std::fill_n(sums.begin(), num_splits, 0);
    std::fill_n(counter.begin(), num_splits, 0);
    findBestSplitValueSmallQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease,
        possible_split_values, sums, counter);
  }
}

void TreeRegression::findBestSplitValueSmallQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
    double& best_value, size_t& best_varID, double& best_decrease, std::vector<double> possible_split_values,
    std::vector<double>& sums, std::vector<size_t>& counter) {

  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    size_t idx = std::lower_bound(possible_split_values.begin(), possible_split_values.end(),
        data->get_x(sampleID, varID)) - possible_split_values.begin();

    sums[idx] += data->get_y(sampleID, 0);
    ++counter[idx];
  }

  size_t n_left = 0;
  double sum_left = 0;

  // Compute decrease of impurity for each split
  for (size_t i = 0; i < possible_split_values.size() - 1; ++i) {

    // Stop if nothing here
    if (counter[i] == 0) {
      continue;
    }

    n_left += counter[i];
    sum_left += sums[i];

    // Stop if right child empty
    size_t n_right = num_samples_node - n_left;
    if (n_right == 0) {
      break;
    }

    // Stop if minimal bucket size reached
    if (n_left < (*min_bucket)[0] || n_right < (*min_bucket)[0]) {
      continue;
    }

    double sum_right = sum_node - sum_left;
    double decrease = sum_left * sum_left / (double) n_left + sum_right * sum_right / (double) n_right;

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

void TreeRegression::findBestSplitValueLargeQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
    double& best_value, size_t& best_varID, double& best_decrease) {

  // Set counters to 0
  size_t num_unique = data->getNumUniqueDataValues(varID);
  std::fill_n(counter.begin(), num_unique, 0);
  std::fill_n(sums.begin(), num_unique, 0);

  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    size_t index = data->getIndex(sampleID, varID);

    sums[index] += data->get_y(sampleID, 0);
    ++counter[index];
  }

  size_t n_left = 0;
  double sum_left = 0;

  // Compute decrease of impurity for each split
  for (size_t i = 0; i < num_unique - 1; ++i) {

    // Stop if nothing here
    if (counter[i] == 0) {
      continue;
    }

    n_left += counter[i];
    sum_left += sums[i];

    // Stop if right child empty
    size_t n_right = num_samples_node - n_left;
    if (n_right == 0) {
      break;
    }

    // Stop if minimal bucket size reached
    if (n_left < (*min_bucket)[0] || n_right < (*min_bucket)[0]) {
      continue;
    }

    double sum_right = sum_node - sum_left;
    double decrease = sum_left * sum_left / (double) n_left + sum_right * sum_right / (double) n_right;

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

void TreeRegression::findBestSplitValueUnordered(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
    double& best_value, size_t& best_varID, double& best_decrease) {

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
    double sum_right = 0;
    size_t n_right = 0;

    // Sum in right child
    for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
      size_t sampleID = sampleIDs[pos];
      double response = data->get_y(sampleID, 0);
      double value = data->get_x(sampleID, varID);
      size_t factorID = floor(value) - 1;

      // If in right child, count
      // In right child, if bitwise splitID at position factorID is 1
      if ((splitID & (1ULL << factorID))) {
        ++n_right;
        sum_right += response;
      }
    }
    size_t n_left = num_samples_node - n_right;

    // Stop if minimal bucket size reached
    if (n_left < (*min_bucket)[0] || n_right < (*min_bucket)[0]) {
      continue;
    }

    // Sum of squares
    double sum_left = sum_node - sum_right;
    double decrease = sum_left * sum_left / (double) n_left + sum_right * sum_right / (double) n_right;

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

bool TreeRegression::findBestSplitMaxstat(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  size_t num_samples_node = end_pos[nodeID] - start_pos[nodeID];

  // Compute ranks
  std::vector<double> response;
  response.reserve(num_samples_node);
  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    response.push_back(data->get_y(sampleID, 0));
  }
  std::vector<double> ranks = rank(response);

  // Save split stats
  std::vector<double> pvalues;
  pvalues.reserve(possible_split_varIDs.size());
  std::vector<double> values;
  values.reserve(possible_split_varIDs.size());
  std::vector<double> candidate_varIDs;
  candidate_varIDs.reserve(possible_split_varIDs.size());
  std::vector<double> test_statistics;
  test_statistics.reserve(possible_split_varIDs.size());

  // Compute p-values
  for (auto& varID : possible_split_varIDs) {

    // Get all observations
    std::vector<double> x;
    x.reserve(num_samples_node);
    for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
      size_t sampleID = sampleIDs[pos];
      x.push_back(data->get_x(sampleID, varID));
    }

    // Order by x
    std::vector<size_t> indices = order(x, false);
    //std::vector<size_t> indices = orderInData(data, sampleIDs[nodeID], varID, false);

    // Compute maximally selected rank statistics
    double best_maxstat;
    double best_split_value;
    maxstat(ranks, x, indices, best_maxstat, best_split_value, minprop, 1 - minprop);
    //maxstatInData(scores, data, sampleIDs[nodeID], varID, indices, best_maxstat, best_split_value, minprop, 1 - minprop);

    if (best_maxstat > -1) {
      // Compute number of samples left of cutpoints
      std::vector<size_t> num_samples_left = numSamplesLeftOfCutpoint(x, indices);
      //std::vector<size_t> num_samples_left = numSamplesLeftOfCutpointInData(data, sampleIDs[nodeID], varID, indices);

      // Compute p-values
      double pvalue_lau92 = maxstatPValueLau92(best_maxstat, minprop, 1 - minprop);
      double pvalue_lau94 = maxstatPValueLau94(best_maxstat, minprop, 1 - minprop, num_samples_node, num_samples_left);

      // Use minimum of Lau92 and Lau94
      double pvalue = std::min(pvalue_lau92, pvalue_lau94);

      // Save split stats
      pvalues.push_back(pvalue);
      values.push_back(best_split_value);
      candidate_varIDs.push_back(varID);
      test_statistics.push_back(best_maxstat);
    }
  }

  double adjusted_best_pvalue = std::numeric_limits<double>::max();
  size_t best_varID = 0;
  double best_value = 0;
  double best_maxstat = 0;

  if (pvalues.size() > 0) {
    // Adjust p-values with Benjamini/Hochberg
    std::vector<double> adjusted_pvalues = adjustPvalues(pvalues);

    // Use smallest p-value
    double min_pvalue = std::numeric_limits<double>::max();
    for (size_t i = 0; i < pvalues.size(); ++i) {
      if (pvalues[i] < min_pvalue) {
        min_pvalue = pvalues[i];
        best_varID = candidate_varIDs[i];
        best_value = values[i];
        adjusted_best_pvalue = adjusted_pvalues[i];
        best_maxstat = test_statistics[i];
      }
    }
  }

  // Stop if no good split found (this is terminal node).
  if (adjusted_best_pvalue > alpha) {
    return true;
  } else {
    // If not terminal node save best values
    split_varIDs[nodeID] = best_varID;
    split_values[nodeID] = best_value;
    
    // Save split statistics
    if (save_node_stats) {
      split_stats[nodeID] = best_maxstat;
    }

    // Compute decrease of impurity for this node and add to variable importance if needed
    if (importance_mode == IMP_GINI || importance_mode == IMP_GINI_CORRECTED) {
      addImpurityImportance(nodeID, best_varID, best_maxstat);
    }

    return false;
  }
}

bool TreeRegression::findBestSplitExtraTrees(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  size_t num_samples_node = end_pos[nodeID] - start_pos[nodeID];
  double best_decrease = -1;
  size_t best_varID = 0;
  double best_value = 0;

  // Compute sum of responses in node
  double sum_node = sumNodeResponse(nodeID);

  // Stop early if no split posssible
  if (num_samples_node >= 2 * (*min_bucket)[0]) {
  
    // For all possible split variables
    for (auto& varID : possible_split_varIDs) {

      // Find best split value, if ordered consider all values as split values, else all 2-partitions
      if (data->isOrderedVariable(varID)) {
        findBestSplitValueExtraTrees(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease);
      } else {
        findBestSplitValueExtraTreesUnordered(nodeID, varID, sum_node, num_samples_node, best_value, best_varID,
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
  
  // Save split statistics
  if (save_node_stats) {
    split_stats[nodeID] = best_decrease;
  }

  // Compute decrease of impurity for this node and add to variable importance if needed
  if (importance_mode == IMP_GINI || importance_mode == IMP_GINI_CORRECTED) {
    addImpurityImportance(nodeID, best_varID, best_decrease);
  }

  // Regularization
  saveSplitVarID(best_varID);

  return false;
}

void TreeRegression::findBestSplitValueExtraTrees(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
    double& best_value, size_t& best_varID, double& best_decrease) {

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
    std::vector<double> sums_right(num_splits);
    std::vector<size_t> n_right(num_splits);
    findBestSplitValueExtraTrees(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease,
        possible_split_values, sums_right, n_right);
  } else {
    std::fill_n(sums.begin(), num_splits, 0);
    std::fill_n(counter.begin(), num_splits, 0);
    findBestSplitValueExtraTrees(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease,
        possible_split_values, sums, counter);
  }
}

void TreeRegression::findBestSplitValueExtraTrees(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
    double& best_value, size_t& best_varID, double& best_decrease, std::vector<double> possible_split_values,
    std::vector<double>& sums_right, std::vector<size_t>& n_right) {
  const size_t num_splits = possible_split_values.size();

  // Sum in right child and possbile split
  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    double value = data->get_x(sampleID, varID);
    double response = data->get_y(sampleID, 0);

    // Count samples until split_value reached
    for (size_t i = 0; i < num_splits; ++i) {
      if (value > possible_split_values[i]) {
        ++n_right[i];
        sums_right[i] += response;
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
    if (n_left < (*min_bucket)[0] || n_right[i] < (*min_bucket)[0]) {
      continue;
    }

    double sum_right = sums_right[i];
    double sum_left = sum_node - sum_right;
    double decrease = sum_left * sum_left / (double) n_left + sum_right * sum_right / (double) n_right[i];

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

void TreeRegression::findBestSplitValueExtraTreesUnordered(size_t nodeID, size_t varID, double sum_node,
    size_t num_samples_node, double& best_value, size_t& best_varID, double& best_decrease) {

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
    double sum_right = 0;
    size_t n_right = 0;

    // Sum in right child
    for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
      size_t sampleID = sampleIDs[pos];
      double response = data->get_y(sampleID, 0);
      double value = data->get_x(sampleID, varID);
      size_t factorID = floor(value) - 1;

      // If in right child, count
      // In right child, if bitwise splitID at position factorID is 1
      if ((splitID & (1ULL << factorID))) {
        ++n_right;
        sum_right += response;
      }
    }
    size_t n_left = num_samples_node - n_right;

    // Stop if minimal bucket size reached
    if (n_left < (*min_bucket)[0] || n_right < (*min_bucket)[0]) {
      continue;
    }

    // Sum of squares
    double sum_left = sum_node - sum_right;
    double decrease = sum_left * sum_left / (double) n_left + sum_right * sum_right / (double) n_right;

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

bool TreeRegression::findBestSplitBeta(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  size_t num_samples_node = end_pos[nodeID] - start_pos[nodeID];
  double best_decrease = -std::numeric_limits<double>::infinity();
  size_t best_varID = 0;
  double best_value = 0;

  // Compute sum of responses in node
  double sum_node = sumNodeResponse(nodeID);

  // Stop early if no split posssible
  if (num_samples_node >= 2 * (*min_bucket)[0]) {

    // For all possible split variables find best split value
    for (auto& varID : possible_split_varIDs) {
      findBestSplitValueBeta(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease);
    }
  }

  // Stop if no good split found
  if (std::isinf(-best_decrease)) {
    return true;
  }

  // Save best values
  split_varIDs[nodeID] = best_varID;
  split_values[nodeID] = best_value;
  
  // Save split statistics
  if (save_node_stats) {
    split_stats[nodeID] = best_decrease;
  }

  // Compute decrease of impurity for this node and add to variable importance if needed
  if (importance_mode == IMP_GINI || importance_mode == IMP_GINI_CORRECTED) {
    addImpurityImportance(nodeID, best_varID, best_decrease);
  }

  // Regularization
  saveSplitVarID(best_varID);

  return false;
}

void TreeRegression::findBestSplitValueBeta(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
    double& best_value, size_t& best_varID, double& best_decrease) {

  // Create possible split values
  std::vector<double> possible_split_values;
  data->getAllValues(possible_split_values, sampleIDs, varID, start_pos[nodeID], end_pos[nodeID]);

  // Try next variable if all equal for this
  if (possible_split_values.size() < 2) {
    return;
  }

  // -1 because no split possible at largest value
  size_t num_splits = possible_split_values.size() - 1;
  if (memory_saving_splitting) {
    std::vector<double> sums_right(num_splits);
    std::vector<size_t> n_right(num_splits);
    findBestSplitValueBeta(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease,
        possible_split_values, sums_right, n_right);
  } else {
    std::fill_n(sums.begin(), num_splits, 0);
    std::fill_n(counter.begin(), num_splits, 0);
    findBestSplitValueBeta(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease,
        possible_split_values, sums, counter);
  }
}

void TreeRegression::findBestSplitValueBeta(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
    double& best_value, size_t& best_varID, double& best_decrease, std::vector<double> possible_split_values,
    std::vector<double>& sums_right, std::vector<size_t>& n_right) {
  // -1 because no split possible at largest value
  const size_t num_splits = possible_split_values.size() - 1;

  // Sum in right child and possbile split
  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    double value = data->get_x(sampleID, varID);
    double response = data->get_y(sampleID, 0);

    // Count samples until split_value reached
    for (size_t i = 0; i < num_splits; ++i) {
      if (value > possible_split_values[i]) {
        ++n_right[i];
        sums_right[i] += response;
      } else {
        break;
      }
    }
  }

  // Compute LogLik of beta distribution for each possible split
  for (size_t i = 0; i < num_splits; ++i) {

    // Stop if one child too small
    size_t n_left = num_samples_node - n_right[i];
    if (n_left < 2 || n_right[i] < 2) {
      continue;
    }

    // Stop if minimal bucket size reached
    if (n_left < (*min_bucket)[0] || n_right[i] < (*min_bucket)[0]) {
      continue;
    }

    // Compute mean
    double sum_right = sums_right[i];
    double mean_right = sum_right / (double) n_right[i];
    double sum_left = sum_node - sum_right;
    double mean_left = sum_left / (double) n_left;

    // Compute variance
    double var_right = 0;
    double var_left = 0;
    for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
      size_t sampleID = sampleIDs[pos];
      double value = data->get_x(sampleID, varID);
      double response = data->get_y(sampleID, 0);

      if (value > possible_split_values[i]) {
        var_right += (response - mean_right) * (response - mean_right);
      } else {
        var_left += (response - mean_left) * (response - mean_left);
      }
    }
    var_right /= (double) n_right[i] - 1;
    var_left /= (double) n_left - 1;

    // Stop if zero variance
    if (var_right < std::numeric_limits<double>::epsilon() || var_left < std::numeric_limits<double>::epsilon()) {
      continue;
    }

    // Compute phi for beta distribution
    double phi_right = mean_right * (1 - mean_right) / var_right - 1;
    double phi_left = mean_left * (1 - mean_left) / var_left - 1;

    // Compute LogLik of beta distribution
    double beta_loglik_right = 0;
    double beta_loglik_left = 0;
    for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
      size_t sampleID = sampleIDs[pos];
      double value = data->get_x(sampleID, varID);
      double response = data->get_y(sampleID, 0);

      if (value > possible_split_values[i]) {
        beta_loglik_right += betaLogLik(response, mean_right, phi_right);
      } else {
        beta_loglik_left += betaLogLik(response, mean_left, phi_left);
      }
    }

    // Split statistic is sum of both log-likelihoods
    double decrease = beta_loglik_right + beta_loglik_left;

    // Stop if no result
    if (std::isnan(decrease)) {
      continue;
    }

    // Regularization (negative values)
    regularizeNegative(decrease, varID);

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
}

bool TreeRegression::findBestSplitPoisson(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {
  
  size_t num_samples_node = end_pos[nodeID] - start_pos[nodeID];
  double best_decrease = -std::numeric_limits<double>::infinity();
  size_t best_varID = 0;
  double best_value = 0;
  
  // Compute sum of responses in node
  double sum_node = sumNodeResponse(nodeID);
  
  // Stop early if no split posssible
  if (num_samples_node >= 2 * (*min_bucket)[0]) {
    
    // For all possible split variables find best split value
    for (auto& varID : possible_split_varIDs) {
      findBestSplitValuePoissonSmallQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID,
                                      best_decrease);
    }
  }
  
  // Stop if no good split found
  if (std::isinf(-best_decrease)) {
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

void TreeRegression::findBestSplitValuePoissonSmallQ(size_t nodeID, size_t varID, double sum_node,
                                                     size_t num_samples_node, double& best_value, size_t& best_varID, double& best_decrease) {
  
  // Create possible split values
  std::vector<double> possible_split_values;
  data->getAllValues(possible_split_values, sampleIDs, varID, start_pos[nodeID], end_pos[nodeID]);
  
  // Try next variable if all equal for this
  if (possible_split_values.size() < 2) {
    return;
  }
  
  // -1 because no split possible at largest value
  const size_t num_splits = possible_split_values.size() - 1;
  if (memory_saving_splitting) {
    std::vector<double> sums_right(num_splits);
    std::vector<size_t> n_right(num_splits);
    findBestSplitValuePoissonSmallQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease,
                                    possible_split_values, sums_right, n_right);
  } else {
    std::fill_n(sums.begin(), num_splits, 0);
    std::fill_n(counter.begin(), num_splits, 0);
    findBestSplitValuePoissonSmallQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease,
                                    possible_split_values, sums, counter);
  }
}

void TreeRegression::findBestSplitValuePoissonSmallQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
                                                     double& best_value, size_t& best_varID, double& best_decrease, std::vector<double> possible_split_values,
                                                     std::vector<double>& sums, std::vector<size_t>& counter) {
  
  // Sum and sample count for possbile splits
  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    size_t idx = std::lower_bound(possible_split_values.begin(), possible_split_values.end(),
                                  data->get_x(sampleID, varID)) - possible_split_values.begin();
    
    sums[idx] += data->get_y(sampleID, 0);
    ++counter[idx];
  }
  
  size_t n_left = 0;
  double sum_left = 0;
  
  // Compute decrease in Poisson deviance for each possible split
  for (size_t i = 0; i < possible_split_values.size() - 1; ++i) {
    
    // Stop if nothing here
    if (counter[i] == 0) {
      continue;
    }
    
    n_left += counter[i];
    sum_left += sums[i];
    
    // Stop if right child empty
    size_t n_right = num_samples_node - n_left;
    if (n_right == 0) {
      break;
    }
    
    // Stop if minimal bucket size reached
    if (n_left < (*min_bucket)[0] || n_right < (*min_bucket)[0]) {
      continue;
    }
    
    // Compute mean
    double sum_right = sum_node - sum_left;
    double mean_right = sum_right / (double) n_right;
    double mean_left = sum_left / (double) n_left;
    
    // Poisson deviance = 2 * (y_true * log(y_true/y_pred) + y_pred - y_true)
    // decrease = - 1/2 * (sum_left(poisson_deviance) + sum_right(poisson_deviance))
    //          = + sum_left(y) * log(mean_left) + sum_right(y) * log(mean_right) + const + 0
    // The smaller the deviance, the better => the larger the decrease, the better.
    double decrease = xlogy(sum_left, mean_left) + xlogy(sum_right, mean_right);
    
    // Stop if no result
    if (std::isnan(decrease)) {
      continue;
    }
    
    // Regularization
    if (decrease > 0) {
      regularize(decrease, varID);
    } else {
      regularizeNegative(decrease, varID);
    }
    
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
}

void TreeRegression::findBestSplitValueNanSmallQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
                                                 double& best_value, size_t& best_varID, double& best_decrease) {
  
  // Create possible split values
  std::vector<double> possible_split_values;
  data->getAllValues(possible_split_values, sampleIDs, varID, start_pos[nodeID], end_pos[nodeID]);
  
  // Try next variable if all equal for this
  if (possible_split_values.size() < 2) {
    return;
  }
  
  const size_t num_splits = possible_split_values.size();
  if (memory_saving_splitting) {
    std::vector<double> sums_right(num_splits);
    std::vector<size_t> n_right(num_splits);
    findBestSplitValueSmallQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease,
                             possible_split_values, sums_right, n_right);
  } else {
    std::fill_n(sums.begin(), num_splits, 0);
    std::fill_n(counter.begin(), num_splits, 0);
    findBestSplitValueSmallQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease,
                             possible_split_values, sums, counter);
  }
}

void TreeRegression::findBestSplitValueNanSmallQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
                                                 double& best_value, size_t& best_varID, double& best_decrease, std::vector<double> possible_split_values,
                                                 std::vector<double>& sums, std::vector<size_t>& counter) {
  
  // Counters without NaNs
  double sum_nan = 0;
  size_t num_samples_node_nan = 0;
  
  size_t last_index = possible_split_values.size() - 1;
  if (std::isnan(possible_split_values[last_index])) {
    for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
      size_t sampleID = sampleIDs[pos];
      
      if (std::isnan(data->get_x(sampleID, varID))) {
        sum_nan += data->get_y(sampleID, 0);
        ++num_samples_node_nan;
      } else {
        size_t idx = std::lower_bound(possible_split_values.begin(), possible_split_values.end(),
                                      data->get_x(sampleID, varID)) - possible_split_values.begin();
        
        sums[idx] += data->get_y(sampleID, 0);
        ++counter[idx];
      }
    }
  } else {
    for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
      size_t sampleID = sampleIDs[pos];
      size_t idx = std::lower_bound(possible_split_values.begin(), possible_split_values.end(),
                                    data->get_x(sampleID, varID)) - possible_split_values.begin();
      
      sums[idx] += data->get_y(sampleID, 0);
      ++counter[idx];
    }
  }
  
  size_t n_left = 0;
  double sum_left = 0;

  // If MIA approach, use one additional split
  size_t num_possible_splits = possible_split_values.size() - 1;
  if (mia) {
    num_possible_splits = possible_split_values.size();
  } 
  
  // Compute decrease of impurity for each split
  for (size_t i = 0; i < num_possible_splits; ++i) {
    
    // Stop if nothing here
    if (counter[i] == 0 && (!mia || i < num_possible_splits - 1)) {
      continue;
    }
    
    n_left += counter[i];
    sum_left += sums[i];
    
    // Stop if right child empty
    size_t n_right = num_samples_node - num_samples_node_nan - n_left;
    if (n_right == 0 && (!mia || num_samples_node_nan == 0)) {
      break;
    }
    
    // Stop if minimal bucket size reached
    if ((n_left < (*min_bucket)[0] || n_right < (*min_bucket)[0]) && 
    (!mia || i < num_possible_splits - 1)) {
      continue;
    }

    // Stop if only NAs
    if ((n_right + n_left) == 0) {
      continue;
    }
    
    double sum_right = sum_node - sum_left - sum_nan;
    double decrease = sum_left * sum_left / (double) n_left + sum_right * sum_right / (double) n_right;
    
    double decrease_nanleft = (sum_left + sum_nan) * (sum_left + sum_nan)  / (double) (n_left + num_samples_node_nan) + sum_right * sum_right / (double) n_right;
    double decrease_nanright = sum_left * sum_left / (double) n_left + (sum_right + sum_nan)  * (sum_right + sum_nan)  / (double) (n_right + num_samples_node_nan);

    // If right child empty, put NAs right
    if (n_right == 0) {
      decrease = decrease_nanright;
      decrease_nanleft = 0;
    }
    
    // Regularization
    regularize(decrease, varID);
    
    // If better than before, use this
    if (decrease > best_decrease) {
      // Use mid-point split
      best_value = (possible_split_values[i] + possible_split_values[i + 1]) / 2;
      best_varID = varID;
      best_decrease = decrease;

      if (n_right == 0) {
        // Split value infinity -> all non-missing values go left
        best_value = std::numeric_limits<double>::infinity(); 
      }
      
      if (decrease_nanright > decrease_nanleft) {
        nan_go_right = true;
      } else {
        nan_go_right = false;
      }
      
      // Use smaller value if average is numerically the same as the larger value
      if (best_value == possible_split_values[i + 1]) {
        best_value = possible_split_values[i];
      }
    }
  }
}

void TreeRegression::findBestSplitValueNanLargeQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
                                                 double& best_value, size_t& best_varID, double& best_decrease) {
  
  // Set counters to 0
  size_t num_unique = data->getNumUniqueDataValues(varID);
  std::fill_n(counter.begin(), num_unique, 0);
  std::fill_n(sums.begin(), num_unique, 0);
  
  // Counters without NaNs
  double sum_nan = 0;
  size_t num_samples_node_nan = 0;
  
  size_t last_index = data->getNumUniqueDataValues(varID) - 1;
  if (std::isnan(data->getUniqueDataValue(varID, last_index))) {
    for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
      size_t sampleID = sampleIDs[pos];
      
      if (std::isnan(data->get_x(sampleID, varID))) {
        sum_nan += data->get_y(sampleID, 0);
        ++num_samples_node_nan;
      } else {
        size_t index = data->getIndex(sampleID, varID);
        sums[index] += data->get_y(sampleID, 0);
        ++counter[index];
      }
    }
  } else {
    for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
      size_t sampleID = sampleIDs[pos];
      size_t index = data->getIndex(sampleID, varID);
      
      sums[index] += data->get_y(sampleID, 0);
      ++counter[index];
    }
  }
  
  
  size_t n_left = 0;
  double sum_left = 0;

  // If MIA approach, use one additional split
  size_t num_possible_splits = num_unique - 1;
  if (mia) {
    num_possible_splits = num_unique;
  } 
  
  // Compute decrease of impurity for each split
  for (size_t i = 0; i < num_possible_splits; ++i) {
    
    // Stop if nothing here
    if (counter[i] == 0 && (!mia || i < num_possible_splits - 1)) {
      continue;
    }
    
    n_left += counter[i];
    sum_left += sums[i];
    
    // Stop if right child empty
    size_t n_right = num_samples_node - num_samples_node_nan - n_left;
    if (n_right == 0 && (!mia || num_samples_node_nan == 0)) {
      break;
    }
    
    // Stop if minimal bucket size reached
    if ((n_left < (*min_bucket)[0] || n_right < (*min_bucket)[0]) && 
    (!mia || i < num_possible_splits - 1)) {
      continue;
    }

    // Stop if only NAs
    if ((n_right + n_left) == 0) {
      continue;
    }
    
    double sum_right = sum_node - sum_left;
    double decrease = sum_left * sum_left / (double) n_left + sum_right * sum_right / (double) n_right;
    
    double decrease_nanleft = (sum_left + sum_nan) * (sum_left + sum_nan)  / (double) (n_left + num_samples_node_nan) + sum_right * sum_right / (double) n_right;
    double decrease_nanright = sum_left * sum_left / (double) n_left + (sum_right + sum_nan)  * (sum_right + sum_nan)  / (double) (n_right + num_samples_node_nan);
    
    // If right child empty, put NAs right
    if (n_right == 0) {
      decrease = decrease_nanright;
      decrease_nanleft = 0;
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

      if (n_right == 0) {
        // Split value infinity -> all non-missing values go left
        best_value = std::numeric_limits<double>::infinity(); 
      }
      
      if (decrease_nanright > decrease_nanleft) {
        nan_go_right = true;
      } else {
        nan_go_right = false;
      }
      
      // Use smaller value if average is numerically the same as the larger value
      if (best_value == data->getUniqueDataValue(varID, j)) {
        best_value = data->getUniqueDataValue(varID, i);
      }
    }
  }
}

void TreeRegression::addImpurityImportance(size_t nodeID, size_t varID, double decrease) {

  size_t num_samples_node = end_pos[nodeID] - start_pos[nodeID];

  double best_decrease = decrease;
  if (splitrule != MAXSTAT) {
    double sum_node = sumNodeResponse(nodeID);
    double impurity_node = (sum_node * sum_node / (double) num_samples_node);

    // Account for the regularization
    regularize(impurity_node, varID);

    best_decrease = decrease - impurity_node;
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

} // namespace ranger
