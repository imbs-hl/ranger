/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>

#include "utility.h"
#include "TreeSurvival.h"
#include "Data.h"

namespace ranger {

TreeSurvival::TreeSurvival(std::vector<double>* unique_timepoints, size_t status_varID,
    std::vector<size_t>* response_timepointIDs) :
    status_varID(status_varID), unique_timepoints(unique_timepoints), response_timepointIDs(response_timepointIDs), num_deaths(
        0), num_samples_at_risk(0) {
  this->num_timepoints = unique_timepoints->size();
}

TreeSurvival::TreeSurvival(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
    std::vector<double>& split_values, std::vector<std::vector<double>> chf, std::vector<double>* unique_timepoints,
    std::vector<size_t>* response_timepointIDs) :
    Tree(child_nodeIDs, split_varIDs, split_values), status_varID(0), unique_timepoints(unique_timepoints), response_timepointIDs(
        response_timepointIDs), chf(chf), num_deaths(0), num_samples_at_risk(0) {
  this->num_timepoints = unique_timepoints->size();
}

void TreeSurvival::allocateMemory() {
  // Number of deaths and samples at risk for each timepoint
  num_deaths.resize(num_timepoints);
  num_samples_at_risk.resize(num_timepoints);
}

void TreeSurvival::appendToFileInternal(std::ofstream& file) {  // #nocov start

  // Convert to vector without empty elements and save
  std::vector<size_t> terminal_nodes;
  std::vector<std::vector<double>> chf_vector;
  for (size_t i = 0; i < chf.size(); ++i) {
    if (!chf[i].empty()) {
      terminal_nodes.push_back(i);
      chf_vector.push_back(chf[i]);
    }
  }
  saveVector1D(terminal_nodes, file);
  saveVector2D(chf_vector, file);
} // #nocov end

void TreeSurvival::createEmptyNodeInternal() {
  chf.push_back(std::vector<double>());
}

void TreeSurvival::computeSurvival(size_t nodeID) {
  std::vector<double> chf_temp;
  chf_temp.reserve(num_timepoints);
  double chf_value = 0;
  for (size_t i = 0; i < num_timepoints; ++i) {
    if (num_samples_at_risk[i] != 0) {
      chf_value += (double) num_deaths[i] / (double) num_samples_at_risk[i];
    }
    chf_temp.push_back(chf_value);
  }
  chf[nodeID] = chf_temp;
}

double TreeSurvival::computePredictionAccuracyInternal() {

  // Compute summed chf for samples
  std::vector<double> sum_chf;
  for (size_t i = 0; i < prediction_terminal_nodeIDs.size(); ++i) {
    size_t terminal_nodeID = prediction_terminal_nodeIDs[i];
    sum_chf.push_back(std::accumulate(chf[terminal_nodeID].begin(), chf[terminal_nodeID].end(), 0.0));
  }

  // Return concordance index
  return computeConcordanceIndex(*data, sum_chf, dependent_varID, status_varID, oob_sampleIDs);
}

bool TreeSurvival::splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  if (splitrule == MAXSTAT) {
    return findBestSplitMaxstat(nodeID, possible_split_varIDs);
  } else if (splitrule == EXTRATREES) {
    return findBestSplitExtraTrees(nodeID, possible_split_varIDs);
  } else {
    return findBestSplit(nodeID, possible_split_varIDs);
  }
}

bool TreeSurvival::findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  double best_decrease = -1;
  size_t num_samples_node = sampleIDs[nodeID].size();
  size_t best_varID = 0;
  double best_value = 0;

  computeDeathCounts(nodeID);

  // Stop early if no split posssible
  if (num_samples_node >= 2 * min_node_size) {

    // For all possible split variables
    for (auto& varID : possible_split_varIDs) {

      // Find best split value, if ordered consider all values as split values, else all 2-partitions
      if (data->isOrderedVariable(varID)) {
        if (splitrule == LOGRANK) {
          findBestSplitValueLogRank(nodeID, varID, best_value, best_varID, best_decrease);
        } else if (splitrule == AUC || splitrule == AUC_IGNORE_TIES) {
          findBestSplitValueAUC(nodeID, varID, best_value, best_varID, best_decrease);
        }
      } else {
        findBestSplitValueLogRankUnordered(nodeID, varID, best_value, best_varID, best_decrease);
      }

    }
  }

  // Stop and save CHF if no good split found (this is terminal node).
  if (best_decrease < 0) {
    computeSurvival(nodeID);
    return true;
  } else {
    // If not terminal node save best values
    split_varIDs[nodeID] = best_varID;
    split_values[nodeID] = best_value;

    // Compute decrease of impurity for this node and add to variable importance if needed
    if (importance_mode == IMP_GINI || importance_mode == IMP_GINI_CORRECTED) {
      addImpurityImportance(nodeID, best_varID, best_decrease);
    }

    return false;
  }
}

bool TreeSurvival::findBestSplitMaxstat(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  size_t num_samples_node = sampleIDs[nodeID].size();

  // Check node size, stop if maximum reached
  if (num_samples_node <= min_node_size) {
    computeDeathCounts(nodeID);
    computeSurvival(nodeID);
    return true;
  }

  // Compute scores
  std::vector<double> time;
  time.reserve(num_samples_node);
  std::vector<double> status;
  status.reserve(num_samples_node);
  for (auto& sampleID : sampleIDs[nodeID]) {
    time.push_back(data->get(sampleID, dependent_varID));
    status.push_back(data->get(sampleID, status_varID));
  }
  std::vector<double> scores = logrankScores(time, status);
  //std::vector<double> scores = logrankScoresData(data, dependent_varID, status_varID, sampleIDs[nodeID]);

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
    for (auto& sampleID : sampleIDs[nodeID]) {
      x.push_back(data->get(sampleID, varID));
    }

    // Order by x
    std::vector<size_t> indices = order(x, false);
    //std::vector<size_t> indices = orderInData(data, sampleIDs[nodeID], varID, false);

    // Compute maximally selected rank statistics
    double best_maxstat;
    double best_split_value;
    maxstat(scores, x, indices, best_maxstat, best_split_value, minprop, 1 - minprop);
    //maxstatInData(scores, data, sampleIDs[nodeID], varID, indices, best_maxstat, best_split_value, minprop, 1 - minprop);

    if (best_maxstat > -1) {
      // Compute number of samples left of cutpoints
      std::vector<size_t> num_samples_left = numSamplesLeftOfCutpoint(x, indices);
      //std::vector<size_t> num_samples_left = numSamplesLeftOfCutpointInData(data, sampleIDs[nodeID], varID, indices);

      // Remove largest cutpoint (all observations left)
      num_samples_left.pop_back();

      // Use unadjusted p-value if only 1 split point
      double pvalue;
      if (num_samples_left.size() == 1) {
        pvalue = maxstatPValueUnadjusted(best_maxstat);
      } else {
        // Compute p-values
        double pvalue_lau92 = maxstatPValueLau92(best_maxstat, minprop, 1 - minprop);
        double pvalue_lau94 = maxstatPValueLau94(best_maxstat, minprop, 1 - minprop, num_samples_node,
            num_samples_left);

        // Use minimum of Lau92 and Lau94
        pvalue = std::min(pvalue_lau92, pvalue_lau94);
      }

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

  // Stop and save CHF if no good split found (this is terminal node).
  if (adjusted_best_pvalue > alpha) {
    computeDeathCounts(nodeID);
    computeSurvival(nodeID);
    return true;
  } else {
    // If not terminal node save best values
    split_varIDs[nodeID] = best_varID;
    split_values[nodeID] = best_value;

    // Compute decrease of impurity for this node and add to variable importance if needed
    if (importance_mode == IMP_GINI || importance_mode == IMP_GINI_CORRECTED) {
      addImpurityImportance(nodeID, best_varID, best_maxstat);
    }

    return false;
  }
}

void TreeSurvival::computeDeathCounts(size_t nodeID) {

  // Initialize
  for (size_t i = 0; i < num_timepoints; ++i) {
    num_deaths[i] = 0;
    num_samples_at_risk[i] = 0;
  }

  for (auto& sampleID : sampleIDs[nodeID]) {
    double survival_time = data->get(sampleID, dependent_varID);

    size_t t = 0;
    while (t < num_timepoints && (*unique_timepoints)[t] < survival_time) {
      ++num_samples_at_risk[t];
      ++t;
    }

    // Now t is the survival time, add to at risk and to death if death
    if (t < num_timepoints) {
      ++num_samples_at_risk[t];
      if (data->get(sampleID, status_varID) == 1) {
        ++num_deaths[t];
      }
    }
  }
}

void TreeSurvival::computeChildDeathCounts(size_t nodeID, size_t varID, std::vector<double>& possible_split_values,
    std::vector<size_t>& num_samples_right_child, std::vector<size_t>& delta_samples_at_risk_right_child,
    std::vector<size_t>& num_deaths_right_child, size_t num_splits) {

  // Count deaths in right child per timepoint and possbile split
  for (auto& sampleID : sampleIDs[nodeID]) {
    double value = data->get(sampleID, varID);
    size_t survival_timeID = (*response_timepointIDs)[sampleID];

    // Count deaths until split_value reached
    for (size_t i = 0; i < num_splits; ++i) {

      if (value > possible_split_values[i]) {
        ++num_samples_right_child[i];
        ++delta_samples_at_risk_right_child[i * num_timepoints + survival_timeID];
        if (data->get(sampleID, status_varID) == 1) {
          ++num_deaths_right_child[i * num_timepoints + survival_timeID];
        }
      } else {
        break;
      }
    }
  }
}

void TreeSurvival::findBestSplitValueLogRank(size_t nodeID, size_t varID, double& best_value, size_t& best_varID,
    double& best_logrank) {

  // Create possible split values
  std::vector<double> possible_split_values;
  data->getAllValues(possible_split_values, sampleIDs[nodeID], varID);

  // Try next variable if all equal for this
  if (possible_split_values.size() < 2) {
    return;
  }

  // -1 because no split possible at largest value
  size_t num_splits = possible_split_values.size() - 1;

  // Initialize
  std::vector<size_t> num_deaths_right_child(num_splits * num_timepoints);
  std::vector<size_t> delta_samples_at_risk_right_child(num_splits * num_timepoints);
  std::vector<size_t> num_samples_right_child(num_splits);

  computeChildDeathCounts(nodeID, varID, possible_split_values, num_samples_right_child,
      delta_samples_at_risk_right_child, num_deaths_right_child, num_splits);

  // Compute logrank test for all splits and use best
  for (size_t i = 0; i < num_splits; ++i) {
    double numerator = 0;
    double denominator_squared = 0;

    // Stop if minimal node size reached
    size_t num_samples_left_child = sampleIDs[nodeID].size() - num_samples_right_child[i];
    if (num_samples_right_child[i] < min_node_size || num_samples_left_child < min_node_size) {
      continue;
    }

    // Compute logrank test statistic for this split
    size_t num_samples_at_risk_right_child = num_samples_right_child[i];
    for (size_t t = 0; t < num_timepoints; ++t) {
      if (num_samples_at_risk[t] < 2 || num_samples_at_risk_right_child < 1) {
        break;
      }

      if (num_deaths[t] > 0) {
        // Numerator and demoninator for log-rank test, notation from Ishwaran et al.
        double di = (double) num_deaths[t];
        double di1 = (double) num_deaths_right_child[i * num_timepoints + t];
        double Yi = (double) num_samples_at_risk[t];
        double Yi1 = (double) num_samples_at_risk_right_child;
        numerator += di1 - Yi1 * (di / Yi);
        denominator_squared += (Yi1 / Yi) * (1.0 - Yi1 / Yi) * ((Yi - di) / (Yi - 1)) * di;
      }

      // Reduce number of samples at risk for next timepoint
      num_samples_at_risk_right_child -= delta_samples_at_risk_right_child[i * num_timepoints + t];

    }
    double logrank = -1;
    if (denominator_squared != 0) {
      logrank = fabs(numerator / sqrt(denominator_squared));
    }

    if (logrank > best_logrank) {
      best_value = (possible_split_values[i] + possible_split_values[i + 1]) / 2;
      best_varID = varID;
      best_logrank = logrank;

      // Use smaller value if average is numerically the same as the larger value
      if (best_value == possible_split_values[i + 1]) {
        best_value = possible_split_values[i];
      }
    }
  }
}

void TreeSurvival::findBestSplitValueLogRankUnordered(size_t nodeID, size_t varID, double& best_value,
    size_t& best_varID, double& best_logrank) {

  // Create possible split values
  std::vector<double> factor_levels;
  data->getAllValues(factor_levels, sampleIDs[nodeID], varID);

  // Try next variable if all equal for this
  if (factor_levels.size() < 2) {
    return;
  }

  // Number of possible splits is 2^num_levels
  size_t num_splits = (1 << factor_levels.size());

  // Compute logrank test statistic for each possible split
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
    std::vector<size_t> num_deaths_right_child(num_timepoints);
    std::vector<size_t> delta_samples_at_risk_right_child(num_timepoints);
    size_t num_samples_right_child = 0;
    double numerator = 0;
    double denominator_squared = 0;

    // Count deaths in right child per timepoint
    for (auto& sampleID : sampleIDs[nodeID]) {
      size_t survival_timeID = (*response_timepointIDs)[sampleID];
      double value = data->get(sampleID, varID);
      size_t factorID = floor(value) - 1;

      // If in right child, count
      // In right child, if bitwise splitID at position factorID is 1
      if ((splitID & (1 << factorID))) {
        ++num_samples_right_child;
        ++delta_samples_at_risk_right_child[survival_timeID];
        if (data->get(sampleID, status_varID) == 1) {
          ++num_deaths_right_child[survival_timeID];
        }
      }

    }

    // Stop if minimal node size reached
    size_t num_samples_left_child = sampleIDs[nodeID].size() - num_samples_right_child;
    if (num_samples_right_child < min_node_size || num_samples_left_child < min_node_size) {
      continue;
    }

    // Compute logrank test statistic for this split
    size_t num_samples_at_risk_right_child = num_samples_right_child;
    for (size_t t = 0; t < num_timepoints; ++t) {
      if (num_samples_at_risk[t] < 2 || num_samples_at_risk_right_child < 1) {
        break;
      }

      if (num_deaths[t] > 0) {
        // Numerator and demoninator for log-rank test, notation from Ishwaran et al.
        double di = (double) num_deaths[t];
        double di1 = (double) num_deaths_right_child[t];
        double Yi = (double) num_samples_at_risk[t];
        double Yi1 = (double) num_samples_at_risk_right_child;
        numerator += di1 - Yi1 * (di / Yi);
        denominator_squared += (Yi1 / Yi) * (1.0 - Yi1 / Yi) * ((Yi - di) / (Yi - 1)) * di;
      }

      // Reduce number of samples at risk for next timepoint
      num_samples_at_risk_right_child -= delta_samples_at_risk_right_child[t];
    }
    double logrank = -1;
    if (denominator_squared != 0) {
      logrank = fabs(numerator / sqrt(denominator_squared));
    }

    if (logrank > best_logrank) {
      best_value = splitID;
      best_varID = varID;
      best_logrank = logrank;
    }
  }
}

void TreeSurvival::findBestSplitValueAUC(size_t nodeID, size_t varID, double& best_value, size_t& best_varID,
    double& best_auc) {

  // Create possible split values
  std::vector<double> possible_split_values;
  data->getAllValues(possible_split_values, sampleIDs[nodeID], varID);

  // Try next variable if all equal for this
  if (possible_split_values.size() < 2) {
    return;
  }

  size_t num_node_samples = sampleIDs[nodeID].size();
  size_t num_splits = possible_split_values.size() - 1;
  size_t num_possible_pairs = num_node_samples * (num_node_samples - 1) / 2;

  // Initialize
  std::vector<double> num_count(num_splits, num_possible_pairs);
  std::vector<double> num_total(num_splits, num_possible_pairs);
  std::vector<size_t> num_samples_left_child(num_splits);

  // For all pairs
  for (size_t k = 0; k < num_node_samples; ++k) {
    size_t sample_k = sampleIDs[nodeID][k];
    double time_k = data->get(sample_k, dependent_varID);
    double status_k = data->get(sample_k, status_varID);
    double value_k = data->get(sample_k, varID);

    // Count samples in left node
    for (size_t i = 0; i < num_splits; ++i) {
      double split_value = possible_split_values[i];
      if (value_k <= split_value) {
        ++num_samples_left_child[i];
      }
    }

    for (size_t l = k + 1; l < num_node_samples; ++l) {
      size_t sample_l = sampleIDs[nodeID][l];
      double time_l = data->get(sample_l, dependent_varID);
      double status_l = data->get(sample_l, status_varID);
      double value_l = data->get(sample_l, varID);

      // Compute split
      computeAucSplit(time_k, time_l, status_k, status_l, value_k, value_l, num_splits, possible_split_values,
          num_count, num_total);
    }
  }

  for (size_t i = 0; i < num_splits; ++i) {
    // Do not consider this split point if fewer than min_node_size samples in one node
    size_t num_samples_right_child = num_node_samples - num_samples_left_child[i];
    if (num_samples_left_child[i] < min_node_size || num_samples_right_child < min_node_size) {
      continue;
    } else {
      double auc = fabs((num_count[i] / 2) / num_total[i] - 0.5);
      if (auc > best_auc) {
        best_value = (possible_split_values[i] + possible_split_values[i + 1]) / 2;
        best_varID = varID;
        best_auc = auc;

        // Use smaller value if average is numerically the same as the larger value
        if (best_value == possible_split_values[i + 1]) {
          best_value = possible_split_values[i];
        }
      }
    }
  }
}

void TreeSurvival::computeAucSplit(double time_k, double time_l, double status_k, double status_l, double value_k,
    double value_l, size_t num_splits, std::vector<double>& possible_split_values, std::vector<double>& num_count,
    std::vector<double>& num_total) {

  bool ignore_pair = false;
  bool do_nothing = false;

  double value_smaller = 0;
  double value_larger = 0;
  double status_smaller = 0;

  if (time_k < time_l) {
    value_smaller = value_k;
    value_larger = value_l;
    status_smaller = status_k;
  } else if (time_l < time_k) {
    value_smaller = value_l;
    value_larger = value_k;
    status_smaller = status_l;
  } else {
    // Tie in survival time
    if (status_k == 0 || status_l == 0) {
      ignore_pair = true;
    } else {
      if (splitrule == AUC_IGNORE_TIES) {
        ignore_pair = true;
      } else {
        if (value_k == value_l) {
          // Tie in survival time and in covariate
          ignore_pair = true;
        } else {
          // Tie in survival time in covariate
          do_nothing = true;
        }
      }
    }
  }

  // Do not count if smaller time censored
  if (status_smaller == 0) {
    ignore_pair = true;
  }

  if (ignore_pair) {
    for (size_t i = 0; i < num_splits; ++i) {
      --num_count[i];
      --num_total[i];
    }
  } else if (do_nothing) {
    // Do nothing
  } else {
    for (size_t i = 0; i < num_splits; ++i) {
      double split_value = possible_split_values[i];

      if (value_smaller <= split_value && value_larger > split_value) {
        ++num_count[i];
      } else if (value_smaller > split_value && value_larger <= split_value) {
        --num_count[i];
      } else if (value_smaller <= split_value && value_larger <= split_value) {
        break;
      }
    }
  }

}

bool TreeSurvival::findBestSplitExtraTrees(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  double best_decrease = -1;
  size_t num_samples_node = sampleIDs[nodeID].size();
  size_t best_varID = 0;
  double best_value = 0;

  computeDeathCounts(nodeID);

  // Stop early if no split posssible
  if (num_samples_node >= 2 * min_node_size) {

    // For all possible split variables
    for (auto& varID : possible_split_varIDs) {

      // Find best split value, if ordered consider all values as split values, else all 2-partitions
      if (data->isOrderedVariable(varID)) {
        findBestSplitValueExtraTrees(nodeID, varID, best_value, best_varID, best_decrease);
      } else {
        findBestSplitValueExtraTreesUnordered(nodeID, varID, best_value, best_varID, best_decrease);
      }

    }
  }

  // Stop and save CHF if no good split found (this is terminal node).
  if (best_decrease < 0) {
    computeSurvival(nodeID);
    return true;
  } else {
    // If not terminal node save best values
    split_varIDs[nodeID] = best_varID;
    split_values[nodeID] = best_value;

    // Compute decrease of impurity for this node and add to variable importance if needed
    if (importance_mode == IMP_GINI || importance_mode == IMP_GINI_CORRECTED) {
      addImpurityImportance(nodeID, best_varID, best_decrease);
    }

    return false;
  }
}

void TreeSurvival::findBestSplitValueExtraTrees(size_t nodeID, size_t varID, double& best_value, size_t& best_varID,
    double& best_logrank) {

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

  size_t num_splits = possible_split_values.size();

  // Initialize
  std::vector<size_t> num_deaths_right_child(num_splits * num_timepoints);
  std::vector<size_t> delta_samples_at_risk_right_child(num_splits * num_timepoints);
  std::vector<size_t> num_samples_right_child(num_splits);

  computeChildDeathCounts(nodeID, varID, possible_split_values, num_samples_right_child,
      delta_samples_at_risk_right_child, num_deaths_right_child, num_splits);

  // Compute logrank test for all splits and use best
  for (size_t i = 0; i < num_splits; ++i) {
    double numerator = 0;
    double denominator_squared = 0;

    // Stop if minimal node size reached
    size_t num_samples_left_child = sampleIDs[nodeID].size() - num_samples_right_child[i];
    if (num_samples_right_child[i] < min_node_size || num_samples_left_child < min_node_size) {
      continue;
    }

    // Compute logrank test statistic for this split
    size_t num_samples_at_risk_right_child = num_samples_right_child[i];
    for (size_t t = 0; t < num_timepoints; ++t) {
      if (num_samples_at_risk[t] < 2 || num_samples_at_risk_right_child < 1) {
        break;
      }

      if (num_deaths[t] > 0) {
        // Numerator and demoninator for log-rank test, notation from Ishwaran et al.
        double di = (double) num_deaths[t];
        double di1 = (double) num_deaths_right_child[i * num_timepoints + t];
        double Yi = (double) num_samples_at_risk[t];
        double Yi1 = (double) num_samples_at_risk_right_child;
        numerator += di1 - Yi1 * (di / Yi);
        denominator_squared += (Yi1 / Yi) * (1.0 - Yi1 / Yi) * ((Yi - di) / (Yi - 1)) * di;
      }

      // Reduce number of samples at risk for next timepoint
      num_samples_at_risk_right_child -= delta_samples_at_risk_right_child[i * num_timepoints + t];

    }
    double logrank = -1;
    if (denominator_squared != 0) {
      logrank = fabs(numerator / sqrt(denominator_squared));
    }

    if (logrank > best_logrank) {
      best_value = possible_split_values[i];
      best_varID = varID;
      best_logrank = logrank;
    }
  }
}

void TreeSurvival::findBestSplitValueExtraTreesUnordered(size_t nodeID, size_t varID, double& best_value,
    size_t& best_varID, double& best_logrank) {

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
    std::vector<size_t> num_deaths_right_child(num_timepoints);
    std::vector<size_t> delta_samples_at_risk_right_child(num_timepoints);
    size_t num_samples_right_child = 0;
    double numerator = 0;
    double denominator_squared = 0;

    // Count deaths in right child per timepoint
    for (auto& sampleID : sampleIDs[nodeID]) {
      size_t survival_timeID = (*response_timepointIDs)[sampleID];
      double value = data->get(sampleID, varID);
      size_t factorID = floor(value) - 1;

      // If in right child, count
      // In right child, if bitwise splitID at position factorID is 1
      if ((splitID & (1 << factorID))) {
        ++num_samples_right_child;
        ++delta_samples_at_risk_right_child[survival_timeID];
        if (data->get(sampleID, status_varID) == 1) {
          ++num_deaths_right_child[survival_timeID];
        }
      }

    }

    // Stop if minimal node size reached
    size_t num_samples_left_child = sampleIDs[nodeID].size() - num_samples_right_child;
    if (num_samples_right_child < min_node_size || num_samples_left_child < min_node_size) {
      continue;
    }

    // Compute logrank test statistic for this split
    size_t num_samples_at_risk_right_child = num_samples_right_child;
    for (size_t t = 0; t < num_timepoints; ++t) {
      if (num_samples_at_risk[t] < 2 || num_samples_at_risk_right_child < 1) {
        break;
      }

      if (num_deaths[t] > 0) {
        // Numerator and demoninator for log-rank test, notation from Ishwaran et al.
        double di = (double) num_deaths[t];
        double di1 = (double) num_deaths_right_child[t];
        double Yi = (double) num_samples_at_risk[t];
        double Yi1 = (double) num_samples_at_risk_right_child;
        numerator += di1 - Yi1 * (di / Yi);
        denominator_squared += (Yi1 / Yi) * (1.0 - Yi1 / Yi) * ((Yi - di) / (Yi - 1)) * di;
      }

      // Reduce number of samples at risk for next timepoint
      num_samples_at_risk_right_child -= delta_samples_at_risk_right_child[t];
    }
    double logrank = -1;
    if (denominator_squared != 0) {
      logrank = fabs(numerator / sqrt(denominator_squared));
    }

    if (logrank > best_logrank) {
      best_value = splitID;
      best_varID = varID;
      best_logrank = logrank;
    }
  }
}

void TreeSurvival::addImpurityImportance(size_t nodeID, size_t varID, double decrease) {

  // No variable importance for no split variables
  size_t tempvarID = data->getUnpermutedVarID(varID);
  for (auto& skip : data->getNoSplitVariables()) {
    if (tempvarID >= skip) {
      --tempvarID;
    }
  }

  // Subtract if corrected importance and permuted variable, else add
  if (importance_mode == IMP_GINI_CORRECTED && varID >= data->getNumCols()) {
    (*variable_importance)[tempvarID] -= decrease;
  } else {
    (*variable_importance)[tempvarID] += decrease;
  }
}

} // namespace ranger
