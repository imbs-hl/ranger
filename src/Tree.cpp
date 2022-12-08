/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#include <iterator>

#include "Tree.h"
#include "utility.h"

namespace ranger {

Tree::Tree() :
    mtry(0), num_samples(0), num_samples_oob(0), min_node_size(0), deterministic_varIDs(0), split_select_weights(0), case_weights(
        0), manual_inbag(0), oob_sampleIDs(0), holdout(false), keep_inbag(false), data(0), regularization_factor(0), regularization_usedepth(
        false), split_varIDs_used(0), variable_importance(0), importance_mode(DEFAULT_IMPORTANCE_MODE), sample_with_replacement(
        true), sample_fraction(0), memory_saving_splitting(false), splitrule(DEFAULT_SPLITRULE), alpha(DEFAULT_ALPHA), minprop(
        DEFAULT_MINPROP), num_random_splits(DEFAULT_NUM_RANDOM_SPLITS), max_depth(DEFAULT_MAXDEPTH), depth(0), last_left_nodeID(
        0) {
}

Tree::Tree(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
    std::vector<double>& split_values) :
    mtry(0), num_samples(0), num_samples_oob(0), min_node_size(0), deterministic_varIDs(0), split_select_weights(0), case_weights(
        0), manual_inbag(0), split_varIDs(split_varIDs), split_values(split_values), child_nodeIDs(child_nodeIDs), oob_sampleIDs(
        0), holdout(false), keep_inbag(false), data(0), regularization_factor(0), regularization_usedepth(false), split_varIDs_used(
        0), variable_importance(0), importance_mode(DEFAULT_IMPORTANCE_MODE), sample_with_replacement(true), sample_fraction(
        0), memory_saving_splitting(false), splitrule(DEFAULT_SPLITRULE), alpha(DEFAULT_ALPHA), minprop(
        DEFAULT_MINPROP), num_random_splits(DEFAULT_NUM_RANDOM_SPLITS), max_depth(DEFAULT_MAXDEPTH), depth(0), last_left_nodeID(
        0) {
}

void Tree::init(const Data* data, uint mtry, size_t num_samples, uint seed, std::vector<size_t>* deterministic_varIDs,
    std::vector<double>* split_select_weights, ImportanceMode importance_mode, uint min_node_size,
    bool sample_with_replacement, bool memory_saving_splitting, SplitRule splitrule, std::vector<double>* case_weights,
    std::vector<size_t>* manual_inbag, bool keep_inbag, std::vector<double>* sample_fraction, double alpha,
    double minprop, bool holdout, uint num_random_splits, uint max_depth, std::vector<double>* regularization_factor,
    bool regularization_usedepth, std::vector<bool>* split_varIDs_used) {

  this->data = data;
  this->mtry = mtry;
  this->num_samples = num_samples;
  this->memory_saving_splitting = memory_saving_splitting;

  // Create root node, assign bootstrap sample and oob samples
  child_nodeIDs.push_back(std::vector<size_t>());
  child_nodeIDs.push_back(std::vector<size_t>());
  createEmptyNode();

  // Initialize random number generator and set seed
  random_number_generator.seed(seed);

  this->deterministic_varIDs = deterministic_varIDs;
  this->split_select_weights = split_select_weights;
  this->importance_mode = importance_mode;
  this->min_node_size = min_node_size;
  this->sample_with_replacement = sample_with_replacement;
  this->splitrule = splitrule;
  this->case_weights = case_weights;
  this->manual_inbag = manual_inbag;
  this->keep_inbag = keep_inbag;
  this->sample_fraction = sample_fraction;
  this->holdout = holdout;
  this->alpha = alpha;
  this->minprop = minprop;
  this->num_random_splits = num_random_splits;
  this->max_depth = max_depth;
  this->regularization_factor = regularization_factor;
  this->regularization_usedepth = regularization_usedepth;
  this->split_varIDs_used = split_varIDs_used;

  // Regularization
  if (regularization_factor->size() > 0) {
    regularization = true;
  } else {
    regularization = false;
  }
}

void Tree::grow(std::vector<double>* variable_importance) {
  // Allocate memory for tree growing
  allocateMemory();

  this->variable_importance = variable_importance;

  // Bootstrap, dependent if weighted or not and with or without replacement
  if (!case_weights->empty()) {
    if (sample_with_replacement) {
      bootstrapWeighted();
    } else {
      bootstrapWithoutReplacementWeighted();
    }
  } else if (sample_fraction->size() > 1) {
    if (sample_with_replacement) {
      bootstrapClassWise();
    } else {
      bootstrapWithoutReplacementClassWise();
    }
  } else if (!manual_inbag->empty()) {
    setManualInbag();
  } else {
    if (sample_with_replacement) {
      bootstrap();
    } else {
      bootstrapWithoutReplacement();
    }
  }

  // Init start and end positions
  start_pos[0] = 0;
  end_pos[0] = sampleIDs.size();

  // While not all nodes terminal, split next node
  size_t num_open_nodes = 1;
  size_t i = 0;
  depth = 0;
  while (num_open_nodes > 0) {
    // Split node
    bool is_terminal_node = splitNode(i);
    if (is_terminal_node) {
      --num_open_nodes;
    } else {
      ++num_open_nodes;
      if (i >= last_left_nodeID) {
        // If new level, increase depth
        // (left_node saves left-most node in current level, new level reached if that node is splitted)
        last_left_nodeID = split_varIDs.size() - 2;
        ++depth;
      }
    }
    ++i;
  }

  // Delete sampleID vector to save memory
  if (importance_mode != IMP_SOBOL_MDA){
    sampleIDs.clear();
    sampleIDs.shrink_to_fit();
  }
  cleanUpInternal();
  
}

void Tree::predict(const Data* prediction_data, bool oob_prediction) {

  size_t num_samples_predict;
  if (oob_prediction) {
    num_samples_predict = num_samples_oob;
  } else {
    num_samples_predict = prediction_data->getNumRows();
  }

  prediction_terminal_nodeIDs.resize(num_samples_predict, 0);

  // For each sample start in root, drop down the tree and return final value
  for (size_t i = 0; i < num_samples_predict; ++i) {
    size_t sample_idx;
    if (oob_prediction) {
      sample_idx = oob_sampleIDs[i];
    } else {
      sample_idx = i;
    }
    size_t nodeID = 0;
    while (1) {

      // Break if terminal node
      if (child_nodeIDs[0][nodeID] == 0 && child_nodeIDs[1][nodeID] == 0) {
        break;
      }

      // Move to child
      size_t split_varID = split_varIDs[nodeID];

      double value = prediction_data->get_x(sample_idx, split_varID);
      if (prediction_data->isOrderedVariable(split_varID)) {
        if (value <= split_values[nodeID]) {
          // Move to left child
          nodeID = child_nodeIDs[0][nodeID];
        } else {
          // Move to right child
          nodeID = child_nodeIDs[1][nodeID];
        }
      } else {
        size_t factorID = floor(value) - 1;
        size_t splitID = floor(split_values[nodeID]);

        // Left if 0 found at position factorID
        if (!(splitID & (1ULL << factorID))) {
          // Move to left child
          nodeID = child_nodeIDs[0][nodeID];
        } else {
          // Move to right child
          nodeID = child_nodeIDs[1][nodeID];
        }
      }
    }

    prediction_terminal_nodeIDs[i] = nodeID;
  }
}

void Tree::computePermutationImportance(std::vector<double>& forest_importance, std::vector<double>& forest_variance,
    std::vector<double>& forest_importance_casewise) {

  size_t num_independent_variables = data->getNumCols();

  // Compute normal prediction accuracy for each tree. Predictions already computed..
  double accuracy_normal;
  std::vector<double> prederr_normal_casewise;
  std::vector<double> prederr_shuf_casewise;
  if (importance_mode == IMP_PERM_CASEWISE) {
    prederr_normal_casewise.resize(num_samples_oob, 0);
    prederr_shuf_casewise.resize(num_samples_oob, 0);
    accuracy_normal = computePredictionAccuracyInternal(&prederr_normal_casewise);
  } else {
    accuracy_normal = computePredictionAccuracyInternal(NULL);
  }

  prediction_terminal_nodeIDs.clear();
  prediction_terminal_nodeIDs.resize(num_samples_oob, 0);

  // Reserve space for permutations, initialize with oob_sampleIDs
  std::vector<size_t> permutations(oob_sampleIDs);

  // Randomly permute for all independent variables
  for (size_t i = 0; i < num_independent_variables; ++i) {

    // Check whether the i-th variable is used in the
	  // tree:
    bool isused = false;
    for (size_t j = 0; j < split_varIDs.size(); ++j) {
      if (split_varIDs[j] == i) {
        isused = true;
        break;
      }
    }
     
	 // Only do permutations if the variable is used in the tree, otherwise variable importance is 0
    if (isused) {
      // Permute and compute prediction accuracy again for this permutation and save difference
      permuteAndPredictOobSamples(i, permutations);
      double accuracy_permuted;
      if (importance_mode == IMP_PERM_CASEWISE) {
        accuracy_permuted = computePredictionAccuracyInternal(&prederr_shuf_casewise);
        for (size_t j = 0; j < num_samples_oob; ++j) {
          size_t pos = i * num_samples + oob_sampleIDs[j];
          forest_importance_casewise[pos] += prederr_shuf_casewise[j] - prederr_normal_casewise[j];
        }
      } else {
        accuracy_permuted = computePredictionAccuracyInternal(NULL);
      }
  
      double accuracy_difference = accuracy_normal - accuracy_permuted;
      forest_importance[i] += accuracy_difference;
  
      // Compute variance
      if (importance_mode == IMP_PERM_BREIMAN) {
        forest_variance[i] += accuracy_difference * accuracy_difference;
      } else if (importance_mode == IMP_PERM_LIAW) {
        forest_variance[i] += accuracy_difference * accuracy_difference * num_samples_oob;
      }
    }
  }
}

// Compute predictions of projected trees for all input variables.
std::vector<std::vector<double>> Tree::predictProjectedTree() {
  
  // Initialize variables.
  size_t p = data->getNumCols(); // Set p as the number of input variables.
  // The vectors nodeIDs_samplesIDs_inb and nodeIDs_samplesIDs_oob are of size p, and respectively store the in-bag and out-of-bag sampleIDs of the projected partitions.
  // The j-th component of these two vectors stores the projected partitions of the tree at all levels when the splits on variable j are ignored (i.e., observations are sent to both children nodes).
  // Each projected partition is a map: the key is a list of nodeIDs, and the value is the list of the sampleIDs falling simultaneously in all nodes of the key.
  std::vector<std::map<std::vector<size_t>, std::vector<size_t>>> nodeIDs_samplesIDs_inb(p); 
  std::vector<std::map<std::vector<size_t>, std::vector<size_t>>> nodeIDs_samplesIDs_oob(p);
  std::map<size_t, std::vector<size_t>> leaves_variables; // Lists of variables used in the sequence of splits to reach each terminal leave.
  std::map<size_t, std::vector<size_t>> leaves_sampleIDs; // Lists of sampleIDs falling in each terminal leave.
  
  // Deduplication of in-bag sample to speed up computations.
  std::vector<size_t> sampleIDs_unique; // Deduplicated sampleIDs.
  std::map<size_t, size_t> sampleIDs_count; // Number of repetitions of each observation of sampleIDs_unique in sampleIDs.
  for (auto & i : sampleIDs){
    std::pair<std::map<size_t, size_t>::iterator, bool> itIDs = sampleIDs_count.insert(std::pair<size_t, size_t> (i, 1));
    if (!itIDs.second){
      itIDs.first->second += 1;
    }else{
      sampleIDs_unique.push_back(i);
    }
  }
  
  // Drop observations in the projected partitions. First consider in-bag sample (s = 0), and then the out-of-bag sample (s = 1).
  for (size_t s = 0; s <= 1; ++s){
    std::vector<size_t> sampleIDs_temp;
    if (s == 0){
      sampleIDs_temp = sampleIDs_unique; // Deduplicated in-bag sample.
    } else {
      sampleIDs_temp = oob_sampleIDs; // Out-of-bag sample.
    }
    // Loop through all observations of the considered sample.
    for (auto & i : sampleIDs_temp) {
      std::vector<size_t> split_variables; // Track the list of variables involved in the splits of the previous tree levels.
      bool is_terminal = false; // Track whether terminal leave is reached.
      size_t nodeIDorigin = 0; // nodeIDorigin is the nodeID where observation i falls at each level of the original tree. Initialized here at the root node.
      // Loop through the tree levels, and stop when observation i has reached the terminal leave of the original tree.
      while (!is_terminal){
        if (child_nodeIDs[0][nodeIDorigin] == 0 && child_nodeIDs[1][nodeIDorigin] == 0) {
          is_terminal = true; // If nodeIDorigin is a terminal leave, break the loop through the levels of the original tree.
        }else{
          std::vector<size_t> nodeIDtemp; // Store the list of nodeIDs at each level when observation i is dropped down the projected tree.
          nodeIDtemp.push_back(nodeIDorigin);
          size_t varIDorigin = split_varIDs[nodeIDorigin]; // varIDorigin is the splitting variable at nodeIDorigin.
          nodeIDorigin = getChildID(i, nodeIDorigin, varIDorigin); // Update nodeIDorigin as the child node where observation i falls at the next level of the original tree.
          // If varIDorigin is NOT already involved in the nodes hit by observation i at the previous tree levels, 
          // then compute the projected cell where observation i falls, by ignoring all splits using varIDorigin down the tree, otherwise send observation i to the next level of the original tree.
          if (find(split_variables.begin(), split_variables.end(), varIDorigin) == split_variables.end()){
            split_variables.push_back(varIDorigin); // Add varIDorigin to the list of variables already met by observation i at the previous tree levels.
            std::vector<size_t> nodeIDnext; // List of inner nodeIDs of the next tree level where observation i falls.
            std::vector<size_t> nodeIDpred; // List of terminal nodeIDs where observation i falls.
            std::vector<size_t> nodeIDpredtemp = nodeIDtemp; // Store projected cell of the previous level of the projected tree.
            bool next_level = true;
            // Loop through the remaining tree levels, ignoring all splits using varIDorigin, i.e.,
            // observation i is sent to both children nodes when a split on varIDorigin is hit, and then, observation i ends up in multiple terminal leaves.
            while (next_level) {
              for (auto & nodeID : nodeIDtemp){
                if (child_nodeIDs[0][nodeID] == 0 && child_nodeIDs[1][nodeID] == 0) {
                  nodeIDpred.push_back(nodeID);
                }else{
                  size_t split_varID = split_varIDs[nodeID];
                  if (split_varID != varIDorigin){
                    // Move to child node if splitting variable is NOT varIDorigin.
                    nodeIDnext.push_back(getChildID(i, nodeID, split_varID));
                  }else{
                    // Move to both left and right children nodes if splitting variable is varIDorigin.
                    nodeIDnext.push_back(child_nodeIDs[0][nodeID]);
                    nodeIDnext.push_back(child_nodeIDs[1][nodeID]);
                  }
                }
              }
              // Store observation i in the projected partition of the current tree level.
              if (s == 0){
                // In-bag observation case.
                if (nodeIDnext.size() > 0 || nodeIDtemp == std::vector<size_t> {0}){
                  std::vector<size_t> nodeIDlevel = nodeIDpred; // nodeIDlevel is the list of nodeIDs defining the projected cell at the current tree level.
                  nodeIDlevel.insert(nodeIDlevel.end(), nodeIDnext.begin(), nodeIDnext.end()); // Assemble nodeIDs of terminal leaves and inner nodes in nodeIDlevel.
                  std::pair<std::map<std::vector<size_t>, std::vector<size_t>>::iterator, bool> it_bool;
                  // Insert the list of nodeIDs defining the projected cell at the current tree level (nodeIDlevel) and the associated sampleIDs.
                  it_bool = nodeIDs_samplesIDs_inb[varIDorigin].insert(std::pair<std::vector<size_t>, std::vector<size_t>> (nodeIDlevel, {i}));
                  if (!it_bool.second){
                    it_bool.first->second.push_back(i);
                  }
                }
                nodeIDtemp = nodeIDnext;
                nodeIDnext = {};
                next_level = nodeIDtemp.size() > 0;
              }else{
                // Out-of-bag observation case.
                std::vector<size_t> nodeIDlevel = nodeIDpred;
                nodeIDlevel.insert(nodeIDlevel.end(), nodeIDnext.begin(), nodeIDnext.end()); // Assemble nodeIDs of terminal leaves and inner nodes.
                // Check whether the projected cell exists in the pre-computed in-bag projected partition. 
                std::map<std::vector<size_t>, std::vector<size_t>>::iterator it = nodeIDs_samplesIDs_inb[varIDorigin].find(nodeIDlevel);
                if (it == nodeIDs_samplesIDs_inb[varIDorigin].end() || nodeIDnext.size() == 0){
                  // If the projected cell is NOT in the in-bag projected partition or the projected cell is empty, 
                  // set the projected cell of the previous tree level (nodeIDpredtemp) as the terminal projected leave.
                  next_level = false;
                  std::pair<std::map<std::vector<size_t>, std::vector<size_t>>::iterator, bool> it_bool;
                  it_bool = nodeIDs_samplesIDs_oob[varIDorigin].insert(std::pair<std::vector<size_t>, std::vector<size_t>> (nodeIDpredtemp, {i}));
                  if (!it_bool.second){
                    it_bool.first->second.push_back(i);
                  }
                }else{
                  // If the projected cell is in the in-bag projected partition and not empty, keep splitting.
                  nodeIDpredtemp = nodeIDlevel;
                  nodeIDtemp = nodeIDnext;
                  nodeIDnext = {};
                }
              }
            }
          }
        }
      }
      // For out-of-bag observations, store terminal leave information.
      if (s == 1){
        // Store list of variables used in the sequence of splits.
        leaves_variables.insert(std::pair<size_t, std::vector<size_t>> (nodeIDorigin, split_variables));
        std::pair<std::map<size_t, std::vector<size_t>>::iterator, bool> it_leaves;
        // Store oob sampleIDs in terminal leaves.
        it_leaves = leaves_sampleIDs.insert(std::pair<size_t, std::vector<size_t>> (nodeIDorigin, {i}));
        if (!it_leaves.second){
          it_leaves.first->second.push_back(i);
        }
      }
    }
  }
  
  // Generate out-of-bag predictions from the obtained projected partitions.
  int nsample = data->getNumRows(); // Sample size.
  std::vector<std::vector<double>> predictions_oob; // Out-of-bag projected tree predictions for all input variables.
  std::map<std::vector<size_t>, double> nodeIDs_response;
  std::pair<std::map<std::vector<size_t>, double>::iterator, bool> it_response;
  for (size_t j = 0; j < p; ++j){
    std::vector<double> predictions_varID(nsample, 0); // Out-of-bag projected tree predictions when variable varID is ignored.
    for (std::map<std::vector<size_t>, std::vector<size_t>>::iterator it = nodeIDs_samplesIDs_oob[j].begin(); it != nodeIDs_samplesIDs_oob[j].end(); ++it){
      double response;
      std::vector<size_t> nodeIDs = it->first;
      it_response = nodeIDs_response.insert(std::pair<std::vector<size_t>, double> (nodeIDs, 0));
      // Compute projected cell output at its first occurence, using in-bag data.
      if (it_response.second){
        std::vector<size_t> sampleIDs_inb = nodeIDs_samplesIDs_inb[j][nodeIDs];
        response = 0;
        size_t inb_size = 0;
        for (auto & i : sampleIDs_inb){
          size_t size_temp = sampleIDs_count[i];
          response += data->get_y(i, 0)*size_temp;
          inb_size += size_temp;
        }
        response = response/inb_size;
        it_response.first->second = response;
      }else{
        response = it_response.first->second;
      }
      // Assign projected cell output to all oob observations falling in this projected cell.
      for (auto & i : it->second){
        predictions_varID[i] = response;
      }
    }
    predictions_oob.push_back(predictions_varID);
  }
  // For all variables not involved in the splits of the original tree path of a given observation, the projected prediction is given by the original tree.
  for (std::map<size_t, std::vector<size_t>>::iterator it = leaves_variables.begin(); it != leaves_variables.end(); ++it){
    double response = split_values[it->first];
    for (size_t j = 0; j < p; ++j){
      if (find(it->second.begin(), it->second.end(), j) == it->second.end()){
        for (auto & i : leaves_sampleIDs[it->first]){
          predictions_oob[j][i] = response;
        }
      }
    }
  }
  
  return(predictions_oob);
  
}

// Return child nodeID.
size_t Tree::getChildID(size_t i, size_t nodeID, size_t split_varID){
  size_t childID;
  double value = data->get_x(i, split_varID);
  if (data->isOrderedVariable(split_varID)) {
    if (value <= split_values[nodeID]) {
      // Move to left child.
      childID = child_nodeIDs[0][nodeID];
    } else {
      // Move to right child.
      childID = child_nodeIDs[1][nodeID];
    }
  } else {
    size_t factorID = floor(value) - 1;
    size_t splitID = floor(split_values[nodeID]);
    // Left if 0 found at position factorID.
    if (!(splitID & (1ULL << factorID))) {
      // Move to left child.
      childID = child_nodeIDs[0][nodeID];
    } else {
      // Move to right child.
      childID = child_nodeIDs[1][nodeID];
    }
  }
  return(childID);
}

// #nocov start
void Tree::appendToFile(std::ofstream& file) {

  // Save general fields
  saveVector2D(child_nodeIDs, file);
  saveVector1D(split_varIDs, file);
  saveVector1D(split_values, file);

  // Call special functions for subclasses to save special fields.
  appendToFileInternal(file);
}
// #nocov end

void Tree::createPossibleSplitVarSubset(std::vector<size_t>& result) {

  size_t num_vars = data->getNumCols();

  // For corrected Gini importance add dummy variables
  if (importance_mode == IMP_GINI_CORRECTED) {
    num_vars += data->getNumCols();
  }

  // Randomly add non-deterministic variables (according to weights if needed)
  if (split_select_weights->empty()) {
    if (deterministic_varIDs->empty()) {
      drawWithoutReplacement(result, random_number_generator, num_vars, mtry);
    } else {
      drawWithoutReplacementSkip(result, random_number_generator, num_vars, (*deterministic_varIDs), mtry);
    }
  } else {
    drawWithoutReplacementWeighted(result, random_number_generator, num_vars, mtry, *split_select_weights);
  }

  // Always use deterministic variables
  std::copy(deterministic_varIDs->begin(), deterministic_varIDs->end(), std::inserter(result, result.end()));
}

bool Tree::splitNode(size_t nodeID) {

  // Select random subset of variables to possibly split at
  std::vector<size_t> possible_split_varIDs;
  createPossibleSplitVarSubset(possible_split_varIDs);

  // Call subclass method, sets split_varIDs and split_values
  bool stop = splitNodeInternal(nodeID, possible_split_varIDs);
  if (stop) {
    // Terminal node
    return true;
  }

  size_t split_varID = split_varIDs[nodeID];
  double split_value = split_values[nodeID];

  // Save non-permuted variable for prediction
  split_varIDs[nodeID] = data->getUnpermutedVarID(split_varID);

  // Create child nodes
  size_t left_child_nodeID = split_varIDs.size();
  child_nodeIDs[0][nodeID] = left_child_nodeID;
  createEmptyNode();
  start_pos[left_child_nodeID] = start_pos[nodeID];

  size_t right_child_nodeID = split_varIDs.size();
  child_nodeIDs[1][nodeID] = right_child_nodeID;
  createEmptyNode();
  start_pos[right_child_nodeID] = end_pos[nodeID];

  // For each sample in node, assign to left or right child
  if (data->isOrderedVariable(split_varID)) {
    // Ordered: left is <= splitval and right is > splitval
    size_t pos = start_pos[nodeID];
    while (pos < start_pos[right_child_nodeID]) {
      size_t sampleID = sampleIDs[pos];
      if (data->get_x(sampleID, split_varID) <= split_value) {
        // If going to left, do nothing
        ++pos;
      } else {
        // If going to right, move to right end
        --start_pos[right_child_nodeID];
        std::swap(sampleIDs[pos], sampleIDs[start_pos[right_child_nodeID]]);
      }
    }
  } else {
    // Unordered: If bit at position is 1 -> right, 0 -> left
    size_t pos = start_pos[nodeID];
    while (pos < start_pos[right_child_nodeID]) {
      size_t sampleID = sampleIDs[pos];
      double level = data->get_x(sampleID, split_varID);
      size_t factorID = floor(level) - 1;
      size_t splitID = floor(split_value);

      // Left if 0 found at position factorID
      if (!(splitID & (1ULL << factorID))) {
        // If going to left, do nothing
        ++pos;
      } else {
        // If going to right, move to right end
        --start_pos[right_child_nodeID];
        std::swap(sampleIDs[pos], sampleIDs[start_pos[right_child_nodeID]]);
      }
    }
  }

  // End position of left child is start position of right child
  end_pos[left_child_nodeID] = start_pos[right_child_nodeID];
  end_pos[right_child_nodeID] = end_pos[nodeID];

  // No terminal node
  return false;
}

void Tree::createEmptyNode() {
  split_varIDs.push_back(0);
  split_values.push_back(0);
  child_nodeIDs[0].push_back(0);
  child_nodeIDs[1].push_back(0);
  start_pos.push_back(0);
  end_pos.push_back(0);

  createEmptyNodeInternal();
}

size_t Tree::dropDownSamplePermuted(size_t permuted_varID, size_t sampleID, size_t permuted_sampleID) {

  // Start in root and drop down
  size_t nodeID = 0;
  while (child_nodeIDs[0][nodeID] != 0 || child_nodeIDs[1][nodeID] != 0) {

    // Permute if variable is permutation variable
    size_t split_varID = split_varIDs[nodeID];
    size_t sampleID_final = sampleID;
    if (split_varID == permuted_varID) {
      sampleID_final = permuted_sampleID;
    }

    // Move to child
    double value = data->get_x(sampleID_final, split_varID);
    if (data->isOrderedVariable(split_varID)) {
      if (value <= split_values[nodeID]) {
        // Move to left child
        nodeID = child_nodeIDs[0][nodeID];
      } else {
        // Move to right child
        nodeID = child_nodeIDs[1][nodeID];
      }
    } else {
      size_t factorID = floor(value) - 1;
      size_t splitID = floor(split_values[nodeID]);

      // Left if 0 found at position factorID
      if (!(splitID & (1ULL << factorID))) {
        // Move to left child
        nodeID = child_nodeIDs[0][nodeID];
      } else {
        // Move to right child
        nodeID = child_nodeIDs[1][nodeID];
      }
    }

  }
  return nodeID;
}

void Tree::permuteAndPredictOobSamples(size_t permuted_varID, std::vector<size_t>& permutations) {

  // Permute OOB sample
  //std::vector<size_t> permutations(oob_sampleIDs);
  std::shuffle(permutations.begin(), permutations.end(), random_number_generator);

  // For each sample, drop down the tree and add prediction
  for (size_t i = 0; i < num_samples_oob; ++i) {
    size_t nodeID = dropDownSamplePermuted(permuted_varID, oob_sampleIDs[i], permutations[i]);
    prediction_terminal_nodeIDs[i] = nodeID;
  }
}

void Tree::bootstrap() {

  // Use fraction (default 63.21%) of the samples
  size_t num_samples_inbag = (size_t) num_samples * (*sample_fraction)[0];

  // Reserve space, reserve a little more to be save)
  sampleIDs.reserve(num_samples_inbag);
  oob_sampleIDs.reserve(num_samples * (exp(-(*sample_fraction)[0]) + 0.1));

  std::uniform_int_distribution<size_t> unif_dist(0, num_samples - 1);

  // Start with all samples OOB
  inbag_counts.resize(num_samples, 0);

  // Draw num_samples samples with replacement (num_samples_inbag out of n) as inbag and mark as not OOB
  for (size_t s = 0; s < num_samples_inbag; ++s) {
    size_t draw = unif_dist(random_number_generator);
    sampleIDs.push_back(draw);
    ++inbag_counts[draw];
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

void Tree::bootstrapWeighted() {

  // Use fraction (default 63.21%) of the samples
  size_t num_samples_inbag = (size_t) num_samples * (*sample_fraction)[0];

  // Reserve space, reserve a little more to be save)
  sampleIDs.reserve(num_samples_inbag);
  oob_sampleIDs.reserve(num_samples * (exp(-(*sample_fraction)[0]) + 0.1));

  std::discrete_distribution<> weighted_dist(case_weights->begin(), case_weights->end());

  // Start with all samples OOB
  inbag_counts.resize(num_samples, 0);

  // Draw num_samples samples with replacement (n out of n) as inbag and mark as not OOB
  for (size_t s = 0; s < num_samples_inbag; ++s) {
    size_t draw = weighted_dist(random_number_generator);
    sampleIDs.push_back(draw);
    ++inbag_counts[draw];
  }

  // Save OOB samples. In holdout mode these are the cases with 0 weight.
  if (holdout) {
    for (size_t s = 0; s < (*case_weights).size(); ++s) {
      if ((*case_weights)[s] == 0) {
        oob_sampleIDs.push_back(s);
      }
    }
  } else {
    for (size_t s = 0; s < inbag_counts.size(); ++s) {
      if (inbag_counts[s] == 0) {
        oob_sampleIDs.push_back(s);
      }
    }
  }
  num_samples_oob = oob_sampleIDs.size();

  if (!keep_inbag) {
    inbag_counts.clear();
    inbag_counts.shrink_to_fit();
  }
}

void Tree::bootstrapWithoutReplacement() {

  // Use fraction (default 63.21%) of the samples
  size_t num_samples_inbag = (size_t) num_samples * (*sample_fraction)[0];
  shuffleAndSplit(sampleIDs, oob_sampleIDs, num_samples, num_samples_inbag, random_number_generator);
  num_samples_oob = oob_sampleIDs.size();

  if (keep_inbag) {
    // All observation are 0 or 1 times inbag
    inbag_counts.resize(num_samples, 1);
    for (size_t i = 0; i < oob_sampleIDs.size(); i++) {
      inbag_counts[oob_sampleIDs[i]] = 0;
    }
  }
}

void Tree::bootstrapWithoutReplacementWeighted() {

  // Use fraction (default 63.21%) of the samples
  size_t num_samples_inbag = (size_t) num_samples * (*sample_fraction)[0];
  drawWithoutReplacementWeighted(sampleIDs, random_number_generator, num_samples - 1, num_samples_inbag, *case_weights);

  // All observation are 0 or 1 times inbag
  inbag_counts.resize(num_samples, 0);
  for (auto& sampleID : sampleIDs) {
    inbag_counts[sampleID] = 1;
  }

  // Save OOB samples. In holdout mode these are the cases with 0 weight.
  if (holdout) {
    for (size_t s = 0; s < (*case_weights).size(); ++s) {
      if ((*case_weights)[s] == 0) {
        oob_sampleIDs.push_back(s);
      }
    }
  } else {
    for (size_t s = 0; s < inbag_counts.size(); ++s) {
      if (inbag_counts[s] == 0) {
        oob_sampleIDs.push_back(s);
      }
    }
  }
  num_samples_oob = oob_sampleIDs.size();

  if (!keep_inbag) {
    inbag_counts.clear();
    inbag_counts.shrink_to_fit();
  }
}

void Tree::bootstrapClassWise() {
  // Empty on purpose (virtual function only implemented in classification and probability)
}

void Tree::bootstrapWithoutReplacementClassWise() {
  // Empty on purpose (virtual function only implemented in classification and probability)
}

void Tree::setManualInbag() {
  // Select observation as specified in manual_inbag vector
  sampleIDs.reserve(manual_inbag->size());
  inbag_counts.resize(num_samples, 0);
  for (size_t i = 0; i < manual_inbag->size(); ++i) {
    size_t inbag_count = (*manual_inbag)[i];
    if ((*manual_inbag)[i] > 0) {
      for (size_t j = 0; j < inbag_count; ++j) {
        sampleIDs.push_back(i);
      }
      inbag_counts[i] = inbag_count;
    } else {
      oob_sampleIDs.push_back(i);
    }
  }
  num_samples_oob = oob_sampleIDs.size();

  // Shuffle samples
  std::shuffle(sampleIDs.begin(), sampleIDs.end(), random_number_generator);

  if (!keep_inbag) {
    inbag_counts.clear();
    inbag_counts.shrink_to_fit();
  }
}

} // namespace ranger
