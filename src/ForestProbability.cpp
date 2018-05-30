/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#include <stdexcept>

#include "utility.h"
#include "ForestProbability.h"
#include "TreeProbability.h"
#include "Data.h"

namespace ranger {

void ForestProbability::loadForest(size_t dependent_varID, size_t num_trees,
    std::vector<std::vector<std::vector<size_t>> >& forest_child_nodeIDs,
    std::vector<std::vector<size_t>>& forest_split_varIDs, std::vector<std::vector<double>>& forest_split_values,
    std::vector<double>& class_values, std::vector<std::vector<std::vector<double>>>& forest_terminal_class_counts,
    std::vector<bool>& is_ordered_variable) {

  this->dependent_varID = dependent_varID;
  this->num_trees = num_trees;
  this->class_values = class_values;
  data->setIsOrderedVariable(is_ordered_variable);

  // Create trees
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    trees.push_back(
        make_unique<TreeProbability>(forest_child_nodeIDs[i], forest_split_varIDs[i], forest_split_values[i],
            &this->class_values, &response_classIDs, forest_terminal_class_counts[i]));
  }

  // Create thread ranges
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);
}

std::vector<std::vector<std::vector<double>>> ForestProbability::getTerminalClassCounts() const {
  std::vector<std::vector<std::vector<double>>> result;
  result.reserve(num_trees);
  for (const auto& tree : trees) {
    const auto& temp = dynamic_cast<const TreeProbability&>(*tree);
    result.push_back(temp.getTerminalClassCounts());
  }
  return result;
}

void ForestProbability::initInternal(std::string status_variable_name) {

  // If mtry not set, use floored square root of number of independent variables.
  if (mtry == 0) {
    unsigned long temp = sqrt((double) (num_variables - 1));
    mtry = std::max((unsigned long) 1, temp);
  }

  // Set minimal node size
  if (min_node_size == 0) {
    min_node_size = DEFAULT_MIN_NODE_SIZE_PROBABILITY;
  }

  // Create class_values and response_classIDs
  if (!prediction_mode) {
    for (size_t i = 0; i < num_samples; ++i) {
      double value = data->get(i, dependent_varID);

      // If classID is already in class_values, use ID. Else create a new one.
      uint classID = find(class_values.begin(), class_values.end(), value) - class_values.begin();
      if (classID == class_values.size()) {
        class_values.push_back(value);
      }
      response_classIDs.push_back(classID);
    }
  }

  // Create sampleIDs_per_class if required
  if (sample_fraction.size() > 1) {
    sampleIDs_per_class.resize(sample_fraction.size());
    for (auto& v : sampleIDs_per_class) {
      v.reserve(num_samples);
    }
    for (size_t i = 0; i < num_samples; ++i) {
      size_t classID = response_classIDs[i];
      sampleIDs_per_class[classID].push_back(i);
    }
  }

  // Set class weights all to 1
  class_weights = std::vector<double>(class_values.size(), 1.0);

  // Sort data if memory saving mode
  if (!memory_saving_splitting) {
    data->sort();
  }
}

void ForestProbability::growInternal() {
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    trees.push_back(
        make_unique<TreeProbability>(&class_values, &response_classIDs, &sampleIDs_per_class, &class_weights));
  }
}

void ForestProbability::allocatePredictMemory() {
  size_t num_prediction_samples = data->getNumRows();
  if (predict_all) {
    predictions = std::vector<std::vector<std::vector<double>>>(num_prediction_samples,
        std::vector<std::vector<double>>(class_values.size(), std::vector<double>(num_trees, 0)));
  } else if (prediction_type == TERMINALNODES) {
    predictions = std::vector<std::vector<std::vector<double>>>(1,
        std::vector<std::vector<double>>(num_prediction_samples, std::vector<double>(num_trees, 0)));
  } else {
    predictions = std::vector<std::vector<std::vector<double>>>(1,
        std::vector<std::vector<double>>(num_prediction_samples, std::vector<double>(class_values.size(), 0)));
  }
}

void ForestProbability::predictInternal(size_t sample_idx) {
  // For each sample compute proportions in each tree
  for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
    if (predict_all) {
      std::vector<double> counts = getTreePrediction(tree_idx, sample_idx);

      for (size_t class_idx = 0; class_idx < counts.size(); ++class_idx) {
        predictions[sample_idx][class_idx][tree_idx] += counts[class_idx];
      }
    } else if (prediction_type == TERMINALNODES) {
      predictions[0][sample_idx][tree_idx] = getTreePredictionTerminalNodeID(tree_idx, sample_idx);
    } else {
      std::vector<double> counts = getTreePrediction(tree_idx, sample_idx);

      for (size_t class_idx = 0; class_idx < counts.size(); ++class_idx) {
        predictions[0][sample_idx][class_idx] += counts[class_idx];
      }
    }
  }

  // Average over trees
  if (!predict_all && prediction_type != TERMINALNODES) {
    for (size_t class_idx = 0; class_idx < predictions[0][sample_idx].size(); ++class_idx) {
      predictions[0][sample_idx][class_idx] /= num_trees;
    }
  }
}

void ForestProbability::computePredictionErrorInternal() {

// For each sample sum over trees where sample is OOB
  std::vector<size_t> samples_oob_count;
  samples_oob_count.resize(num_samples, 0);
  predictions = std::vector<std::vector<std::vector<double>>>(1,
      std::vector<std::vector<double>>(num_samples, std::vector<double>(class_values.size(), 0)));

  for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
    for (size_t sample_idx = 0; sample_idx < trees[tree_idx]->getNumSamplesOob(); ++sample_idx) {
      size_t sampleID = trees[tree_idx]->getOobSampleIDs()[sample_idx];
      std::vector<double> counts = getTreePrediction(tree_idx, sample_idx);

      for (size_t class_idx = 0; class_idx < counts.size(); ++class_idx) {
        predictions[0][sampleID][class_idx] += counts[class_idx];
      }
      ++samples_oob_count[sampleID];
    }
  }

// MSE with predicted probability and true data
  size_t num_predictions = 0;
  overall_prediction_error = 0;
  for (size_t i = 0; i < predictions[0].size(); ++i) {
    if (samples_oob_count[i] > 0) {
      ++num_predictions;
      for (size_t j = 0; j < predictions[0][i].size(); ++j) {
        predictions[0][i][j] /= (double) samples_oob_count[i];
      }
      size_t real_classID = response_classIDs[i];
      double predicted_value = predictions[0][i][real_classID];
      overall_prediction_error += (1 - predicted_value) * (1 - predicted_value);
    } else {
      for (size_t j = 0; j < predictions[0][i].size(); ++j) {
        predictions[0][i][j] = NAN;
      }
    }
  }

  overall_prediction_error /= (double) num_predictions;
}

// #nocov start
void ForestProbability::writeOutputInternal() {
  if (verbose_out) {
    *verbose_out << "Tree type:                         " << "Probability estimation" << std::endl;
  }
}

void ForestProbability::writeConfusionFile() {

// Open confusion file for writing
  std::string filename = output_prefix + ".confusion";
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to confusion file: " + filename + ".");
  }

// Write confusion to file
  outfile << "Overall OOB prediction error (MSE): " << overall_prediction_error << std::endl;

  outfile.close();
  if (verbose_out)
    *verbose_out << "Saved prediction error to file " << filename << "." << std::endl;
}

void ForestProbability::writePredictionFile() {

  // Open prediction file for writing
  std::string filename = output_prefix + ".prediction";
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to prediction file: " + filename + ".");
  }

  // Write
  outfile << "Class predictions, one sample per row." << std::endl;
  for (auto& class_value : class_values) {
    outfile << class_value << " ";
  }
  outfile << std::endl << std::endl;

  if (predict_all) {
    for (size_t k = 0; k < num_trees; ++k) {
      outfile << "Tree " << k << ":" << std::endl;
      for (size_t i = 0; i < predictions.size(); ++i) {
        for (size_t j = 0; j < predictions[i].size(); ++j) {
          outfile << predictions[i][j][k] << " ";
        }
        outfile << std::endl;
      }
      outfile << std::endl;
    }
  } else {
    for (size_t i = 0; i < predictions.size(); ++i) {
      for (size_t j = 0; j < predictions[i].size(); ++j) {
        for (size_t k = 0; k < predictions[i][j].size(); ++k) {
          outfile << predictions[i][j][k] << " ";
        }
        outfile << std::endl;
      }
    }
  }

  if (verbose_out)
    *verbose_out << "Saved predictions to file " << filename << "." << std::endl;
}

void ForestProbability::saveToFileInternal(std::ofstream& outfile) {

// Write num_variables
  outfile.write((char*) &num_variables, sizeof(num_variables));

// Write treetype
  TreeType treetype = TREE_PROBABILITY;
  outfile.write((char*) &treetype, sizeof(treetype));

// Write class_values
  saveVector1D(class_values, outfile);
}

void ForestProbability::loadFromFileInternal(std::ifstream& infile) {

// Read number of variables
  size_t num_variables_saved;
  infile.read((char*) &num_variables_saved, sizeof(num_variables_saved));

// Read treetype
  TreeType treetype;
  infile.read((char*) &treetype, sizeof(treetype));
  if (treetype != TREE_PROBABILITY) {
    throw std::runtime_error("Wrong treetype. Loaded file is not a probability estimation forest.");
  }

// Read class_values
  readVector1D(class_values, infile);

  for (size_t i = 0; i < num_trees; ++i) {

    // Read data
    std::vector<std::vector<size_t>> child_nodeIDs;
    readVector2D(child_nodeIDs, infile);
    std::vector<size_t> split_varIDs;
    readVector1D(split_varIDs, infile);
    std::vector<double> split_values;
    readVector1D(split_values, infile);

    // Read Terminal node class counts
    std::vector<size_t> terminal_nodes;
    readVector1D(terminal_nodes, infile);
    std::vector<std::vector<double>> terminal_class_counts_vector;
    readVector2D(terminal_class_counts_vector, infile);

    // Convert Terminal node class counts to vector with empty elemtents for non-terminal nodes
    std::vector<std::vector<double>> terminal_class_counts;
    terminal_class_counts.resize(child_nodeIDs[0].size(), std::vector<double>());
    for (size_t j = 0; j < terminal_nodes.size(); ++j) {
      terminal_class_counts[terminal_nodes[j]] = terminal_class_counts_vector[j];
    }

    // If dependent variable not in test data, change variable IDs accordingly
    if (num_variables_saved > num_variables) {
      for (auto& varID : split_varIDs) {
        if (varID >= dependent_varID) {
          --varID;
        }
      }
    }

    // Create tree
    trees.push_back(
        make_unique<TreeProbability>(child_nodeIDs, split_varIDs, split_values, &class_values, &response_classIDs,
            terminal_class_counts));
  }
}

const std::vector<double>& ForestProbability::getTreePrediction(size_t tree_idx, size_t sample_idx) const {
  const auto& tree = dynamic_cast<const TreeProbability&>(*trees[tree_idx]);
  return tree.getPrediction(sample_idx);
}

size_t ForestProbability::getTreePredictionTerminalNodeID(size_t tree_idx, size_t sample_idx) const {
  const auto& tree = dynamic_cast<const TreeProbability&>(*trees[tree_idx]);
  return tree.getPredictionTerminalNodeID(sample_idx);
}

// #nocov end

}// namespace ranger
