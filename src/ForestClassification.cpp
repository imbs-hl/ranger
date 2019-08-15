/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#include <unordered_map>
#include <algorithm>
#include <iterator>
#include <random>
#include <stdexcept>
#include <cmath>
#include <string>

#include "utility.h"
#include "ForestClassification.h"
#include "TreeClassification.h"
#include "Data.h"

namespace ranger {

void ForestClassification::loadForest(size_t dependent_varID, size_t num_trees,
    std::vector<std::vector<std::vector<size_t>> >& forest_child_nodeIDs,
    std::vector<std::vector<size_t>>& forest_split_varIDs, std::vector<std::vector<double>>& forest_split_values,
    std::vector<double>& class_values, std::vector<bool>& is_ordered_variable) {

  this->dependent_varID = dependent_varID;
  this->num_trees = num_trees;
  this->class_values = class_values;
  data->setIsOrderedVariable(is_ordered_variable);

  // Create trees
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    trees.push_back(
        make_unique<TreeClassification>(forest_child_nodeIDs[i], forest_split_varIDs[i], forest_split_values[i],
            &this->class_values, &response_classIDs));
  }

  // Create thread ranges
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);
}

void ForestClassification::initInternal(std::string status_variable_name) {

  // If mtry not set, use floored square root of number of independent variables.
  if (mtry == 0) {
    unsigned long temp = sqrt((double) (num_variables - 1));
    mtry = std::max((unsigned long) 1, temp);
  }

  // Set minimal node size
  if (min_node_size == 0) {
    min_node_size = DEFAULT_MIN_NODE_SIZE_CLASSIFICATION;
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

void ForestClassification::growInternal() {
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    trees.push_back(
        make_unique<TreeClassification>(&class_values, &response_classIDs, &sampleIDs_per_class, &class_weights));
  }
}

void ForestClassification::allocatePredictMemory() {
  size_t num_prediction_samples = data->getNumRows();
  if (predict_all || prediction_type == TERMINALNODES) {
    predictions = std::vector<std::vector<std::vector<double>>>(1,
        std::vector<std::vector<double>>(num_prediction_samples, std::vector<double>(num_trees)));
  } else {
    predictions = std::vector<std::vector<std::vector<double>>>(1,
        std::vector<std::vector<double>>(1, std::vector<double>(num_prediction_samples)));
  }
}

void ForestClassification::predictInternal(size_t sample_idx) {
  if (predict_all || prediction_type == TERMINALNODES) {
    // Get all tree predictions
    for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
      if (prediction_type == TERMINALNODES) {
        predictions[0][sample_idx][tree_idx] = getTreePredictionTerminalNodeID(tree_idx, sample_idx);
      } else {
        predictions[0][sample_idx][tree_idx] = getTreePrediction(tree_idx, sample_idx);
      }
    }
  } else {
    // Count classes over trees and save class with maximum count
    std::unordered_map<double, size_t> class_count;
    for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
      ++class_count[getTreePrediction(tree_idx, sample_idx)];
    }
    predictions[0][0][sample_idx] = mostFrequentValue(class_count, random_number_generator);
  }
}

void ForestClassification::computePredictionErrorInternal() {

  // Class counts for samples
  std::vector<std::unordered_map<double, size_t>> class_counts;
  class_counts.reserve(num_samples);
  for (size_t i = 0; i < num_samples; ++i) {
    class_counts.push_back(std::unordered_map<double, size_t>());
  }

  // For each tree loop over OOB samples and count classes
  for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
    for (size_t sample_idx = 0; sample_idx < trees[tree_idx]->getNumSamplesOob(); ++sample_idx) {
      size_t sampleID = trees[tree_idx]->getOobSampleIDs()[sample_idx];
      ++class_counts[sampleID][getTreePrediction(tree_idx, sample_idx)];
    }
  }

  // Compute majority vote for each sample
  predictions = std::vector<std::vector<std::vector<double>>>(1,
      std::vector<std::vector<double>>(1, std::vector<double>(num_samples)));
  for (size_t i = 0; i < num_samples; ++i) {
    if (!class_counts[i].empty()) {
      predictions[0][0][i] = mostFrequentValue(class_counts[i], random_number_generator);
    } else {
      predictions[0][0][i] = NAN;
    }
  }

  // Compare predictions with true data
  size_t num_missclassifications = 0;
  size_t num_predictions = 0;
  for (size_t i = 0; i < predictions[0][0].size(); ++i) {
    double predicted_value = predictions[0][0][i];
    if (!std::isnan(predicted_value)) {
      ++num_predictions;
      double real_value = data->get(i, dependent_varID);
      if (predicted_value != real_value) {
        ++num_missclassifications;
      }
      ++classification_table[std::make_pair(real_value, predicted_value)];
    }
  }
  overall_prediction_error = (double) num_missclassifications / (double) num_predictions;
}

// #nocov start
void ForestClassification::writeOutputInternal() {
  if (verbose_out) {
    *verbose_out << "Tree type:                         " << "Classification" << std::endl;
  }
}

void ForestClassification::writeConfusionFile() {

  // Open confusion file for writing
  std::string filename = output_prefix + ".confusion";
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to confusion file: " + filename + ".");
  }

  // Write confusion to file
  outfile << "Overall OOB prediction error (Fraction missclassified): " << overall_prediction_error << std::endl;
  outfile << std::endl;
  outfile << "Class specific prediction errors:" << std::endl;
  outfile << "           ";
  for (auto& class_value : class_values) {
    outfile << "     " << class_value;
  }
  outfile << std::endl;
  for (auto& predicted_value : class_values) {
    outfile << "predicted " << predicted_value << "     ";
    for (auto& real_value : class_values) {
      size_t value = classification_table[std::make_pair(real_value, predicted_value)];
      outfile << value;
      if (value < 10) {
        outfile << "     ";
      } else if (value < 100) {
        outfile << "    ";
      } else if (value < 1000) {
        outfile << "   ";
      } else if (value < 10000) {
        outfile << "  ";
      } else if (value < 100000) {
        outfile << " ";
      }
    }
    outfile << std::endl;
  }

  outfile.close();
  if (verbose_out)
    *verbose_out << "Saved confusion matrix to file " << filename << "." << std::endl;
}

void ForestClassification::writePredictionFile() {

  // Open prediction file for writing
  std::string filename = output_prefix + ".prediction";
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to prediction file: " + filename + ".");
  }

  // Write
  outfile << "Predictions: " << std::endl;
  if (predict_all) {
    for (size_t k = 0; k < num_trees; ++k) {
      outfile << "Tree " << k << ":" << std::endl;
      for (size_t i = 0; i < predictions.size(); ++i) {
        for (size_t j = 0; j < predictions[i].size(); ++j) {
          outfile << predictions[i][j][k] << std::endl;
        }
      }
      outfile << std::endl;
    }
  } else {
    for (size_t i = 0; i < predictions.size(); ++i) {
      for (size_t j = 0; j < predictions[i].size(); ++j) {
        for (size_t k = 0; k < predictions[i][j].size(); ++k) {
          outfile << predictions[i][j][k] << std::endl;
        }
      }
    }
  }

  if (verbose_out)
    *verbose_out << "Saved predictions to file " << filename << "." << std::endl;
}

void ForestClassification::saveToFileInternal(std::ofstream& outfile) {

  // Write num_variables
  outfile.write((char*) &num_variables, sizeof(num_variables));

  // Write treetype
  TreeType treetype = TREE_CLASSIFICATION;
  outfile.write((char*) &treetype, sizeof(treetype));

  // Write class_values
  saveVector1D(class_values, outfile);
}

void ForestClassification::loadFromFileInternal(std::ifstream& infile) {

  // Read number of variables
  size_t num_variables_saved;
  infile.read((char*) &num_variables_saved, sizeof(num_variables_saved));

  // Read treetype
  TreeType treetype;
  infile.read((char*) &treetype, sizeof(treetype));
  if (treetype != TREE_CLASSIFICATION) {
    throw std::runtime_error("Wrong treetype. Loaded file is not a classification forest.");
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
        make_unique<TreeClassification>(child_nodeIDs, split_varIDs, split_values, &class_values, &response_classIDs));
  }
}

double ForestClassification::getTreePrediction(size_t tree_idx, size_t sample_idx) const {
  const auto& tree = dynamic_cast<const TreeClassification&>(*trees[tree_idx]);
  return tree.getPrediction(sample_idx);
}

size_t ForestClassification::getTreePredictionTerminalNodeID(size_t tree_idx, size_t sample_idx) const {
  const auto& tree = dynamic_cast<const TreeClassification&>(*trees[tree_idx]);
  return tree.getPredictionTerminalNodeID(sample_idx);
}

// #nocov end

}// namespace ranger
