/*-------------------------------------------------------------------------------
This file is part of ranger.

Copyright (c) [2014-2018] [Marvin N. Wright]

This software may be modified and distributed under the terms of the MIT license.

Please note that the C++ core of ranger is distributed under MIT license and the
R package "ranger" under GPL3 license.
#-------------------------------------------------------------------------------*/

#ifndef DATA_H_
#define DATA_H_

#include <vector>
#include <iostream>
#include <numeric>
#include <random>
#include <algorithm>

#include "globals.h"

namespace ranger {

class Data {
public:
  Data();
  virtual ~Data();

  virtual double get(size_t row, size_t col) const = 0;

  size_t getVariableID(std::string variable_name);

  virtual void reserveMemory() = 0;
  virtual void set(size_t col, size_t row, double value, bool& error) = 0;

  void addSnpData(unsigned char* snp_data, size_t num_cols_snp);

  bool loadFromFile(std::string filename);
  bool loadFromFileWhitespace(std::ifstream& input_file, std::string header_line);
  bool loadFromFileOther(std::ifstream& input_file, std::string header_line, char seperator);

  void getAllValues(std::vector<double>& all_values, std::vector<size_t>& sampleIDs, size_t varID);

  void getMinMaxValues(double& min, double&max, std::vector<size_t>& sampleIDs, size_t varID);

  size_t getIndex(size_t row, size_t col) const {
    // Use permuted data for corrected impurity importance
    if (col >= num_cols) {
      col = getUnpermutedVarID(col);
      row = getPermutedSampleID(row);
    }

    if (col < num_cols_no_snp) {
      return index_data[col * num_rows + row];
    } else {
      // Get data out of snp storage. -1 because of GenABEL coding.
      size_t idx = (col - num_cols_no_snp) * num_rows_rounded + row;
      size_t result = (((snp_data[idx / 4] & mask[idx % 4]) >> offset[idx % 4]) - 1);

      // TODO: Better way to treat missing values?
      if (result > 2) {
        return 0;
      } else {
        return result;
      }
    }
  }

  double getUniqueDataValue(size_t varID, size_t index) const {
    // Use permuted data for corrected impurity importance
    if (varID >= num_cols) {
      varID = getUnpermutedVarID(varID);
    }

    if (varID < num_cols_no_snp) {
      return unique_data_values[varID][index];
    } else {
      // For GWAS data the index is the value
      return (index);
    }
  }

  size_t getNumUniqueDataValues(size_t varID) const {
    // Use permuted data for corrected impurity importance
    if (varID >= num_cols) {
      varID = getUnpermutedVarID(varID);
    }

    if (varID < num_cols_no_snp) {
      return unique_data_values[varID].size();
    } else {
      // For GWAS data 0,1,2
      return (3);
    }
  }

  void sort();

  const std::vector<std::string>& getVariableNames() const {
    return variable_names;
  }
  size_t getNumCols() const {
    return num_cols;
  }
  size_t getNumRows() const {
    return num_rows;
  }

  size_t getMaxNumUniqueValues() const {
    if (snp_data == 0 || max_num_unique_values > 3) {
      // If no snp data or one variable with more than 3 unique values, return that value
      return max_num_unique_values;
    } else {
      // If snp data and no variable with more than 3 unique values, return 3
      return 3;
    }
  }

  std::vector<size_t>& getNoSplitVariables() {
    return no_split_variables;
  }

  void addNoSplitVariable(size_t varID) {
    no_split_variables.push_back(varID);
    std::sort(no_split_variables.begin(), no_split_variables.end());
  }

  std::vector<bool>& getIsOrderedVariable() {
    return is_ordered_variable;
  }

  void setIsOrderedVariable(std::vector<std::string>& unordered_variable_names) {
    is_ordered_variable.resize(num_cols, true);
    for (auto& variable_name : unordered_variable_names) {
      size_t varID = getVariableID(variable_name);
      is_ordered_variable[varID] = false;
    }
  }

  void setIsOrderedVariable(std::vector<bool>& is_ordered_variable) {
    this->is_ordered_variable = is_ordered_variable;
  }

  const bool isOrderedVariable(size_t varID) const {
    // Use permuted data for corrected impurity importance
    if (varID >= num_cols) {
      varID = getUnpermutedVarID(varID);
    }
    return is_ordered_variable[varID];
  }

  void permuteSampleIDs(std::mt19937_64 random_number_generator) {
    permuted_sampleIDs.resize(num_rows);
    std::iota(permuted_sampleIDs.begin(), permuted_sampleIDs.end(), 0);
    std::shuffle(permuted_sampleIDs.begin(), permuted_sampleIDs.end(), random_number_generator);
  }

  const size_t getPermutedSampleID(size_t sampleID) const {
    return permuted_sampleIDs[sampleID];
  }

  const size_t getUnpermutedVarID(size_t varID) const {
    if (varID >= num_cols) {
      varID -= num_cols;

      for (auto& skip : no_split_variables) {
        if (varID >= skip) {
          ++varID;
        }
      }
    }
    return varID;
  }

protected:
  std::vector<std::string> variable_names;
  size_t num_rows;
  size_t num_rows_rounded;
  size_t num_cols;

  unsigned char* snp_data;
  size_t num_cols_no_snp;

  bool externalData;

  size_t* index_data;
  std::vector<std::vector<double>> unique_data_values;
  size_t max_num_unique_values;

  // Variable to not split at (only dependent_varID for non-survival trees)
  std::vector<size_t> no_split_variables;

  // For each varID true if ordered
  std::vector<bool> is_ordered_variable;

  // Permuted samples for corrected impurity importance
  std::vector<size_t> permuted_sampleIDs;

private:
  DISALLOW_COPY_AND_ASSIGN(Data);
};

} // namespace ranger

#endif /* DATA_H_ */
