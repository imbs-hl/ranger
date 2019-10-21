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

  Data(const Data&) = delete;
  Data& operator=(const Data&) = delete;

  virtual ~Data() = default;

  virtual double get_x(size_t row, size_t col) const = 0;
  virtual double get_y(size_t row, size_t col) const = 0;

  size_t getVariableID(const std::string& variable_name) const;

  virtual void reserveMemory(size_t y_cols) = 0;

  virtual void set_x(size_t col, size_t row, double value, bool& error) = 0;
  virtual void set_y(size_t col, size_t row, double value, bool& error) = 0;

  void addSnpData(unsigned char* snp_data, size_t num_cols_snp);

  bool loadFromFile(std::string filename, std::vector<std::string>& dependent_variable_names);
  bool loadFromFileWhitespace(std::ifstream& input_file, std::string header_line,
      std::vector<std::string>& dependent_variable_names);
  bool loadFromFileOther(std::ifstream& input_file, std::string header_line,
      std::vector<std::string>& dependent_variable_names, char seperator);

  void getAllValues(std::vector<double>& all_values, std::vector<size_t>& sampleIDs, size_t varID, size_t start,
      size_t end) const;

  void getMinMaxValues(double& min, double&max, std::vector<size_t>& sampleIDs, size_t varID, size_t start,
      size_t end) const;

  size_t getIndex(size_t row, size_t col) const {
    // Use permuted data for corrected impurity importance
    size_t col_permuted = col;
    if (col >= num_cols) {
      col = getUnpermutedVarID(col);
      row = getPermutedSampleID(row);
    }

    if (col < num_cols_no_snp) {
      return index_data[col * num_rows + row];
    } else {
      return getSnp(row, col, col_permuted);
    }
  }

  // #nocov start (cannot be tested anymore because GenABEL not on CRAN)
  size_t getSnp(size_t row, size_t col, size_t col_permuted) const {
    // Get data out of snp storage. -1 because of GenABEL coding.
    size_t idx = (col - num_cols_no_snp) * num_rows_rounded + row;
    size_t result = ((snp_data[idx / 4] & mask[idx % 4]) >> offset[idx % 4]) - 1;

    // TODO: Better way to treat missing values?
    if (result > 2) {
      result = 0;
    }

    // Order SNPs
    if (order_snps) {
      if (col_permuted >= num_cols) {
        result = snp_order[col_permuted - 2 * num_cols_no_snp][result];
      } else {
        result = snp_order[col - num_cols_no_snp][result];
      }
    }
    return result;
  }
  // #nocov end

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

  void orderSnpLevels(bool corrected_importance);

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

  std::vector<bool>& getIsOrderedVariable() noexcept {
    return is_ordered_variable;
  }

  void setIsOrderedVariable(const std::vector<std::string>& unordered_variable_names) {
    is_ordered_variable.resize(num_cols, true);
    for (auto& variable_name : unordered_variable_names) {
      size_t varID = getVariableID(variable_name);
      is_ordered_variable[varID] = false;
    }
  }

  void setIsOrderedVariable(std::vector<bool>& is_ordered_variable) {
    this->is_ordered_variable = is_ordered_variable;
  }

  bool isOrderedVariable(size_t varID) const {
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

  size_t getPermutedSampleID(size_t sampleID) const {
    return permuted_sampleIDs[sampleID];
  }

  size_t getUnpermutedVarID(size_t varID) const {
    if (varID >= num_cols) {
      varID -= num_cols;
    }
    return varID;
  }

  // #nocov start (cannot be tested anymore because GenABEL not on CRAN)
  const std::vector<std::vector<size_t>>& getSnpOrder() const {
    return snp_order;
  }

  void setSnpOrder(std::vector<std::vector<size_t>>& snp_order) {
    this->snp_order = snp_order;
    order_snps = true;
  }
  // #nocov end

protected:
  std::vector<std::string> variable_names;
  size_t num_rows;
  size_t num_rows_rounded;
  size_t num_cols;

  unsigned char* snp_data;
  size_t num_cols_no_snp;

  bool externalData;

  std::vector<size_t> index_data;
  std::vector<std::vector<double>> unique_data_values;
  size_t max_num_unique_values;

  // For each varID true if ordered
  std::vector<bool> is_ordered_variable;

  // Permuted samples for corrected impurity importance
  std::vector<size_t> permuted_sampleIDs;

  // Order of 0/1/2 for ordered splitting
  std::vector<std::vector<size_t>> snp_order;
  bool order_snps;
};

} // namespace ranger

#endif /* DATA_H_ */
