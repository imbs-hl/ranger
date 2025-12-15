/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <cstring>

#include "Data.h"
#include "utility.h"

namespace ranger {

Data::Data() :
    num_rows(0), num_rows_rounded(0), num_cols(0), snp_data(0), num_cols_no_snp(0), externalData(true), index_data(0), max_num_unique_values(
        0), order_snps(false), any_na(false) {
}

size_t Data::getVariableID(const std::string& variable_name) const {
  auto it = std::find(variable_names.cbegin(), variable_names.cend(), variable_name);
  if (it == variable_names.cend()) {
    throw std::runtime_error("Variable " + variable_name + " not found.");
  }
  return (std::distance(variable_names.cbegin(), it));
}

// #nocov start (cannot be tested anymore because GenABEL not on CRAN)
void Data::addSnpData(unsigned char* snp_data, size_t num_cols_snp) {
  num_cols = num_cols_no_snp + num_cols_snp;
  num_rows_rounded = roundToNextMultiple(num_rows, 4);
  this->snp_data = snp_data;
}
// #nocov end

// #nocov start
bool Data::loadFromFile(std::string filename, std::vector<std::string>& dependent_variable_names) {

  bool result;

  // Open input file
  std::ifstream input_file;
  input_file.open(filename);
  if (!input_file.good()) {
    throw std::runtime_error("Could not open input file.");
  }

  // Count number of rows
  size_t line_count = 0;
  std::string line;
  while (getline(input_file, line)) {
    ++line_count;
  }
  num_rows = line_count - 1;
  input_file.close();
  input_file.open(filename);

  // Check if comma, semicolon or whitespace separated
  std::string header_line;
  getline(input_file, header_line);

  // Find out if comma, semicolon or whitespace separated and call appropriate method
  if (header_line.find(',') != std::string::npos) {
    result = loadFromFileOther(input_file, header_line, dependent_variable_names, ',');
  } else if (header_line.find(';') != std::string::npos) {
    result = loadFromFileOther(input_file, header_line, dependent_variable_names, ';');
  } else {
    result = loadFromFileWhitespace(input_file, header_line, dependent_variable_names);
  }

  externalData = false;
  input_file.close();
  return result;
}

bool Data::loadFromFileWhitespace(std::ifstream& input_file, std::string header_line,
    std::vector<std::string>& dependent_variable_names) {

  size_t num_dependent_variables = dependent_variable_names.size();
  std::vector<size_t> dependent_varIDs;
  dependent_varIDs.resize(num_dependent_variables);

  // Read header
  std::string header_token;
  std::stringstream header_line_stream(header_line);
  size_t col = 0;
  while (header_line_stream >> header_token) {
    bool is_dependent_var = false;
    for (size_t i = 0; i < dependent_variable_names.size(); ++i) {
      if (header_token == dependent_variable_names[i]) {
        dependent_varIDs[i] = col;
        is_dependent_var = true;
      }
    }
    if (!is_dependent_var) {
      variable_names.push_back(header_token);
    }
    ++col;
  }

  num_cols = variable_names.size();
  num_cols_no_snp = num_cols;

  // Read body
  reserveMemory(num_dependent_variables);
  bool error = false;
  std::string line;
  size_t row = 0;
  while (getline(input_file, line)) {
    double token;
    std::stringstream line_stream(line);
    size_t column = 0;
    while (readFromStream(line_stream, token)) {
      size_t column_x = column;
      bool is_dependent_var = false;
      for (size_t i = 0; i < dependent_varIDs.size(); ++i) {
        if (column == dependent_varIDs[i]) {
          set_y(i, row, token, error);
          is_dependent_var = true;
          break;
        } else if (column > dependent_varIDs[i]) {
          --column_x;
        }
      }
      if (!is_dependent_var) {
        set_x(column_x, row, token, error);
      }
      ++column;
    }
    if (column > (num_cols + num_dependent_variables)) {
      throw std::runtime_error(
          std::string("Could not open input file. Too many columns in row ") + std::to_string(row) + std::string("."));
    } else if (column < (num_cols + num_dependent_variables)) {
      throw std::runtime_error(
          std::string("Could not open input file. Too few columns in row ") + std::to_string(row)
              + std::string(". Are all values numeric?"));
    }
    ++row;
  }
  num_rows = row;
  return error;
}

bool Data::loadFromFileOther(std::ifstream& input_file, std::string header_line,
    std::vector<std::string>& dependent_variable_names, char separator) {

  size_t num_dependent_variables = dependent_variable_names.size();
  std::vector<size_t> dependent_varIDs;
  dependent_varIDs.resize(num_dependent_variables);

  // Read header
  std::string header_token;
  std::stringstream header_line_stream(header_line);
  size_t col = 0;
  while (getline(header_line_stream, header_token, separator)) {
    bool is_dependent_var = false;
    for (size_t i = 0; i < dependent_variable_names.size(); ++i) {
      if (header_token == dependent_variable_names[i]) {
        dependent_varIDs[i] = col;
        is_dependent_var = true;
      }
    }
    if (!is_dependent_var) {
      variable_names.push_back(header_token);
    }
    ++col;
  }

  num_cols = variable_names.size();
  num_cols_no_snp = num_cols;

  // Read body
  reserveMemory(num_dependent_variables);
  bool error = false;
  std::string line;
  size_t row = 0;
  while (getline(input_file, line)) {
    std::string token_string;
    double token;
    std::stringstream line_stream(line);
    size_t column = 0;
    while (getline(line_stream, token_string, separator)) {
      std::stringstream token_stream(token_string);
      readFromStream(token_stream, token);

      size_t column_x = column;
      bool is_dependent_var = false;
      for (size_t i = 0; i < dependent_varIDs.size(); ++i) {
        if (column == dependent_varIDs[i]) {
          set_y(i, row, token, error);
          is_dependent_var = true;
          break;
        } else if (column > dependent_varIDs[i]) {
          --column_x;
        }
      }
      if (!is_dependent_var) {
        set_x(column_x, row, token, error);
      }
      ++column;
    }
    ++row;
  }
  num_rows = row;
  return error;
}

void Data::loadSnpsFromFilePlink(std::ifstream& bed_file, std::ifstream& fam_file, std::ifstream& bim_file) {
  if (!bed_file) {
    throw std::runtime_error("Cannot open .bed file.");
  }
  
  uint8_t header[3];
  bed_file.read(reinterpret_cast<char*>(header), 3);
  
  if (header[0] != 0x6C || header[1] != 0x1B || header[2] != 0x01) {
    throw std::runtime_error("Invalid or unsupported .bed file");
  }

  // Get dimensions
  size_t n_samples = count_fam_samples(fam_file);
  size_t n_snps = count_bim_snps(bim_file);

  if (n_samples == 0 || n_snps == 0) {
    throw std::runtime_error("Empty .fam or .bim file.");
  }
  // if ((bed_file.tellg() - 3) != n_snps * ((n_samples + 3) / 4)) {
  //   throw std::runtime_error("BED/FAM/BIM dimension mismatch.");
  // }
  if (n_samples != num_rows) {
    throw std::runtime_error("Geno/Pheno sample size mismatch.");
  }

  const size_t bytes_per_snp = (n_samples + 3) / 4;
  std::vector<uint8_t> buffer(bytes_per_snp);

  // Reserve memory
  const size_t total_genotypes = n_samples * n_snps;
  const size_t total_bytes = (total_genotypes + 3) / 4;
  snp_data = new unsigned char[total_bytes];
  std::memset(snp_data, 0, total_bytes);

  for (size_t snp = 0; snp < n_snps; ++snp) {
      bed_file.read(reinterpret_cast<char*>(buffer.data()), bytes_per_snp);

      for (size_t i = 0; i < n_samples; ++i) {
          // Decode PLINK
          uint8_t byte = buffer[i >> 2];
          uint8_t plink_bits = (byte >> ((i & 3) << 1)) & 0x03;

          uint8_t genabel_bits = plink2genabel[plink_bits];

          // Store in GenABEL format
          size_t idx = i * n_snps + snp; 

          size_t byte_idx = idx >> 2;
          size_t bit_pos  = 6 - ((idx & 3) << 1);

          snp_data[byte_idx] |= (genabel_bits << bit_pos);
      }
  }
  num_cols = num_cols_no_snp + n_snps;
  num_rows_rounded = roundToNextMultiple(num_rows, 4);
  owns_snp_data = true;
}
// #nocov end

void Data::getAllValues(std::vector<double>& all_values, std::vector<size_t>& sampleIDs, size_t varID, size_t start,
    size_t end) const {

  // All values for varID (no duplicates) for given sampleIDs
  if (getUnpermutedVarID(varID) < num_cols_no_snp) {
    
    all_values.reserve(end - start);
    for (size_t pos = start; pos < end; ++pos) {
      all_values.push_back(get_x(sampleIDs[pos], varID));
    }
    if (any_na) {
      std::sort(all_values.begin(), all_values.end(), less_nan<double>);
    } else {
      std::sort(all_values.begin(), all_values.end());
    }
    all_values.erase(std::unique(all_values.begin(), all_values.end()), all_values.end());
    
    // Keep only one NaN value
    if (any_na) {
      while (all_values.size() >= 2 && std::isnan(all_values[all_values.size() - 2])) {
        all_values.pop_back();
      }
    }
  } else {
    // If GWA data just use 0, 1, 2
    all_values = std::vector<double>( { 0, 1, 2 });
  }
}

void Data::getMinMaxValues(double& min, double&max, std::vector<size_t>& sampleIDs, size_t varID, size_t start,
    size_t end) const {
  if (sampleIDs.size() > 0) {
    min = get_x(sampleIDs[start], varID);
    max = min;
  }
  for (size_t pos = start; pos < end; ++pos) {
    double value = get_x(sampleIDs[pos], varID);
    if (value < min) {
      min = value;
    }
    if (value > max) {
      max = value;
    }
  }
}

void Data::sort() {

  // Reserve memory
  index_data.resize(num_cols_no_snp * num_rows);

  // For all columns, get unique values and save index for each observation
  for (size_t col = 0; col < num_cols_no_snp; ++col) {

    // Get all unique values
    std::vector<double> unique_values(num_rows);
    for (size_t row = 0; row < num_rows; ++row) {
      unique_values[row] = get_x(row, col);
    }
    
    if (any_na) {
      std::sort(unique_values.begin(), unique_values.end(), less_nan<double>);
    } else {
      std::sort(unique_values.begin(), unique_values.end());
    }
    unique_values.erase(unique(unique_values.begin(), unique_values.end()), unique_values.end());

    // Get index of unique value
    for (size_t row = 0; row < num_rows; ++row) {
      size_t idx;
      if (any_na) {
        idx = std::lower_bound(unique_values.begin(), unique_values.end(), get_x(row, col), less_nan<double>) - unique_values.begin();
      } else {
        idx = std::lower_bound(unique_values.begin(), unique_values.end(), get_x(row, col)) - unique_values.begin();
      }
      index_data[col * num_rows + row] = idx;
    }
    
    // Save unique values (keep NaN)
    if (any_na) {
      while (unique_values.size() >= 2 && std::isnan(unique_values[unique_values.size() - 2])) {
        unique_values.pop_back();
      }
    }
    unique_data_values.push_back(unique_values);
    if (unique_values.size() > max_num_unique_values) {
      max_num_unique_values = unique_values.size();
    }
  }
}

// TODO: Implement ordering for multiclass and survival
// #nocov start (cannot be tested anymore because GenABEL not on CRAN)
void Data::orderSnpLevels(bool corrected_importance) {
  // Stop if now SNP data
  if (snp_data == 0) {
    return;
  }

  size_t num_snps;
  if (corrected_importance) {
    num_snps = 2 * (num_cols - num_cols_no_snp);
  } else {
    num_snps = num_cols - num_cols_no_snp;
  }

  // Reserve space
  snp_order.resize(num_snps, std::vector<size_t>(3));

  // For each SNP
  for (size_t i = 0; i < num_snps; ++i) {
    size_t col = i;
    if (i >= (num_cols - num_cols_no_snp)) {
      // Get unpermuted SNP ID
      col = i - num_cols + num_cols_no_snp;
    }

    // Order by mean response
    std::vector<double> means(3, 0);
    std::vector<double> counts(3, 0);
    for (size_t row = 0; row < num_rows; ++row) {
      size_t row_permuted = row;
      if (i >= (num_cols - num_cols_no_snp)) {
        row_permuted = getPermutedSampleID(row);
      }
      size_t idx = col * num_rows_rounded + row_permuted;
      size_t value = (((snp_data[idx / 4] & mask[idx % 4]) >> offset[idx % 4]) - 1);

      // TODO: Better way to treat missing values?
      if (value > 2) {
        value = 0;
      }

      means[value] += get_y(row, 0);
      ++counts[value];
    }

    for (size_t value = 0; value < 3; ++value) {
      means[value] /= counts[value];
    }

    // Save order
    snp_order[i] = order(means, false);
  }

  order_snps = true;
}
// #nocov end

} // namespace ranger

