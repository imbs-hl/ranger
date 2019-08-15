/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

// Ignore in coverage report (not used in R package)
// #nocov start
#ifndef DATACHAR_H_
#define DATACHAR_H_

#include <vector>
#include <limits.h>

#include "globals.h"
#include "Data.h"

namespace ranger {

class DataChar: public Data {
public:
  DataChar() = default;
  DataChar(double* data_double, std::vector<std::string> variable_names, size_t num_rows, size_t num_cols, bool& error);

  DataChar(const DataChar&) = delete;
  DataChar& operator=(const DataChar&) = delete;

  virtual ~DataChar() override = default;

  double get(size_t row, size_t col) const override {
    // Use permuted data for corrected impurity importance
    size_t col_permuted = col;
    if (col >= num_cols) {
      col = getUnpermutedVarID(col);
      row = getPermutedSampleID(row);
    }

    if (col < num_cols_no_snp) {
      return data[col * num_rows + row];
    } else {
      return getSnp(row, col, col_permuted);
    }
  }

  void reserveMemory() override {
    data.resize(num_cols * num_rows);
  }

  void set(size_t col, size_t row, double value, bool& error) override {
    if (value > CHAR_MAX || value < CHAR_MIN) {
      error = true;
    }
    if (floor(value) != ceil(value)) {
      error = true;
    }
    data[col * num_rows + row] = value;
  }

private:
  std::vector<char> data;
};

} // namespace ranger

#endif /* DATACHAR_H_ */
// #nocov end
