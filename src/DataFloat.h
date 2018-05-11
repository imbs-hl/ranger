/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

// Ignore in coverage report (not used in R package)
// #nocov start
#ifndef DATAFLOAT_H_
#define DATAFLOAT_H_

#include <vector>

#include "globals.h"
#include "Data.h"

namespace ranger {

class DataFloat: public Data {
public:
  DataFloat() = default;
  DataFloat(double* data_double, std::vector<std::string> variable_names, size_t num_rows, size_t num_cols);

  DataFloat(const DataFloat&) = delete;
  DataFloat& operator=(const DataFloat&) = delete;

  virtual ~DataFloat() override = default;

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
    data[col * num_rows + row] = value;
  }

private:
  std::vector<float> data;
};

} // namespace ranger

#endif /* DATAFLOAT_H_ */
// #nocov end

