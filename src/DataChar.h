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

#include <limits.h>

#include "globals.h"
#include "Data.h"

class DataChar: public Data {
public:
  DataChar();
  DataChar(double* data_double, std::vector<std::string> variable_names, size_t num_rows, size_t num_cols, bool& error);
  virtual ~DataChar();

  double get(size_t row, size_t col) const {
    // Use permuted data for corrected impurity importance
    if (col >= num_cols) {
      col = getUnpermutedVarID(col);
      row = getPermutedSampleID(row);
    }

    if (col < num_cols_no_snp) {
      return data[col * num_rows + row];
    } else {
      // Get data out of snp storage. -1 because of GenABEL coding.
      size_t idx = (col - num_cols_no_snp) * num_rows_rounded + row;
      size_t result = ((snp_data[idx / 4] & mask[idx % 4]) >> offset[idx % 4]) - 1;

      // TODO: Better way to treat missing values?
      if (result > 2) {
        result = 0;
      }

      // Order SNPs
      if (snp_order.empty()) {
        return result;
      } else {
        return snp_order[col - num_cols_no_snp][result];
      }
    }
  }

  void reserveMemory() {
    data = new char[num_cols * num_rows];
  }

  void set(size_t col, size_t row, double value, bool& error) {
    if (value > CHAR_MAX || value < CHAR_MIN) {
      error = true;
    }
    if (floor(value) != ceil(value)) {
      error = true;
    }
    data[col * num_rows + row] = (char) value;
  }

private:
  char* data;

  DISALLOW_COPY_AND_ASSIGN(DataChar);
};

#endif /* DATACHAR_H_ */
// #nocov end
