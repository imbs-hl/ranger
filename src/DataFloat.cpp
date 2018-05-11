/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

// Ignore in coverage report (not used in R package)
// #nocov start
#include <iostream>

#include "DataFloat.h"

namespace ranger {

DataFloat::DataFloat(double* data_double, std::vector<std::string> variable_names, size_t num_rows, size_t num_cols) {
  this->variable_names = variable_names;
  this->num_rows = num_rows;
  this->num_cols = num_cols;
  this->num_cols_no_snp = num_cols;

  reserveMemory();
  for (size_t i = 0; i < num_cols; ++i) {
    for (size_t j = 0; j < num_rows; ++j) {
      data[i * num_rows + j] = (float) data_double[i * num_rows + j];
    }
  }
}

// #nocov end

}// namespace ranger
