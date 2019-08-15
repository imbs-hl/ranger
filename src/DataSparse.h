/*-------------------------------------------------------------------------------
 This file is part of Ranger.

 Ranger is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Ranger is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with Ranger. If not, see <http://www.gnu.org/licenses/>.

 Written by:

 Marvin N. Wright
 Institut für Medizinische Biometrie und Statistik
 Universität zu Lübeck
 Ratzeburger Allee 160
 23562 Lübeck

 http://www.imbs-luebeck.de
 #-------------------------------------------------------------------------------*/

#ifndef DATASPARSE_H_
#define DATASPARSE_H_

#include <RcppEigen.h>

#include "globals.h"
#include "utility.h"
#include "Data.h"

namespace ranger {

class DataSparse: public Data {
public:
  DataSparse() = default;

  DataSparse(Eigen::SparseMatrix<double>& data, std::vector<std::string> variable_names, size_t num_rows,
      size_t num_cols);

  DataSparse(const DataSparse&) = delete;
  DataSparse& operator=(const DataSparse&) = delete;

  virtual ~DataSparse() override = default;

  double get(size_t row, size_t col) const override {
    // Use permuted data for corrected impurity importance
    size_t col_permuted = col;
    if (col >= num_cols) {
      col = getUnpermutedVarID(col);
      row = getPermutedSampleID(row);
    }
    return data.coeff(row, col);
  }

  void reserveMemory() override {
    data.resize(num_rows, num_cols);
  }

  void set(size_t col, size_t row, double value, bool& error) override {
    data.coeffRef(row, col) = value;
  }

private:
  Eigen::SparseMatrix<double> data;
};

} // namespace ranger

#endif /* DATASPARSE_H_ */
