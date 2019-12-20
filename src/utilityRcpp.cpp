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
Germany

http://www.imbs-luebeck.de
#-------------------------------------------------------------------------------*/

#include <Rcpp.h>

// Count number of elements in reference smaller than values
//[[Rcpp::export]]
Rcpp::IntegerVector numSmaller(Rcpp::NumericVector values, Rcpp::NumericVector reference) {
  std::sort(reference.begin(), reference.end());
  Rcpp::IntegerVector result(values.size());
  for (int i = 0; i < values.size(); ++i)
    result[i] = std::lower_bound(reference.begin(), reference.end(), values[i]) - reference.begin();
  return result;
}

// Get random other obs. in same terminal node
//[[Rcpp::export]]
Rcpp::NumericMatrix randomObsNode(Rcpp::IntegerMatrix groups, Rcpp::NumericVector y, Rcpp::IntegerMatrix inbag_counts) {
  Rcpp::NumericMatrix result(groups.nrow(), groups.ncol());
  
  // Loop through trees
  for (size_t i = 0; i < groups.ncol(); ++i) {
    // Loop through observations
    for (size_t j = 0; j < groups.nrow(); ++j) {
      result(j, i) = NA_REAL;
      
      if (inbag_counts(j, i) > 0) {
        continue;
      }
      
      // Search for other observations with same group
      Rcpp::IntegerVector others;
      for (size_t k = 0; k < groups.nrow(); ++k) {
        if (j != k) {
          if (groups(j, i) == groups(k, i)) {
            others.push_back(k);
          }
        }
      }
      
      // Randomly select one
      if (others.size() > 0) {
        result(j, i) = y(Rcpp::sample(others, 1, false)[0]);
      }
    }
  }
  return result;
}

