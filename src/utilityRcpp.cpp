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
    // Init result with NA
    for (size_t j = 0; j < groups.nrow(); ++j) {
      result(j, i) = NA_REAL;
    }
    
    // Order by terminal node ID
    std::vector<size_t> idx(groups.nrow());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(std::begin(idx), std::end(idx), [&](size_t j1, size_t j2) {return groups(j1, i) < groups(j2, i);});

    // Loop through change points (next node)
    size_t j = 0;
    while(j < idx.size()) {
      // Find next change point
      size_t k = j;
      while (k < idx.size() && groups(idx[j], i) == groups(idx[k], i)) {
        ++k;
      }
      
      // If other observation in same node
      if (k - j >= 2) {
        // Loop through observations between change points
        for (size_t l = j; l < k; ++l) {
          // Only OOB observations
          if (inbag_counts(idx[l], i) > 0) {
            continue;
          }
          
          // Select random observation in same terminal node, retry if same obs. selected
          size_t rnd = l;
          while (rnd == l) {
            rnd = j - 1 + Rcpp::sample(k - j, 1, false)[0];
          }
          result(idx[l], i) = y(idx[rnd]);
        }
      }

      // Next change point
      j = k;
    }
  }
  return result;
}

