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

  for (size_t i = 0; i < groups.ncol(); ++i) {
    // Loop through observations
    std::vector<int> nodes;
    std::vector<int> idx;
    for (size_t j = 0; j < groups.nrow(); ++j) {
      result(j, i) = NA_REAL;
      nodes.push_back(groups(j, i));
      idx.push_back(j);
    }
    // sort the relevant array by nodes
    std::sort(idx.begin(), idx.end(),
               [&](int n1, int n2){ return nodes[n1] < nodes[n2]; });

    for (int j = 0; j < idx.size();) {
      int k = j;
      while (k != idx.size() && nodes[idx[j]] == nodes[idx[k]]) ++k;
      for (int l = j; l < k; ++l) {
        if (inbag_counts(idx[l], i) > 0) continue;
        if (k-j > 1) {
          int rnd = l;
          while (rnd == l) rnd = j - 1 + Rcpp::sample(k-j, 1, false)[0];
          result(idx[l], i) = y(idx[rnd]);
        }
      }
      j = k;
    }
  }
  return result;
}

