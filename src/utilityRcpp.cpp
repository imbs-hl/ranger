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

// Recursive function for hierarchical shrinkage (regression)
//[[Rcpp::export]]
void hshrink_regr(Rcpp::IntegerVector& left_children, Rcpp::IntegerVector& right_children, 
                  Rcpp::IntegerVector& num_samples_nodes, Rcpp::NumericVector& node_predictions, 
                  Rcpp::NumericVector& split_values, double lambda,
                  size_t nodeID, size_t parent_n, double parent_pred, double cum_sum) {
  if (nodeID == 0) {
    // In the root, just use the prediction
    cum_sum = node_predictions[nodeID];
  } else {
    // If not root, use shrinkage formula
    cum_sum += (node_predictions[nodeID] - parent_pred) / (1 + lambda/parent_n);
  }
  
  if (left_children[nodeID] == 0) {
    // If leaf, change node prediction in split_values (used for prediction) 
    split_values[nodeID] = cum_sum;
  } else {
    // If not leaf, give weighted prediction to child nodes
    hshrink_regr(left_children, right_children, num_samples_nodes, node_predictions, split_values, 
                    lambda, left_children[nodeID], num_samples_nodes[nodeID], node_predictions[nodeID],
                    cum_sum);
    hshrink_regr(left_children, right_children, num_samples_nodes, node_predictions, split_values, 
                    lambda, right_children[nodeID], num_samples_nodes[nodeID], node_predictions[nodeID],
                    cum_sum);
  }
}

// Recursive function for hierarchical shrinkage (probability)
//[[Rcpp::export]]
void hshrink_prob(Rcpp::IntegerVector& left_children, Rcpp::IntegerVector& right_children, 
                  Rcpp::IntegerVector& num_samples_nodes, 
                  Rcpp::NumericMatrix& class_freq, double lambda,
                  size_t nodeID, size_t parent_n, Rcpp::NumericVector parent_pred, Rcpp::NumericVector cum_sum) {
  
  if (nodeID == 0) {
    // In the root, just use the prediction
    cum_sum = class_freq(nodeID, Rcpp::_);
  } else {
    // If not root, use shrinkage formula
    cum_sum += (class_freq(nodeID, Rcpp::_) - parent_pred) / (1 + lambda/parent_n);
  }
  
  if (left_children[nodeID] == 0) {
    // If leaf, change node prediction in split_values (used for prediction) 
    class_freq(nodeID, Rcpp::_) = cum_sum;
  } else {
    // If not leaf, give weighted prediction to child nodes
    hshrink_prob(left_children, right_children, num_samples_nodes, class_freq, lambda, 
                 left_children[nodeID], num_samples_nodes[nodeID], class_freq(nodeID, Rcpp::_), clone(cum_sum));
    hshrink_prob(left_children, right_children, num_samples_nodes, class_freq, lambda, 
                 right_children[nodeID], num_samples_nodes[nodeID], class_freq(nodeID, Rcpp::_), clone(cum_sum));
  }
}

// Replace class counts list(vector) with values from matrix
//[[Rcpp::export]]
void replace_class_counts(Rcpp::List& class_counts_old, Rcpp::NumericMatrix& class_counts_new) {
  for (size_t i = 0; i < class_counts_old.size(); ++i) {
    class_counts_old[i] = class_counts_new(i, Rcpp::_);
  }
}

