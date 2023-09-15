# -------------------------------------------------------------------------------
#   This file is part of Ranger.
#
# Ranger is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Ranger is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Ranger. If not, see <http://www.gnu.org/licenses/>.
#
# Written by:
#
#   Marvin N. Wright
# Institut fuer Medizinische Biometrie und Statistik
# Universitaet zu Luebeck
# Ratzeburger Allee 160
# 23562 Luebeck
# Germany
#
# http://www.imbs-luebeck.de
# -------------------------------------------------------------------------------


#' Horizontal shrinkage
#' 
#' Apply horizontal shrinkage to a ranger object. 
#' Horizontal shrinkage is a regularization technique that recursively shrinks node predictions towards parent node predictions. 
#' For details see Agarwal et al. (2022).
#'
#' @param rf ranger object, created with \code{node.stats = TRUE}. 
#' @param lambda Non-negative shrinkage parameter. 
#'
#' @return The ranger object is modified in-place. 
#'
#' @examples
##' @references
##' \itemize{
##'   \item Agarwal, A., Tan, Y.S., Ronen, O., Singh, C. & Yu, B. (2022). Hierarchical Shrinkage: Improving the accuracy and interpretability of tree-based models. Proceedings of the 39th International Conference on Machine Learning, PMLR 162:111-135.
##'   }
#' @author Marvin N. Wright
#' @export
hshrink <- function(rf, lambda) {
  if (is.null(rf$forest$num.samples.nodes)) {
    stop("Horizontal shrinkage needs node statistics, set node.stats=TRUE in ranger() call.")
  }
  if (lambda < 0) {
    stop("Shrinkage parameter lambda has to be non-negative.")
  }
  
  if (rf$treetype == "Regression") {
    invisible(lapply(1:rf$num.trees, function(treeID) {
      hshrink_regr(
        rf$forest$child.nodeIDs[[treeID]][[1]], rf$forest$child.nodeIDs[[treeID]][[2]], 
        rf$forest$num.samples.nodes[[treeID]], rf$forest$node.predictions[[treeID]], 
        rf$forest$split.values[[treeID]], lambda, 0, 0, 0, 0
      )
    }))
  } else if (rf$treetype == "Probability estimation") {
    invisible(lapply(1:rf$num.trees, function(treeID) {
      # Create temporary class frequency matrix
      class_freq <- t(simplify2array(rf$forest$terminal.class.counts[[treeID]]))
      
      parent_pred <- rep(0, length(rf$forest$class.values))
      cum_sum <- rep(0, length(rf$forest$class.values))
      hshrink_prob(
        rf$forest$child.nodeIDs[[treeID]][[1]], rf$forest$child.nodeIDs[[treeID]][[2]], 
        rf$forest$num.samples.nodes[[treeID]], class_freq, 
        lambda, 0, 0, parent_pred, cum_sum 
      )
      
      # Assign temporary matrix values back to ranger object
      replace_class_counts(rf$forest$terminal.class.counts[[treeID]], class_freq)
    }))
  } else if (rf$treetype == "Classification") {
    stop("To apply horizontal shrinkage to classification forests, use probability=TRUE in the ranger() call.")
  } else if (rf$treetype == "Survival") {
    stop("Horizontal shrinkage not yet implemented for survival.")
  } else {
    stop("Unknown treetype.")
  }
  
}


