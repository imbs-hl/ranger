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


#' Tree information in human readable format
#' 
#' Extract tree information of a \code{ranger} object. 
#' 
#' Node and variable ID's are 0-indexed, i.e., node 0 is the root node. 
#' If the formula interface is used in the \code{ranger} call, the variable ID's are usually different to the original data used to grow the tree. 
#' Refer to the variable name instead to be sure.
#' 
#' @param object \code{ranger} object.
#' @param tree Number of the tree of interest.
#' @return A data.frame with the columns
#' \tabular{ll}{
#'       \code{A} \tab A. \cr
#'       \code{B} \tab B. \cr
#'       \code{C} \tab C. \cr
#'   }
#' @examples
#' require(ranger)
#' rf <- ranger(Species ~ ., data = iris)
#' treeInfo(rf, 1)
#' @seealso \code{\link{ranger}}
#' @author Marvin N. Wright
#' @export
treeInfo <- function(object, tree = 1) {
  if (class(object) != "ranger" & class(object) != "holdoutRF") {
    stop("Error: Invalid class of input object.")
  } 
  forest <- object$forest
  if (is.null(forest)) {
    stop("Error: No saved forest in ranger object. Please set write.forest to TRUE when calling ranger.")
  }
  if (is.null(forest$dependent.varID) || is.null(forest$num.trees) ||
      is.null(forest$child.nodeIDs) || is.null(forest$split.varIDs) ||
      is.null(forest$split.values) || is.null(forest$independent.variable.names) ||
      is.null(forest$treetype)) {
    stop("Error: Invalid forest object.")
  }
  if (forest$treetype == "Survival" && (is.null(forest$status.varID)  ||
                                        is.null(forest$chf) || is.null(forest$unique.death.times))) {
    stop("Error: Invalid forest object.")
  }
  if (length(forest$child.nodeIDs) != forest$num.trees || length(forest$child.nodeIDs[[1]]) != 2) {
    stop("Error: Invalid forest object. Is the forest grown in ranger version <0.3.9? Try with the same version the forest was grown.")
  }
  if (tree > forest$num.trees) {
    stop("Error: Requesting tree ", tree, ", but forest has only ", forest$num.trees, " trees.")
  }
  
  result <- data.frame(nodeID = 0:(length(forest$split.values[[tree]]) - 1),
                       left = forest$child.nodeIDs[[tree]][[1]], 
                       right = forest$child.nodeIDs[[tree]][[2]], 
                       splitvarID = forest$split.varIDs[[tree]], 
                       splitvarName = "X",
                       splitval = forest$split.values[[tree]], 
                       prediction = forest$split.values[[tree]], 
                       terminal = FALSE)
  
  result$left[result$left == 0] <- NA
  result$right[result$right == 0] <- NA
  result$terminal[is.na(result$left)] <- TRUE
  result$splitvarID[result$terminal] <- NA
  result$splitvarName[result$terminal] <- NA
  result$splitval[result$terminal] <- NA
  result$prediction[!result$terminal] <- NA
  
  result
  # TODO: Tests
  # TODO: Variable names
  # TODO: Prediction for factors
  # TODO: Probability?
  # TODO: Survival?
  # TODO: Differences formula/alternative interface?
  # TODO: Unordered splitting (different methods)?
}
