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
#' Splitting at unordered factors (nominal variables) depends on the option \code{respect.unordered.factors} in the \code{ranger} call. 
#' For the "ignore" and "order" approaches, all values smaller or equal the \code{splitval} value go to the left and all values larger go to the right, as usual. 
#' However, with "order" the values correspond to the order in \code{object$forest$covariate.levels} instead of the original order (usually alphabetical).
#' In the "partition" mode, the \code{splitval} values for unordered factor are comma separated lists of values, representing the factor levels (in the original order) going to the right.
#' 
#' @param object \code{ranger} object.
#' @param tree Number of the tree of interest.
#' @return A data.frame with the columns
#' \tabular{ll}{
#'       \code{nodeID} \tab The nodeID, 0-indexed. \cr
#'       \code{leftChild} \tab ID of the left child node, 0-indexed. \cr
#'       \code{rightChild} \tab ID of the right child node, 0-indexed. \cr
#'       \code{splitvarID} \tab ID of the splitting variable, 0-indexed. Caution, the variable order changes if the formula interface is used. \cr
#'       \code{splitvarName} \tab Name of the splitting variable. \cr
#'       \code{splitval} \tab The splitting value. For numeric or ordinal variables, all values smaller or equal go to the left, larger values to the right. For unordered factor variables see above. \cr
#'       \code{terminal} \tab Logical, TRUE for terminal nodes. \cr
#'       \code{prediction} \tab One column with the predicted class (factor) for classification and the predicted numerical value for regression. One probability per class for probability estimation in several columns. Nothing for survival, refer to \code{object$forest$chf} for the CHF node predictions. \cr
#'       \code{numSamples} \tab Number of samples in the node (only if ranger called with \code{node.stats = TRUE}). \cr
#'       \code{splitStat} \tab Split statistics, i.e., value of the splitting criterion (only if ranger called with \code{node.stats = TRUE}). \cr
#'   }
#' @examples
#' rf <- ranger(Species ~ ., data = iris)
#' treeInfo(rf, 1)
#' @seealso \code{\link{ranger}}
#' @author Marvin N. Wright
#' @export
treeInfo <- function(object, tree = 1) {
  if (!inherits(object, "ranger")) {
    stop("Error: Invalid class of input object.")
  }
  forest <- object$forest
  if (is.null(forest)) {
    stop("Error: No saved forest in ranger object. Please set write.forest to TRUE when calling ranger.")
  }
  if (is.null(forest$num.trees) ||
      is.null(forest$child.nodeIDs) || is.null(forest$split.varIDs) ||
      is.null(forest$split.values) || is.null(forest$independent.variable.names) ||
      is.null(forest$treetype)) {
    stop("Error: Invalid forest object.")
  }
  if (forest$treetype == "Survival" && (is.null(forest$chf) || is.null(forest$unique.death.times))) {
    stop("Error: Invalid forest object.")
  }
  if (length(forest$child.nodeIDs) != forest$num.trees || length(forest$child.nodeIDs[[1]]) < 2 || length(forest$child.nodeIDs[[1]]) > 3) {
    stop("Error: Invalid forest object. Is the forest grown in ranger version <0.3.9? Try with the same version the forest was grown.")
  }
  if (!is.null(forest$dependent.varID)) {
    forest <- convert.pre.xy(forest, trees = tree)
  }
  if (tree > forest$num.trees) {
    stop("Error: Requesting tree ", tree, ", but forest has only ", forest$num.trees, " trees.")
  }
  
  result <- data.frame(nodeID = 0:(length(forest$split.values[[tree]]) - 1),
                       leftChild = forest$child.nodeIDs[[tree]][[1]], 
                       rightChild = forest$child.nodeIDs[[tree]][[2]], 
                       splitvarID = forest$split.varIDs[[tree]], 
                       splitvarName = "X",
                       splitval = forest$split.values[[tree]], 
                       terminal = FALSE)
  
  result$leftChild[result$leftChild == 0] <- NA
  result$rightChild[result$rightChild == 0] <- NA
  result$terminal[is.na(result$leftChild)] <- TRUE
  result$splitvarID[result$terminal] <- NA
  result$splitvarName[result$terminal] <- NA
  result$splitval[result$terminal] <- NA
  result$splitvarName <- forest$independent.variable.names[result$splitvarID + 1]

  ## Unordered splitting
  idx.unordered <- !result$terminal & !forest$is.ordered[result$splitvarID + 1]
  if (any(idx.unordered)) {
    if (any(result$splitval[idx.unordered] > (2^31 - 1))) {
      warning("Unordered splitting levels can only be shown for up to 31 levels.")
      result$splitval[idx.unordered] <- NA
    } else {
      result$splitval[idx.unordered] <- sapply(result$splitval[idx.unordered], function(x) {
        paste(which(as.logical(intToBits(x))), collapse = ",")
      })
    }
  }
  
  ## Prediction
  if (forest$treetype == "Classification") {
    if (is.null(forest$num.samples.nodes)) {
      # split.stats=FALSE
      result$prediction <- forest$split.values[[tree]]
      result$prediction[!result$terminal] <- NA
    } else {
      # split.stats=TRUE
      result$prediction <- forest$node.predictions[[tree]]
    }
    if (!is.null(forest$levels)) {
      result$prediction <- integer.to.factor(result$prediction, labels = forest$levels)
    }
  } else if (forest$treetype == "Regression") {
    if (is.null(forest$num.samples.nodes)) {
      # split.stats=FALSE
      result$prediction <- forest$split.values[[tree]]
      result$prediction[!result$terminal] <- NA
    } else {
      # split.stats=TRUE
      result$prediction <- forest$node.predictions[[tree]]
    }
  } else if (forest$treetype == "Probability estimation") {
    predictions <- matrix(nrow = nrow(result), ncol = length(forest$class.values))
    if (is.null(forest$num.samples.nodes)) {
      # split.stats=FALSE
      predictions[result$terminal, ] <- do.call(rbind, forest$terminal.class.counts[[tree]])
    } else {
      # split.stats=TRUE
      predictions <- do.call(rbind, forest$terminal.class.counts[[tree]])
    }
    if (!is.null(forest$levels)) {
      colnames(predictions) <- forest$levels[forest$class.values]
      predictions <- predictions[, forest$levels[sort(forest$class.values)], drop = FALSE]
    } else {
      colnames(predictions) <- forest$class.values
    }
    colnames(predictions) <- paste0("pred.", colnames(predictions))
    result <- data.frame(result, predictions)
  } else if (forest$treetype == "Survival") {
    # No prediction for survival (CHF too large?)
  } else {
    stop("Error: Unknown tree type.")
  }
  
  ## Node statistics
  if (!is.null(forest$num.samples.nodes)) {
    result$numSamples <- forest$num.samples.nodes[[tree]]
  }
  if (!is.null(forest$split.stats)) {
    result$splitStat <- forest$split.stats[[tree]]
    result$splitStat[result$terminal] <- NA
  }
  
  result
}
