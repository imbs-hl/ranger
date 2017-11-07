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

##' Prediction with new data and a saved forest from Ranger.
##' 
##' For \code{type = 'response'} (the default), the predicted classes (classification), predicted numeric values (regression), predicted probabilities (probability estimation) or survival probabilities (survival) are returned. 
##' For \code{type = 'se'}, the standard error of the predictions are returned (regression only). The jackknife-after-bootstrap is used to estimate the standard errors based on out-of-bag predictions. See Wager et al. (2014) for details.
##' For \code{type = 'terminalNodes'}, the IDs of the terminal node in each tree for each observation in the given dataset are returned.
##' 
##' For classification and \code{predict.all = TRUE}, a factor levels are returned as numerics.
##' To retrieve the corresponding factor levels, use \code{rf$forest$levels}, if \code{rf} is the ranger object.
##'
##' @title Ranger prediction
##' @param object Ranger \code{ranger.forest} object.
##' @param data New test data of class \code{data.frame} or \code{gwaa.data} (GenABEL).
##' @param predict.all Return individual predictions for each tree instead of aggregated predictions for all trees. Return a matrix (sample x tree) for classification and regression, a 3d array for probability estimation (sample x class x tree) and survival (sample x time x tree).
##' @param num.trees Number of trees used for prediction. The first \code{num.trees} in the forest are used.
##' @param type Type of prediction. One of 'response', 'se', 'terminalNodes' with default 'response'. See below for details.
##' @param seed Random seed. Default is \code{NULL}, which generates the seed from \code{R}. Set to \code{0} to ignore the \code{R} seed. The seed is used in case of ties in classification mode.
##' @param num.threads Number of threads. Default is number of CPUs available.
##' @param verbose Verbose output on or off.
##' @param inbag.counts Number of times the observations are in-bag in the trees.
##' @param ... further arguments passed to or from other methods.
##' @return Object of class \code{ranger.prediction} with elements
##'   \tabular{ll}{
##'       \code{predictions}    \tab Predicted classes/values (only for classification and regression)  \cr
##'       \code{unique.death.times} \tab Unique death times (only for survival). \cr
##'       \code{chf} \tab Estimated cumulative hazard function for each sample (only for survival). \cr
##'       \code{survival} \tab Estimated survival function for each sample (only for survival). \cr
##'       \code{num.trees}   \tab Number of trees. \cr
##'       \code{num.independent.variables} \tab Number of independent variables. \cr
##'       \code{treetype}    \tab Type of forest/tree. Classification, regression or survival. \cr
##'       \code{num.samples}     \tab Number of samples.
##'   }
##' @references
##' \itemize{
##'   \item Wright, M. N. & Ziegler, A. (2017). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. J Stat Softw 77:1-17. \url{http://dx.doi.org/10.18637/jss.v077.i01}.
##'   \item Wager, S., Hastie T., & Efron, B. (2014). Confidence Intervals for Random Forests: The Jackknife and the Infinitesimal Jackknife. J Mach Learn Res 15:1625-1651. \url{http://jmlr.org/papers/v15/wager14a.html}.
##'   }
##' @seealso \code{\link{ranger}}
##' @author Marvin N. Wright
##' @importFrom Matrix Matrix
##' @export
predict.ranger.forest <- function(object, data, predict.all = FALSE,
                                  num.trees = object$num.trees, 
                                  type = "response",
                                  seed = NULL, num.threads = NULL,
                                  verbose = TRUE, inbag.counts = NULL,...) {

  ## GenABEL GWA data
  if ("gwaa.data" %in% class(data)) {
    snp.names <- snp.names(data)
    snp.data <- data@gtdata@gtps@.Data
    data <- data@phdata[, -1]
    gwa.mode <- TRUE
    variable.names <- c(names(data), snp.names)
  } else {
    snp.data <- as.matrix(0)
    gwa.mode <- FALSE
    variable.names <- colnames(data)
  }

  ## Check forest argument
  if (class(object) != "ranger.forest") {
    stop("Error: Invalid class of input object.")
  } else {
    forest <- object
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
  
  ## Check for old ranger version
  if (length(forest$child.nodeIDs) != forest$num.trees || length(forest$child.nodeIDs[[1]]) != 2) {
    stop("Error: Invalid forest object. Is the forest grown in ranger version <0.3.9? Try to predict with the same version the forest was grown.")
  }
  
  ## Prediction type
  if (type == "response" || type == "se") {
    prediction.type <- 1
  } else if (type == "terminalNodes") {
    prediction.type <- 2
  } else {
    stop("Error: Invalid value for 'type'. Use 'response' or 'terminalNodes'.")
  }
  
  ## Type "se" only for regression
  if (type == "se" && forest$treetype != "Regression") {
    stop("Error: Standard error prediction currently only available for regression.")
  }
  
  ## Type "se" requires keep.inbag=TRUE
  if (type == "se" && is.null(inbag.counts)) {
    stop("Error: No saved inbag counts in ranger object. Please set keep.inbag=TRUE when calling ranger.")
  }
  
  ## Set predict.all if type is "se"
  if (type == "se") {
    predict.all <- TRUE
  }

  ## Create final data
  if (forest$treetype == "Survival") {
    if (forest$dependent.varID > 0 && forest$status.varID > 1) {
      if (ncol(data) == length(forest$independent.variable.names)+2) {
        ## If alternative interface used and same data structure, don't subset data
        data.used <- data
      } else if (ncol(data) == length(forest$independent.variable.names)) {
        data.selected <- data[, forest$independent.variable.names, drop = FALSE]
        data.used <- cbind(0, 0, data.selected)
        variable.names <- c("time", "status", forest$independent.variable.names)
        forest$dependent.varID <- 0
        forest$status.varID <- 1
      } else {
        stop("Invalid prediction data. Include both time and status variable or none.")
      }
    } else {
      ## If formula interface used, subset data
      data.selected <- data[, forest$independent.variable.names, drop = FALSE]

      ## Arange data as in original data
      data.used <- cbind(0, 0, data.selected)
      variable.names <- c("time", "status", forest$independent.variable.names)
    }

  ## Index of no-recode variables
  idx.norecode <- c(-(forest$dependent.varID+1), -(forest$status.varID+1))

  } else {
    ## No survival
    if (ncol(data) == length(forest$independent.variable.names)+1 && forest$dependent.varID > 0) {
      ## If alternative interface used and same data structure, don't subset data
      data.used <- data
    } else {
      ## If formula interface used, subset data
      data.selected <- data[, forest$independent.variable.names, drop = FALSE]

      ## Arange data as in original data
      if (forest$dependent.varID == 0) {
        data.used <- cbind(0, data.selected)
        variable.names <- c("dependent", forest$independent.variable.names)
      } else if (forest$dependent.varID >= ncol(data)) {
        data.used <- cbind(data.selected, 0)
        variable.names <- c(forest$independent.variable.names, "dependent")
      } else {
        data.used <- cbind(data.selected[, 1:forest$dependent.varID],
                           0,
                           data.selected[, (forest$dependent.varID+1):ncol(data.selected)])
        variable.names <- c(forest$independent.variable.names[1:forest$dependent.varID],
                            "dependent",
                            forest$independent.variable.names[(forest$dependent.varID+1):length(forest$independent.variable.names)])
      }
    }

    ## Index of no-recode variables
    idx.norecode <- -(forest$dependent.varID+1)
  }

  ## Recode characters
  if (!is.matrix(data.used) && !inherits(data.used, "Matrix")) {
    char.columns <- sapply(data.used, is.character)
    data.used[char.columns] <- lapply(data.used[char.columns], factor)
  }

  ## Recode factors if forest grown 'order' mode
  if (!is.null(forest$covariate.levels) && !all(sapply(forest$covariate.levels, is.null))) {
    data.used[, idx.norecode] <- mapply(function(x, y) {
      if(is.null(y)) {
        x
      } else {
        new.levels <- setdiff(levels(x), y)
        factor(x, levels = c(y, new.levels))
      }
    }, data.used[, idx.norecode], forest$covariate.levels, SIMPLIFY = !is.data.frame(data.used[, idx.norecode]))
  }

  ## Convert to data matrix
  if (is.matrix(data.used) || inherits(data.used, "Matrix")) {
    data.final <- data.used
  } else {
    data.final <- data.matrix(data.used)
  }
  

  ## If gwa mode, add snp variable names
  if (gwa.mode) {
    variable.names <- c(variable.names, snp.names)
  }

  ## Check missing values
  if (any(is.na(data.final))) {
    offending_columns <- colnames(data.final)[colSums(is.na(data.final)) > 0]
    stop("Missing data in columns: ",
         paste0(offending_columns, collapse = ", "), ".", call. = FALSE)
  }

  if (sum(!(forest$independent.variable.names %in% variable.names)) > 0) {
    stop("Error: One or more independent variables not found in data.")
  }

  ## Num threads
  ## Default 0 -> detect from system in C++.
  if (is.null(num.threads)) {
    num.threads = 0
  } else if (!is.numeric(num.threads) || num.threads < 0) {
    stop("Error: Invalid value for num.threads")
  }

  ## Seed
  if (is.null(seed)) {
    seed <- runif(1 , 0, .Machine$integer.max)
  }

  if (forest$treetype == "Classification") {
    treetype <- 1
  } else if (forest$treetype == "Regression") {
    treetype <- 3
  } else if (forest$treetype == "Survival") {
    treetype <- 5
  } else if (forest$treetype == "Probability estimation") {
    treetype <- 9
  } else {
    stop("Error: Unknown tree type.")
  }

  ## Defaults for variables not needed
  dependent.variable.name <- "none"
  mtry <- 0
  importance <- 0
  min.node.size <- 0
  split.select.weights <- list(c(0, 0))
  use.split.select.weights <- FALSE
  always.split.variables <- c("0", "0")
  use.always.split.variables <- FALSE
  status.variable.name <- "status"
  prediction.mode <- TRUE
  write.forest <- FALSE
  replace <- TRUE
  probability <- FALSE
  unordered.factor.variables <- c("0", "0")
  use.unordered.factor.variables <- FALSE
  save.memory <- FALSE
  splitrule <- 1
  alpha <- 0
  minprop <- 0
  case.weights <- c(0, 0)
  use.case.weights <- FALSE
  keep.inbag <- FALSE
  sample.fraction <- 1
  holdout <- FALSE
  num.random.splits <- 1
  
  ## Use sparse matrix
  if ("dgCMatrix" %in% class(data.final)) {
    sparse.data <- data.final
    data.final <- matrix(c(0, 0))
    use.sparse.data <- TRUE
  } else {
    sparse.data <- Matrix(matrix(c(0, 0)))
    use.sparse.data <- FALSE
  }
  
  ## Call Ranger
  result <- rangerCpp(treetype, dependent.variable.name, data.final, variable.names, mtry,
                      num.trees, verbose, seed, num.threads, write.forest, importance,
                      min.node.size, split.select.weights, use.split.select.weights,
                      always.split.variables, use.always.split.variables,
                      status.variable.name, prediction.mode, forest, snp.data, replace, probability,
                      unordered.factor.variables, use.unordered.factor.variables, save.memory, splitrule,
                      case.weights, use.case.weights, predict.all, keep.inbag, sample.fraction,
                      alpha, minprop, holdout, prediction.type, num.random.splits, sparse.data, use.sparse.data)

  if (length(result) == 0) {
    stop("User interrupt or internal error.")
  }

  ## Prepare results
  result$num.samples <- nrow(data.final)
  result$treetype <- forest$treetype

  if (predict.all) {
    if (forest$treetype %in% c("Classification", "Regression")) {
      if (is.list(result$predictions)) {
        result$predictions <- do.call(rbind, result$predictions)
      } else {
        result$predictions <- array(result$predictions, dim = c(1, length(result$predictions)))
      }
    } else {
      ## TODO: Better solution for this?
      result$predictions <- aperm(array(unlist(result$predictions), 
                                        dim = rev(c(length(result$predictions), 
                                                    length(result$predictions[[1]]), 
                                                    length(result$predictions[[1]][[1]])))))
    }
  } else {
    if (is.list(result$predictions)) {
      result$predictions <- do.call(rbind, result$predictions)
    } 
  }
  
  if (type == "response") {
    if (forest$treetype == "Classification" && !is.null(forest$levels)) {
      if (!predict.all) {
        result$predictions <- integer.to.factor(result$predictions, forest$levels)
      }
    } else if (forest$treetype == "Regression") {
      ## Empty
    } else if (forest$treetype == "Survival") {
      result$unique.death.times <- forest$unique.death.times
      result$chf <- result$predictions
      result$predictions <- NULL
      result$survival <- exp(-result$chf)
    } else if (forest$treetype == "Probability estimation" && !is.null(forest$levels)) {
      if (!predict.all) {
        if (is.vector(result$predictions)) {
          result$predictions <- matrix(result$predictions, nrow = 1)
        }
        
        ## Set colnames and sort by levels
        colnames(result$predictions) <- forest$levels[forest$class.values]
        result$predictions <- result$predictions[, forest$levels, drop = FALSE]
      }
    }
  } 

  ## Compute Jackknife
  if (type == "se") {
    ## Aggregated predictions
    yhat <- rowMeans(result$predictions)

    ## Get inbag counts, keep only observations that are OOB at least once
    inbag.counts <- simplify2array(inbag.counts) 
    if (is.vector(inbag.counts)) {
      inbag.counts <- t(as.matrix(inbag.counts))
    }
    inbag.counts <- inbag.counts[rowSums(inbag.counts == 0) > 0, , drop = FALSE] 
    n <- nrow(inbag.counts)
    oob <- inbag.counts == 0
    
    if (all(!oob)) {
      stop("Error: No OOB observations found, consider increasing num.trees or reducing sample.fraction.")
    }

    ## Compute Jackknife
    jack.n <- sweep(tcrossprod(result$predictions, oob), 
                    2, rowSums(oob), "/", check.margin = FALSE)
    if (is.vector(jack.n)) {
      jack.n <- t(as.matrix(jack.n))
    }
    jack <- (n - 1) / n * rowSums((jack.n - yhat)^2)
    bias <- (exp(1) - 1) * n / result$num.trees^2 * rowSums((result$predictions - yhat)^2)
    jab <- pmax(jack - bias, 0)
    result$se <- sqrt(jab)
    
    ## Response as predictions
    result$predictions <- yhat
  }

  class(result) <- "ranger.prediction"
  return(result)
}

##' Prediction with new data and a saved forest from Ranger.
##' 
##' For \code{type = 'response'} (the default), the predicted classes (classification), predicted numeric values (regression), predicted probabilities (probability estimation) or survival probabilities (survival) are returned. 
##' For \code{type = 'se'}, the standard error of the predictions are returned (regression only). The jackknife-after-bootstrap is used to estimate the standard errors based on out-of-bag predictions. See Wager et al. (2014) for details.
##' For \code{type = 'terminalNodes'}, the IDs of the terminal node in each tree for each observation in the given dataset are returned.
##' 
##' For classification and \code{predict.all = TRUE}, a factor levels are returned as numerics.
##' To retrieve the corresponding factor levels, use \code{rf$forest$levels}, if \code{rf} is the ranger object.
##'
##' @title Ranger prediction
##' @param object Ranger \code{ranger} object.
##' @param data New test data of class \code{data.frame} or \code{gwaa.data} (GenABEL).
##' @param predict.all Return individual predictions for each tree instead of aggregated predictions for all trees. Return a matrix (sample x tree) for classification and regression, a 3d array for probability estimation (sample x class x tree) and survival (sample x time x tree).
##' @param num.trees Number of trees used for prediction. The first \code{num.trees} in the forest are used.
##' @param type Type of prediction. One of 'response', 'se', 'terminalNodes' with default 'response'. See below for details.
##' @param seed Random seed. Default is \code{NULL}, which generates the seed from \code{R}. Set to \code{0} to ignore the \code{R} seed. The seed is used in case of ties in classification mode.
##' @param num.threads Number of threads. Default is number of CPUs available.
##' @param verbose Verbose output on or off.
##' @param ... further arguments passed to or from other methods.
##' @return Object of class \code{ranger.prediction} with elements
##'   \tabular{ll}{
##'       \code{predictions}    \tab Predicted classes/values (only for classification and regression)  \cr
##'       \code{unique.death.times} \tab Unique death times (only for survival). \cr
##'       \code{chf} \tab Estimated cumulative hazard function for each sample (only for survival). \cr
##'       \code{survival} \tab Estimated survival function for each sample (only for survival). \cr
##'       \code{num.trees}   \tab Number of trees. \cr
##'       \code{num.independent.variables} \tab Number of independent variables. \cr
##'       \code{treetype}    \tab Type of forest/tree. Classification, regression or survival. \cr
##'       \code{num.samples}     \tab Number of samples.
##'   }
##' @references
##' \itemize{
##'   \item Wright, M. N. & Ziegler, A. (2017). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. J Stat Softw 77:1-17. \url{http://dx.doi.org/10.18637/jss.v077.i01}.
##'   \item Wager, S., Hastie T., & Efron, B. (2014). Confidence Intervals for Random Forests: The Jackknife and the Infinitesimal Jackknife. J Mach Learn Res 15:1625-1651. \url{http://jmlr.org/papers/v15/wager14a.html}.
##'   }
##' @seealso \code{\link{ranger}}
##' @author Marvin N. Wright
##' @export
predict.ranger <- function(object, data, predict.all = FALSE,
                           num.trees = object$num.trees,
                           type = "response",
                           seed = NULL, num.threads = NULL,
                           verbose = TRUE, ...) {
  forest <- object$forest
  if (is.null(forest)) {
    stop("Error: No saved forest in ranger object. Please set write.forest to TRUE when calling ranger.")
  }
  if (object$importance.mode %in% c("impurity_corrected", "impurity_unbiased")) {
    warning("Forest was grown with 'impurity_corrected' variable importance. For prediction it is STRONGLY advised to grow another forest without this importance setting.")
  }
  
  predict(forest, data, predict.all, num.trees, type, seed, num.threads, verbose, object$inbag.counts, ...)
}
