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
##' For \code{type = 'se'}, the standard error of the predictions are returned (regression only). The jackknife-after-bootstrap or infinitesimal jackknife for bagging is used to estimate the standard errors based on out-of-bag predictions. See Wager et al. (2014) for details.
##' For \code{type = 'terminalNodes'}, the IDs of the terminal node in each tree for each observation in the given dataset are returned.
##' 
##' If \code{type = 'se'} is selected, the method to estimate the variances can be chosen with \code{se.method}. Set \code{se.method = 'jack'} for jackknife after bootstrap and \code{se.method = 'infjack'} for the infinitesimal jackknife for bagging.
##' 
##' For classification and \code{predict.all = TRUE}, a factor levels are returned as numerics.
##' To retrieve the corresponding factor levels, use \code{rf$forest$levels}, if \code{rf} is the ranger object.
##'
##' @title Ranger prediction
##' @param object Ranger \code{ranger.forest} object.
##' @param data New test data of class \code{data.frame} or \code{gwaa.data} (GenABEL). 
##' @param predict.all Return individual predictions for each tree instead of aggregated predictions for all trees. Return a matrix (sample x tree) for classification and regression, a 3d array for probability estimation (sample x class x tree) and survival (sample x time x tree).
##' @param num.trees Number of trees used for prediction. The first \code{num.trees} in the forest are used.
##' @param type Type of prediction. One of 'response', 'se', 'terminalNodes', 'quantiles' with default 'response'. See below for details.
##' @param se.method Method to compute standard errors. One of 'jack', 'infjack' with default 'infjack'. Only applicable if type = 'se'. See below for details.
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
##'   \item Wright, M. N. & Ziegler, A. (2017). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. J Stat Softw 77:1-17. \doi{10.18637/jss.v077.i01}.
##'   \item Wager, S., Hastie T., & Efron, B. (2014). Confidence Intervals for Random Forests: The Jackknife and the Infinitesimal Jackknife. J Mach Learn Res 15:1625-1651. \url{https://jmlr.org/papers/v15/wager14a.html}.
##'   }
##' @seealso \code{\link{ranger}}
##' @author Marvin N. Wright
##' @importFrom Matrix Matrix
##' @export
predict.ranger.forest <- function(object, data, predict.all = FALSE,
                                  num.trees = object$num.trees, 
                                  type = "response", se.method = "infjack",
                                  seed = NULL, num.threads = NULL,
                                  verbose = TRUE, inbag.counts = NULL, ...) {

  ## GenABEL GWA data
  if (inherits(data, "gwaa.data")) {
    snp.names <- snp.names(data)
    snp.data <- data@gtdata@gtps@.Data
    data <- data@phdata[, -1, drop = FALSE]
    gwa.mode <- TRUE
  } else {
    snp.data <- as.matrix(0)
    gwa.mode <- FALSE
  }

  ## Check forest argument
  if (!inherits(object, "ranger.forest")) {
    stop("Error: Invalid class of input object.")
  } else {
    forest <- object
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
  
  ## Check for old ranger version
  if (length(forest$child.nodeIDs) != forest$num.trees || length(forest$child.nodeIDs[[1]]) != 2) {
    stop("Error: Invalid forest object. Is the forest grown in ranger version <0.3.9? Try to predict with the same version the forest was grown.")
  }
  if (!is.null(forest$dependent.varID)) {
    warning("Forest grown in ranger version <0.11.5, converting ...")
    forest <- convert.pre.xy(forest)
  }
  
  ## Prediction type
  if (type == "response" || type == "se") {
    prediction.type <- 1
  } else if (type == "terminalNodes") {
    prediction.type <- 2
  } else if (type == "quantiles") {
    stop("Error: Apply predict() to the ranger object instead of the $forest object to predict quantiles.")
  } else {
    stop("Error: Invalid value for 'type'. Use 'response', 'se', 'terminalNodes', or 'quantiles'.")
  }
  
  ## Type "se" only for certain tree types
  if (type == "se" && se.method == "jack" && forest$treetype != "Regression") {
    stop("Error: Jackknife standard error prediction currently only available for regression.")
  }
  if (type == "se" && se.method == "infjack") {
   if (forest$treetype == "Survival") {
     stop("Error: Infinitesimal jackknife standard error prediction not yet available for survival.")
   } else if (forest$treetype == "Classification") {
     stop("Error: Not a probability forest. Set probability=TRUE to use the infinitesimal jackknife standard error prediction for classification.")
   }
  }
  
  ## Type "se" requires keep.inbag=TRUE
  if (type == "se" && is.null(inbag.counts)) {
    stop("Error: No saved inbag counts in ranger object. Please set keep.inbag=TRUE when calling ranger.")
  }
  
  ## Set predict.all if type is "se"
  if (type == "se") {
    predict.all <- TRUE
  }
  
  x <- data
  
  if (sum(!(forest$independent.variable.names %in% colnames(x))) > 0) {
    stop("Error: One or more independent variables not found in data.")
  }

  ## Subset to same column as in training if necessary
  if (length(colnames(x)) != length(forest$independent.variable.names) || any(colnames(x) != forest$independent.variable.names)) {
    x <- x[, forest$independent.variable.names, drop = FALSE]
  }

  ## Recode characters
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    char.columns <- sapply(x, is.character)
    if (length(char.columns) > 0) {
      x[char.columns] <- lapply(x[char.columns], factor)
    }
  }
  
  ## Recode factors if forest grown 'order' mode
  if (!is.null(forest$covariate.levels) && !all(sapply(forest$covariate.levels, is.null)) && !is.matrix(x)) {
    x <- mapply(function(xx, yy) {
      if(is.null(yy)) {
        xx
      } else {
        new.levels <- setdiff(levels(xx), yy)
        factor(xx, levels = c(yy, new.levels), exclude = NULL)
      }
    }, x, forest$covariate.levels, SIMPLIFY = !is.data.frame(x))
  }
  if (is.list(x) && !is.data.frame(x)) {
    x <- as.data.frame(x)
  }

  ## Convert to data matrix
  if (!is.matrix(x) & !inherits(x, "Matrix")) {
    x <- data.matrix(x)
  }

  ## Check missing values
  if (any(is.na(x))) {
    offending_columns <- colnames(x)[colSums(is.na(x)) > 0]
    stop("Missing data in columns: ",
         paste0(offending_columns, collapse = ", "), ".", call. = FALSE)
  }

  ## Num threads
  ## Default 0 -> detect from system in C++.
  if (is.null(num.threads)) {
    num.threads <- as.integer(Sys.getenv("R_RANGER_NUM_THREADS", getOption("ranger.num.threads", getOption("Ncpus", 2L))))
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
  mtry <- 0
  importance <- 0
  min.node.size <- 0
  min.bucket <- 0
  split.select.weights <- list(c(0, 0))
  use.split.select.weights <- FALSE
  always.split.variables <- c("0", "0")
  use.always.split.variables <- FALSE
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
  class.weights <- c(0, 0)
  keep.inbag <- FALSE
  sample.fraction <- 1
  holdout <- FALSE
  num.random.splits <- 1
  order.snps <- FALSE
  oob.error <- FALSE
  max.depth <- 0
  inbag <- list(c(0,0))
  use.inbag <- FALSE
  y <- matrix(c(0, 0))
  regularization.factor <- c(0, 0)
  use.regularization.factor <- FALSE
  regularization.usedepth <- FALSE
  node.stats <- FALSE
  time.interest <- c(0, 0)
  use.time.interest <- FALSE
  
  ## Use sparse matrix
  if (inherits(x, "dgCMatrix")) {
    sparse.x <- x
    x <- matrix(c(0, 0))
    use.sparse.data <- TRUE
  } else {
    sparse.x <- Matrix(matrix(c(0, 0)))
    use.sparse.data <- FALSE
    x <- data.matrix(x)
  }
  
  ## Call Ranger
  result <- rangerCpp(treetype, x, y, forest$independent.variable.names, mtry,
                      num.trees, verbose, seed, num.threads, write.forest, importance,
                      min.node.size, min.bucket, split.select.weights, use.split.select.weights,
                      always.split.variables, use.always.split.variables,
                      prediction.mode, forest, snp.data, replace, probability,
                      unordered.factor.variables, use.unordered.factor.variables, save.memory, splitrule,
                      case.weights, use.case.weights, class.weights, 
                      predict.all, keep.inbag, sample.fraction, alpha, minprop, holdout, 
                      prediction.type, num.random.splits, sparse.x, use.sparse.data,
                      order.snps, oob.error, max.depth, inbag, use.inbag, 
                      regularization.factor, use.regularization.factor, regularization.usedepth, 
                      node.stats, time.interest, use.time.interest)

  if (length(result) == 0) {
    stop("User interrupt or internal error.")
  }

  ## Prepare results
  result$num.samples <- nrow(x)
  result$treetype <- forest$treetype

  if (predict.all) {
    if (forest$treetype %in% c("Classification", "Regression")) {
      if (is.list(result$predictions)) {
        result$predictions <- do.call(rbind, result$predictions)
      } else {
        result$predictions <- array(result$predictions, dim = c(1, length(result$predictions)))
      }
    } else {
      if (is.list(result$predictions) & length(result$predictions) >= 1 & is.numeric(result$predictions[[1]])) {
        # Fix for single test observation
        result$predictions <- list(result$predictions)
      }
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
    } else if (forest$treetype == "Probability estimation") {
      if (predict.all) {
        ## Set colnames and sort by levels
        if (!is.null(forest$levels)) {
          colnames(result$predictions) <- forest$levels[forest$class.values]
          result$predictions <- result$predictions[, forest$levels[sort(forest$class.values)], , drop = FALSE]
        }
      } else {
        if (is.vector(result$predictions)) {
          result$predictions <- matrix(result$predictions, nrow = 1)
        }
        
        ## Set colnames and sort by levels
        if (!is.null(forest$levels)) {
          colnames(result$predictions) <- forest$levels[forest$class.values]
          result$predictions <- result$predictions[, forest$levels[sort(forest$class.values)], drop = FALSE]
        }
      }
    }
  } else if (type == "terminalNodes") {
    if (is.vector(result$predictions)) {
      result$predictions <- matrix(result$predictions, nrow = 1)
    }
  }

  ## Compute Jackknife
  if (type == "se") {
    ## Aggregated predictions
    if (length(dim(result$predictions)) > 2) {
      yhat <- apply(result$predictions, c(1, 2), mean)
    } else {
      yhat <- rowMeans(result$predictions)
    }

    ## Get inbag counts, keep only observations that are OOB at least once
    inbag.counts <- simplify2array(inbag.counts) 
    if (is.vector(inbag.counts)) {
      inbag.counts <- t(as.matrix(inbag.counts))
    }
    inbag.counts <- inbag.counts[rowSums(inbag.counts == 0) > 0, , drop = FALSE] 
    n <- nrow(inbag.counts)
    oob <- inbag.counts == 0
    if (num.trees != object$num.trees) {
      oob <- oob[, 1:num.trees]
    }
    
    if (all(!oob)) {
      stop("Error: No OOB observations found, consider increasing num.trees or reducing sample.fraction.")
    }

    if (se.method == "jack") {
      ## Compute Jackknife
      oob.count <- rowSums(oob)
      jack.n <- sweep(tcrossprod(result$predictions, oob), 
                      2, oob.count, "/", check.margin = FALSE)
      if (is.vector(jack.n)) {
        jack.n <- t(as.matrix(jack.n))
      }
      if (any(oob.count == 0)) {
        n <- sum(oob.count > 0)
        jack.n <- jack.n[, oob.count > 0]
      } 
      jack <- (n - 1) / n * rowSums((jack.n - yhat)^2)
      bias <- (exp(1) - 1) * n / result$num.trees^2 * rowSums((result$predictions - yhat)^2)
      jab <- pmax(jack - bias, 0)
      result$se <- sqrt(jab)
    } else if (se.method == "infjack") {
      if (forest$treetype == "Regression") {
        infjack <- rInfJack(pred = result$predictions, inbag = inbag.counts, used.trees = 1:num.trees)
        result$se <- sqrt(infjack$var.hat)
      } else if (forest$treetype == "Probability estimation") {
        infjack <- apply(result$predictions, 2, function(x) {
          rInfJack(x, inbag.counts)$var.hat
        })
        result$se <- sqrt(infjack)
      } 
    } else {
      stop("Error: Unknown standard error method (se.method).")
    }
    
    ## Response as predictions
    result$predictions <- yhat
    
    if (forest$treetype == "Probability estimation") {
      ## Set colnames and sort by levels
      colnames(result$predictions) <- forest$levels[forest$class.values]
      result$predictions <- result$predictions[, forest$levels, drop = FALSE]
      
      if (!is.matrix(result$se)) {
        result$se <- matrix(result$se, ncol = length(forest$levels))
      }
      colnames(result$se) <- forest$levels[forest$class.values]
      result$se <- result$se[, forest$levels, drop = FALSE]
    }
  }

  class(result) <- "ranger.prediction"
  return(result)
}

##' Prediction with new data and a saved forest from Ranger.
##' 
##' For \code{type = 'response'} (the default), the predicted classes (classification), predicted numeric values (regression), predicted probabilities (probability estimation) or survival probabilities (survival) are returned. 
##' For \code{type = 'se'}, the standard error of the predictions are returned (regression only). The jackknife-after-bootstrap or infinitesimal jackknife for bagging is used to estimate the standard errors based on out-of-bag predictions. See Wager et al. (2014) for details.
##' For \code{type = 'terminalNodes'}, the IDs of the terminal node in each tree for each observation in the given dataset are returned.
##' For \code{type = 'quantiles'}, the selected quantiles for each observation are estimated. See Meinshausen (2006) for details.
##' 
##' If \code{type = 'se'} is selected, the method to estimate the variances can be chosen with \code{se.method}. Set \code{se.method = 'jack'} for jackknife-after-bootstrap and \code{se.method = 'infjack'} for the infinitesimal jackknife for bagging.
##' 
##' For classification and \code{predict.all = TRUE}, a factor levels are returned as numerics.
##' To retrieve the corresponding factor levels, use \code{rf$forest$levels}, if \code{rf} is the ranger object.
##'
##' @title Ranger prediction
##' @param object Ranger \code{ranger} object.
##' @param data New test data of class \code{data.frame} or \code{gwaa.data} (GenABEL).
##' @param predict.all Return individual predictions for each tree instead of aggregated predictions for all trees. Return a matrix (sample x tree) for classification and regression, a 3d array for probability estimation (sample x class x tree) and survival (sample x time x tree).
##' @param num.trees Number of trees used for prediction. The first \code{num.trees} in the forest are used.
##' @param type Type of prediction. One of 'response', 'se', 'terminalNodes', 'quantiles' with default 'response'. See below for details.
##' @param se.method Method to compute standard errors. One of 'jack', 'infjack' with default 'infjack'. Only applicable if type = 'se'. See below for details.
##' @param quantiles Vector of quantiles for quantile prediction. Set \code{type = 'quantiles'} to use.
##' @param what User specified function for quantile prediction used instead of \code{quantile}. Must return numeric vector, see examples.
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
##' @examples
##' ## Classification forest
##' ranger(Species ~ ., data = iris)
##' train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
##' iris.train <- iris[train.idx, ]
##' iris.test <- iris[-train.idx, ]
##' rg.iris <- ranger(Species ~ ., data = iris.train)
##' pred.iris <- predict(rg.iris, data = iris.test)
##' table(iris.test$Species, pred.iris$predictions)
##' 
##' ## Quantile regression forest
##' rf <- ranger(mpg ~ ., mtcars[1:26, ], quantreg = TRUE)
##' pred <- predict(rf, mtcars[27:32, ], type = "quantiles", quantiles = c(0.1, 0.5, 0.9))
##' pred$predictions
##' 
##' ## Quantile regression forest with user-specified function
##' rf <- ranger(mpg ~ ., mtcars[1:26, ], quantreg = TRUE)
##' pred <- predict(rf, mtcars[27:32, ], type = "quantiles", 
##'                 what = function(x) sample(x, 10, replace = TRUE))
##' pred$predictions
##' 
##' @references
##' \itemize{
##'   \item Wright, M. N. & Ziegler, A. (2017). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. J Stat Softw 77:1-17. \doi{10.18637/jss.v077.i01}.
##'   \item Wager, S., Hastie T., & Efron, B. (2014). Confidence Intervals for Random Forests: The Jackknife and the Infinitesimal Jackknife. J Mach Learn Res 15:1625-1651. \url{https://jmlr.org/papers/v15/wager14a.html}.
##'   \item Meinshausen (2006). Quantile Regression Forests. J Mach Learn Res 7:983-999. \url{https://www.jmlr.org/papers/v7/meinshausen06a.html}.  
##'   }
##' @seealso \code{\link{ranger}}
##' @author Marvin N. Wright
##' @export
predict.ranger <- function(object, data = NULL, predict.all = FALSE,
                           num.trees = object$num.trees,
                           type = "response", se.method = "infjack",
                           quantiles = c(0.1, 0.5, 0.9), 
                           what = NULL,
                           seed = NULL, num.threads = NULL,
                           verbose = TRUE, ...) {
  forest <- object$forest
  if (is.null(forest)) {
    stop("Error: No saved forest in ranger object. Please set write.forest to TRUE when calling ranger.")
  }
  if (object$importance.mode %in% c("impurity_corrected", "impurity_unbiased")) {
    warning("Forest was grown with 'impurity_corrected' variable importance. For prediction it is advised to grow another forest without this importance setting.")
  }
  
  if (type == "quantiles") {
    ## Quantile prediction
    if (object$treetype != "Regression") {
      stop("Error: Quantile prediction implemented only for regression outcomes.")
    }
    if (is.null(object$random.node.values)) {
      stop("Error: Set quantreg=TRUE in ranger(...) for quantile prediction.")
    }
    
    if (is.null(data)) {
      ## OOB prediction
      if (is.null(object$random.node.values.oob)) {
        stop("Error: Set keep.inbag=TRUE in ranger(...) for out-of-bag quantile prediction or provide new data in predict(...).")
      }
      node.values <- object$random.node.values.oob
    } else {
      ## New data prediction
      terminal.nodes <- predict(object, data, num.threads = num.threads, type = "terminalNodes")$predictions + 1
      node.values <- 0 * terminal.nodes
      for (tree in 1:num.trees) {
        node.values[, tree] <- object$random.node.values[terminal.nodes[, tree], tree]
      }
    }
    
    ## Prepare results
    result <- list(num.samples = nrow(node.values),
                   treetype = object$treetype,
                   num.independent.variables = object$num.independent.variables,
                   num.trees = num.trees)
    class(result) <- "ranger.prediction"

    if (is.null(what)) {
      ## Compute quantiles of distribution
      result$predictions <- t(apply(node.values, 1, quantile, quantiles, na.rm=TRUE))
      if (nrow(result$predictions) != result$num.samples) {
        ## Fix result for single quantile
        result$predictions <- t(result$predictions)
      }
      colnames(result$predictions) <- paste("quantile=", quantiles)
    } else {
      ## User function
      if (!is.function(what)) {
        stop("Error: Argument 'what' is not a function.")
      }
      result$predictions <- t(apply(node.values, 1, what))
    }
   
    result
  } else {
    ## Non-quantile prediction
    if (is.null(data)) {
     stop("Error: Argument 'data' is required for non-quantile prediction.") 
    }
    predict(forest, data, predict.all, num.trees, type, se.method, seed, num.threads, verbose, object$inbag.counts, ...)
  }
}
