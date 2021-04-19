# -------------------------------------------------------------------------------
#   This file is part of rangerts.
#
# rangerts is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# -------------------------------------------------------------------------------

##' @export
importance <- function(x, ...)  UseMethod("importance")

##' Extract variable importance of rangerts object.
##'
##'
##' @title rangerts variable importance
##' @param x rangerts object.
##' @param ... Further arguments passed to or from other methods.
##' @return Variable importance measures.
##' @seealso \code{\link{rangerts}}
##' @author Marvin N. Wright
##' @aliases importance
##' @export
importance.rangerts <- function(x, ...) {
  if (!inherits(x, "rangerts")) {
    stop("Object ist no rangerts object.")
  }
  if (is.null(x$variable.importance) || length(x$variable.importance) < 1) {
    stop("No variable importance found. Please use 'importance' option when growing the forest.")
  }
  return(x$variable.importance)
}

##' Compute variable importance with p-values.
##' For high dimensional data, the fast method of Janitza et al. (2016) can be used.
##' The permutation approach of Altmann et al. (2010) is computationally intensive but can be used with all kinds of data.
##' See below for details.
##'
##' The method of Janitza et al. (2016) uses a clever trick:
##' With an unbiased variable importance measure, the importance values of non-associated variables vary randomly around zero.
##' Thus, all non-positive importance values are assumed to correspond to these non-associated variables and they are used to construct a distribution of the importance under the null hypothesis of no association to the response.
##' Since only the non-positive values of this distribution can be observed, the positive values are created by mirroring the negative distribution.
##' See Janitza et al. (2016) for details.
##'
##' The method of Altmann et al. (2010) uses a simple permutation test:
##' The distribution of the importance under the null hypothesis of no association to the response is created by several replications of permuting the response, growing an RF and computing the variable importance.
##' The authors recommend 50-100 permutations.
##' However, much larger numbers have to be used to estimate more precise p-values.
##' We add 1 to the numerator and denominator to avoid zero p-values.
##'
##' @title rangerts variable importance p-values
##' @param x \code{rangerts} or \code{holdoutRF} object.
##' @param method Method to compute p-values. Use "janitza" for the method by Janitza et al. (2016) or "altmann" for the non-parametric method by Altmann et al. (2010).
##' @param num.permutations Number of permutations. Used in the "altmann" method only.
##' @param formula Object of class formula or character describing the model to fit. Used in the "altmann" method only.
##' @param data Training data of class data.frame or matrix. Used in the "altmann" method only.
##' @param ... Further arguments passed to \code{rangerts()}. Used in the "altmann" method only.
##' @return Variable importance and p-value for each variable.
##' @examples
##' ## Janitza's p-values with corrected Gini importance
##' n <- 50
##' p <- 400
##' dat <- data.frame(y = factor(rbinom(n, 1, .5)), replicate(p, runif(n)))
##' rf.sim <- rangerts(y ~ ., dat, importance = "impurity_corrected")
##' importance_pvalues(rf.sim, method = "janitza")
##'
##' ## Permutation p-values
##' \dontrun{
##' rf.iris <- rangerts(Species ~ ., data = iris, importance = 'permutation')
##' importance_pvalues(rf.iris, method = "altmann", formula = Species ~ ., data = iris)
##' }
##' @seealso \code{\link{rangerts}}
##' @author Marvin N. Wright
##' @references
##'   Janitza, S., Celik, E. & Boulesteix, A.-L., (2016). A computationally fast variable importance test for random forests for high-dimensional data. Adv Data Anal Classif \url{https://doi.org/10.1007/s11634-016-0276-4}. \cr
##'   Altmann, A., Tolosi, L., Sander, O. & Lengauer, T. (2010). Permutation importance: a corrected feature importance measure, Bioinformatics 26:1340-1347.
##' @export
importance_pvalues <- function(x, method = c("janitza", "altmann"), num.permutations = 100, formula = NULL, data = NULL, ...) {
  method <- match.arg(method)
  if (!inherits(x, c("rangerts", "holdoutRF"))) {
    stop("Object is no rangerts or holdoutRF object.")
  }
  if (x$importance.mode == "none" || is.null(x$variable.importance) || length(x$variable.importance) < 1) {
    stop("No variable importance found. Please use 'importance' option when growing the forest.")
  }

  if (method == "janitza") {
    if (x$importance.mode == "impurity") {
      stop("Impurity variable importance found. Please use (hold-out) permutation importance or corrected impurity importance to use this method.")
    }
    if (!inherits(x, "holdoutRF") && x$importance.mode == "permutation") {
      warning("Permutation variable importance found, inaccurate p-values. Please use hold-out permutation importance or corrected impurity importance to use this method.")
    }
    if (x$treetype != "Classification") {
      warning("This method is tested for classification only, use with care.")
    }

    ## Mirrored VIMP
    m1 <- x$variable.importance[x$variable.importance < 0]
    m2 <- x$variable.importance[x$variable.importance == 0]
    vimp <- c(m1, -m1, m2)

    ## Compute p-value
    ## Note: ecdf is smaller or equal, problems with 0 importance values
    pval <- 1 - numSmaller(x$variable.importance, vimp) / length(vimp)

    ## TODO: 100 ok? increase?
    if (length(m1) == 0) {
      stop("No negative importance values found. Consider the 'altmann' approach.")
    }
    if (length(m1) < 100) {
      warning("Only few negative importance values found, inaccurate p-values. Consider the 'altmann' approach.")
    }
  } else if (method == "altmann") {
    if (!inherits(x, "rangerts")) {
      stop("Altmann method not available for holdoutRF objects.")
    }
    if (is.null(formula) || is.null(data)) {
      stop("Formula and data required for the 'altmann' method.")
    }
    if (is.character(formula)) {
      formula <- formula(formula)
    }

    ## Permute and compute importance again
    if (x$treetype == "Survival") {
      dependent.variable.name <- all.vars(formula)[1:2]
    } else {
      dependent.variable.name <- all.vars(formula)[1]
    }
    vimp <- sapply(1:num.permutations, function(i) {
      dat <- data
      dat[, dependent.variable.name] <- dat[sample(nrow(dat)), dependent.variable.name]
      rangerts(formula, dat, num.trees = x$num.trees, mtry = x$mtry, min.node.size = x$min.node.size,
             importance = x$importance.mode, replace = x$replace, ...)$variable.importance
    })

    ## Compute p-value
    pval <- sapply(1:nrow(vimp), function(i) {
      (sum(vimp[i, ] >= x$variable.importance[i]) + 1)/(ncol(vimp) + 1)
    })

  } else {
    stop("Unknown p-value method. Available methods are: 'janitza' and 'altmann'.")
  }

  ## Return VIMP and p-values
  res <- cbind(x$variable.importance, pval)
  colnames(res) <- c("importance", "pvalue")
  return(res)
}
