# -------------------------------------------------------------------------------
#   This file is part of rangerts.
#
# rangerts is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# -------------------------------------------------------------------------------

##' Print contents of rangerts object.
##'
##'
##' @title Print rangerts
##' @param x Object of class 'rangerts'.
##' @param ... Further arguments passed to or from other methods.
##' @seealso \code{\link{rangerts}}
##' @author Marvin N. Wright
##' @export
print.rangerts <- function(x, ...) {
  cat("rangerts result\n\n")
  cat("Call:\n", deparse(x$call), "\n\n")
  cat("Type:                            ", x$treetype, "\n")
  cat("Number of trees:                 ", x$num.trees, "\n")
  cat("Sample size:                     ", x$num.samples, "\n")
  cat("Number of independent variables: ", x$num.independent.variables, "\n")
  cat("Mtry:                            ", x$mtry, "\n")
  cat("Target node size:                ", x$min.node.size, "\n")
  cat("Variable importance mode:        ", x$importance.mode, "\n")
  cat("Splitrule:                       ", x$splitrule, "\n")
  if (x$treetype == "Survival") {
    cat("Number of unique death times:    ", length(x$unique.death.times), "\n")
  }
  if (!is.null(x$splitrule) && x$splitrule == "extratrees" && !is.null(x$num.random.splits)) {
    cat("Number of random splits:         ", x$num.random.splits, "\n")
  }
  if (x$treetype == "Classification") {
    cat("OOB prediction error:            ", sprintf("%1.2f %%", 100*x$prediction.error), "\n")
  } else if (x$treetype == "Regression") {
    cat("OOB prediction error (MSE):      ", x$prediction.error, "\n")
  } else if (x$treetype == "Survival") {
    cat("OOB prediction error (1-C):      ", x$prediction.error, "\n")
  } else if (x$treetype == "Probability estimation") {
    cat("OOB prediction error (Brier s.): ", x$prediction.error, "\n")
  } else {
    cat("OOB prediction error:            ", x$prediction.error, "\n")
  }
  if (x$treetype == "Regression") {
    cat("R squared (OOB):                 ", x$r.squared, "\n")
  }
}

##' Print contents of rangerts forest object.
##'
##'
##' @title Print rangerts forest
##' @param x Object of class 'rangerts.forest'.
##' @param ... further arguments passed to or from other methods.
##' @author Marvin N. Wright
##' @export
print.rangerts.forest <- function(x, ...) {
  cat("rangerts forest object\n\n")
  cat("Type:                         ", x$treetype, "\n")
  cat("Number of trees:              ", x$num.trees, "\n")
  if (x$treetype == "Survival") {
    cat("Number of unique death times: ", length(x$unique.death.times), "\n")
  }
}

##' Print contents of rangerts prediction object.
##'
##'
##' @title Print rangerts prediction
##' @param x Object of class 'rangerts.prediction'.
##' @param ... further arguments passed to or from other methods.
##' @author Marvin N. Wright
##' @export
print.rangerts.prediction <- function(x, ...) {
  cat("rangerts prediction\n\n")
  cat("Type:                            ", x$treetype, "\n")
  cat("Sample size:                     ", x$num.samples, "\n")
  cat("Number of independent variables: ", x$num.independent.variables, "\n")
  if (x$treetype == "Survival") {
    cat("Number of unique death times:    ", length(x$unique.death.times), "\n")
  }
}

str.rangerts.forest <- function(object, max.level = 2, ...) {
  class(object) <- "list"
  str(object, max.level = max.level, ...)
}

str.rangerts <- function(object, max.level = 2, ...) {
  class(object) <- "list"
  str(object, max.level = max.level, ...)
}
