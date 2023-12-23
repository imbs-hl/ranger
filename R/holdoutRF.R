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

##' Grow two random forests on two cross-validation folds. 
##' Instead of out-of-bag data, the other fold is used to compute permutation importance.
##' Related to the novel permutation variable importance by Janitza et al. (2015).
##'
##' @title Hold-out random forests
##' @param ... All arguments are passed to \code{\link{ranger}()} (except \code{importance}, \code{case.weights}, \code{replace} and \code{holdout}.). 
##' @return Hold-out random forests with variable importance.
##' @seealso \code{\link{ranger}}
##' @author Marvin N. Wright
##' @references
##'   Janitza, S., Celik, E. & Boulesteix, A.-L., (2015). A computationally fast variable importance test for random forests for high-dimensional data. Adv Data Anal Classif \doi{10.1007/s11634-016-0276-4}. \cr
##' @export 
holdoutRF <- function(...) {
  
  ## Get data from arguments
  args <- list(...)
  if ("data" %in% names(args)) {
    data <- args$data
  } else {
    data <- args[[2]]
  }
  
  ## Split data
  if (inherits(data, "gwaa.data")) {
    n <- nrow(data@phdata) 
  } else {
    n <- nrow(data)
  }
  weights <- rbinom(n, 1, 0.5)
  
  ## Check args
  if ("case.weights" %in% names(args)) {
    stop("Error: Argument 'case.weights' not supported in holdoutRF.")
  }
  if ("holdout" %in% names(args)) {
    stop("Error: Argument 'holdout' not supported in holdoutRF.")
  }
  if ("importance" %in% names(args)) {
    stop("Error: Argument 'importance' not supported in holdoutRF. Always set to 'permutation'.")
  }
  if ("replace" %in% names(args)) {
    stop("Error: Argument 'replace' not supported in holdoutRF.")
  }
  
  ## Grow RFs
  res <- list(
    rf1 = ranger(..., importance = "permutation",  
                 case.weights = weights, replace = FALSE, holdout = TRUE),
    rf2 = ranger(..., importance = "permutation",
                 case.weights = 1-weights, replace = FALSE, holdout = TRUE)
  )
  
  ## Compute importance
  res$variable.importance <- (res$rf1$variable.importance + res$rf2$variable.importance)/2
  res$treetype <- res$rf1$treetype
  res$importance.mode <- res$rf1$importance.mode
  class(res) <- "holdoutRF"
  
  res
}
