# -------------------------------------------------------------------------------
#   This file is part of rangerts.
#
# rangerts is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# -------------------------------------------------------------------------------

##' @export
timepoints <- function(x, ...)  UseMethod("timepoints")

##' Extract unique death times of rangerts Survival prediction object.
##'
##'
##' @title rangerts timepoints
##' @param x rangerts Survival prediction object.
##' @param ... Further arguments passed to or from other methods.
##' @return Unique death times
##' @seealso \code{\link{rangerts}}
##' @author Marvin N. Wright
##' @export
timepoints.rangerts.prediction <- function(x, ...) {
  if (!inherits(x, "rangerts.prediction")) {
    stop("Object ist no rangerts.prediction object.")
  }
  if (x$treetype != "Survival") {
    stop("No timepoints found. Object is no Survival prediction object.")
  }
  if (is.null(x$unique.death.times)) {
    stop("No timepoints found.")
  }
  return(x$unique.death.times)
}

##' Extract unique death times of rangerts Survival forest
##'
##'
##' @title rangerts timepoints
##' @param x rangerts Survival forest object.
##' @param ... Further arguments passed to or from other methods.
##' @return Unique death times
##' @seealso \code{\link{rangerts}}
##' @aliases timepoints
##' @export
timepoints.rangerts <- function(x, ...) {
  if (!inherits(x, "rangerts")) {
    stop("Object ist no rangerts object.")
  }
  if (x$treetype != "Survival") {
    stop("No timepoints found. Object is no Survival forest.")
  }
  if (is.null(x$unique.death.times)) {
    stop("No timepoints found.")
  }
  return(x$unique.death.times)
}
