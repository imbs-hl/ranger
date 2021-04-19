# -------------------------------------------------------------------------------
#   This file is part of rangerts.
#
# rangerts is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# -------------------------------------------------------------------------------

##' This function is deprecated.
##' Please use predict() with \code{type = "terminalNodes"} instead.
##' This function calls predict() now.
##'
##' @title Get terminal node IDs (deprecated)
##' @param rf \code{rangerts} object.
##' @param dat New dataset. Terminal node IDs for this dataset are obtained.
##'
##' @return Matrix with terminal nodeIDs for all observations in dataset and trees.
##'
##' @examples
##' rf <- rangerts(Species ~ ., data = iris, num.trees = 5, write.forest = TRUE)
##' getTerminalNodeIDs(rf, iris)
##' @export
getTerminalNodeIDs <- function(rf, dat) {
  warning("Function getTerminalNodeIDs() deprecated, calling predict().")
  predict(rf, dat, type = "terminalNodes")$predictions
}
