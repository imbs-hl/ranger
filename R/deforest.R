#' Deforesting a random forest
#' 
#' The main purpose of this function is to allow for post-processing of 
#' ensembles via L2 regularized regression (i.e., the LASSO), as described in 
#' Friedman and Popescu (2003). The basic idea is to use the LASSO to 
#' post-process the predictions from the individual base learners in an ensemble 
#' (i.e., decision trees) in the hopes of producing a much smaller model without 
#' sacrificing much in the way of accuracy, and in some cases, improving it. 
#' Friedman and Popescu (2003) describe conditions under which tree-based 
#' ensembles, like random forest, can potentially benefit from such 
#' post-processing (e.g., using shallower trees trained on much smaller samples 
#' of the training data without replacement). However, the computational 
#' benefits of such post-processing can only be realized if the base learners 
#' "zeroed out" by the LASSO can actually be removed from the original ensemble, 
#' hence the purpose of this function. A complete example using 
#' \code{\link{ranger}} can be found at 
#' \url{https://github.com/imbs-hl/ranger/issues/568}.
#' 
#' @param object A fitted random forest (e.g., a \code{\link{ranger}}
#' object).
#' 
#' @param which.trees Vector giving the indices of the trees to remove.
#' 
#' @param warn Logical indicating whether or not to warn users that some of the 
#' standard output of a typical \code{\link{ranger}} object or no longer 
#' available after deforestation. Default is \code{TRUE}.
#' 
#' @param ... Additional (optional) arguments. (Currently ignored.)
#' 
#' @return An object of class \code{"deforest.ranger"}; essentially, a 
#' \code{\link{ranger}} object with certain components replaced with 
#' \code{NA}s (e.g., out-of-bag (OOB) predictions, variable importance scores 
#' (if requested), and OOB-based error metrics).
#' 
#' @note This function is a generic and can be extended by other packages.
#' 
#' @references 
#' Friedman, J. and Popescu, B. (2003). Importance sampled learning ensembles, 
#' Technical report, Stanford University, Department of Statistics.
#' \url{https://jerryfriedman.su.domains/ftp/isle.pdf}.
#' 
#' @rdname deforest
#' 
#' @export
#' 
#' @author Brandon M. Greenwell
#' 
#' @examples
#' ## Example of deforesting a random forest
#' rfo <- ranger(Species ~ ., data = iris, probability = TRUE, num.trees = 100)
#' dfo <- deforest(rfo, which.trees = c(1, 3, 5))
#' dfo  # same as `rfo` but with trees 1, 3, and 5 removed
#' 
#' ## Sanity check
#' preds.rfo <- predict(rfo, data = iris, predict.all = TRUE)$predictions
#' preds.dfo <- predict(dfo, data = iris, predict.all = TRUE)$predictions
#' identical(preds.rfo[, , -c(1, 3, 5)], y = preds.dfo)
deforest <- function(object, which.trees = NULL, ...) {
  UseMethod("deforest")
}


#' @rdname deforest
#' 
#' @export
deforest.ranger <- function(object, which.trees = NULL, warn = TRUE, ...) {
  
  # Warn users about `predictions` and `prediction.error` components
  if (isTRUE(warn)) {
    warning("Many of the components of a typical \"ranger\" object are ",
            "not available after deforestation and are instead replaced with ",
            "`NA` (e.g., out-of-bag (OOB) predictions, variable importance ",
            "scores (if requested), and OOB-based error metrics).", 
            call. = FALSE)
  }
  
  # "Remove trees" by removing necessary components from `forest` object
  object$forest$child.nodeIDs[which.trees] <- NULL
  object$forest$split.values[which.trees] <- NULL
  object$forest$split.varIDs[which.trees] <- NULL
  object$forest$terminal.class.counts[which.trees] <- NULL  # for prob forests
  object$forest$chf[which.trees] <- NULL  # for survival forests

  # Update `num.trees` components so `predict.ranger()` works
  object$forest$num.trees <- object$num.trees <- 
    length(object$forest$child.nodeIDs)
  
  # Coerce other components to `NA` as needed
  if (!is.null(object$prediction.error)) {
    object$prediction.error <- NA
  } 
  if (!is.null(object$predictions)) {  # classification and regression
    object$predictions[] <- NA
  }
  if (!is.null(object$r.squared)) {  # regression
    object$r.squared <- NA
  }
  if (!is.null(object$chf)) {  # survival forests
    object$chf[] <- NA
  }
  if (!is.null(object$survival)) {  # survival forests
    object$survival[] <- NA
  }
  if (object$importance.mode != "none") {  # variable importance
    object$importance.mode <- NA
    object$variable.importance[] <- NA
  }
  
  # Return "deforested" forest
  class(object) <- c("deforest.ranger", class(object))
  object
  
}


#' Print deforested ranger summary
#' 
#' Print basic information about a deforested \code{\link{ranger}} object.
#' 
#' @param x A \code{\link{deforest}} object (i.e., an object that inherits from
#' class \code{"deforest.ranger"}).
#' 
#' @param ... Further arguments passed to or from other methods.
#' 
#' @note Many of the components of a typical \code{\link{ranger}} object are not 
#' available after deforestation and are instead replaced with \code{NA} (e.g., 
#' out-of-bag (OOB) predictions, variable importance scores (if requested), and 
#' OOB-based error metrics).
#' 
#' @seealso \code{\link{deforest}}.
#' 
#' @author Brandon M. Greenwell
#' 
#' @export
print.deforest.ranger <- function (x, ...) {
  cat("Ranger (deforested) result\n\n")
  cat("Note that many of the components of a typical \"ranger\" object are",
      "not available after deforestation and are instead replaced with `NA`",
      "(e.g., out-of-bag (OOB) predictions, variable importance scores (if",
      "requested), and OOB-based error metrics)", 
      "\n\n")
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
  if (!is.null(x$splitrule) && x$splitrule == "extratrees" && 
      !is.null(x$num.random.splits)) {
    cat("Number of random splits:         ", x$num.random.splits, "\n")
  }
  if (x$treetype == "Classification") {
    cat("OOB prediction error:            ", x$prediction.error, "\n")
  }
  else if (x$treetype == "Regression") {
    cat("OOB prediction error (MSE):      ", x$prediction.error, "\n")
  }
  else if (x$treetype == "Survival") {
    cat("OOB prediction error (1-C):      ", x$prediction.error, "\n")
  }
  else if (x$treetype == "Probability estimation") {
    cat("OOB prediction error (Brier s.): ", x$prediction.error, "\n")
  }
  else {
    cat("OOB prediction error:            ", x$prediction.error, "\n")
  }
  if (x$treetype == "Regression") {
    cat("R squared (OOB):                 ", x$r.squared, "\n")
  }
}
