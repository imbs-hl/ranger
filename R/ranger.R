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

##' Ranger is a fast implementation of random forests (Breiman 2001) or recursive partitioning, particularly suited for high dimensional data.
##' Classification, regression, and survival forests are supported.
##' Classification and regression forests are implemented as in the original Random Forest (Breiman 2001), survival forests as in Random Survival Forests (Ishwaran et al. 2008).
##' Includes implementations of extremely randomized trees (Geurts et al. 2006) and quantile regression forests (Meinshausen 2006). 
##'
##' The tree type is determined by the type of the dependent variable.
##' For factors classification trees are grown, for numeric values regression trees and for survival objects survival trees.
##' The Gini index is used as default splitting rule for classification.
##' For regression, the estimated response variances or maximally selected rank statistics (Wright et al. 2016) can be used.
##' For Survival the log-rank test, a C-index based splitting rule (Schmid et al. 2015) and maximally selected rank statistics (Wright et al. 2016) are available.
##' For all tree types, forests of extremely randomized trees (Geurts et al. 2006) can be grown.
##'
##' With the \code{probability} option and factor dependent variable a probability forest is grown.
##' Here, the node impurity is used for splitting, as in classification forests.
##' Predictions are class probabilities for each sample.
##' In contrast to other implementations, each tree returns a probability estimate and these estimates are averaged for the forest probability estimate.
##' For details see Malley et al. (2012).
##'
##' Note that for classification and regression nodes with size smaller than \code{min.node.size} can occur, as in original Random Forests.
##' For survival all nodes contain at \code{min.node.size} samples. 
##' Variables selected with \code{always.split.variables} are tried additionally to the mtry variables randomly selected.
##' In \code{split.select.weights}, weights do not need to sum up to 1, they will be normalized later. 
##' The weights are assigned to the variables in the order they appear in the formula or in the data if no formula is used.
##' Names of the \code{split.select.weights} vector are ignored.
##' The usage of \code{split.select.weights} can increase the computation times for large forests.
##'
##' Unordered factor covariates can be handled in 3 different ways by using \code{respect.unordered.factors}: 
##' For 'ignore' all factors are regarded ordered, for 'partition' all possible 2-partitions are considered for splitting. 
##' For 'order' and 2-class classification the factor levels are ordered by their proportion falling in the second class, for regression by their mean response, as described in Hastie et al. (2009), chapter 9.2.4.
##' For multiclass classification the factor levels are ordered by the first principal component of the weighted covariance matrix of the contingency table (Coppersmith et al. 1999), for survival by the median survival (or the largest available quantile if the median is not available).
##' The use of 'order' is recommended, as it computationally fast and can handle an unlimited number of factor levels. 
##' Note that the factors are only reordered once and not again in each split. 
##'
##' The 'impurity_corrected' importance measure is unbiased in terms of the number of categories and category frequencies and is almost as fast as the standard impurity importance.
##' It is a modified version of the method by Sandri & Zuccolotto (2008), which is faster and more memory efficient. 
##' See Nembrini et al. (2018) for details.
##' This importance measure can be combined with the methods to estimate p-values in \code{\link{importance_pvalues}}.
##'
##' Regularization works by penalizing new variables by multiplying the splitting criterion by a factor, see Deng & Runger (2012) for details.  
##' If \code{regularization.usedepth=TRUE}, \eqn{f^d} is used, where \emph{f} is the regularization factor and \emph{d} the depth of the node.
##' If regularization is used, multithreading is deactivated because all trees need access to the list of variables that are already included in the model.
##'
##' For a large number of variables and data frames as input data the formula interface can be slow or impossible to use.
##' Alternatively \code{dependent.variable.name} (and \code{status.variable.name} for survival) or \code{x} and \code{y} can be used.
##' Use \code{x} and \code{y} with a matrix for \code{x} to avoid conversions and save memory.
##' Consider setting \code{save.memory = TRUE} if you encounter memory problems for very large datasets, but be aware that this option slows down the tree growing. 
##' 
##' For GWAS data consider combining \code{ranger} with the \code{GenABEL} package. 
##' See the Examples section below for a demonstration using \code{Plink} data.
##' All SNPs in the \code{GenABEL} object will be used for splitting. 
##' To use only the SNPs without sex or other covariates from the phenotype file, use \code{0} on the right hand side of the formula. 
##' Note that missing values are treated as an extra category while splitting.
##' 
##' See \url{https://github.com/imbs-hl/ranger} for the development version.
##' 
##' With recent R versions, multithreading on Windows platforms should just work. 
##' If you compile yourself, the new RTools toolchain is required.
##' 
##' @title Ranger
##' @param formula Object of class \code{formula} or \code{character} describing the model to fit. Interaction terms supported only for numerical variables.
##' @param data Training data of class \code{data.frame}, \code{matrix}, \code{dgCMatrix} (Matrix) or \code{gwaa.data} (GenABEL).
##' @param num.trees Number of trees.
##' @param mtry Number of variables to possibly split at in each node. Default is the (rounded down) square root of the number variables. Alternatively, a single argument function returning an integer, given the number of independent variables.
##' @param importance Variable importance mode, one of 'none', 'impurity', 'impurity_corrected', 'permutation'. The 'impurity' measure is the Gini index for classification, the variance of the responses for regression and the sum of test statistics (see \code{splitrule}) for survival. 
##' @param write.forest Save \code{ranger.forest} object, required for prediction. Set to \code{FALSE} to reduce memory usage if no prediction intended.
##' @param probability Grow a probability forest as in Malley et al. (2012). 
##' @param min.node.size Minimal node size. Default 1 for classification, 5 for regression, 3 for survival, and 10 for probability.
##' @param max.depth Maximal tree depth. A value of NULL or 0 (the default) corresponds to unlimited depth, 1 to tree stumps (1 split per tree).
##' @param replace Sample with replacement. 
##' @param sample.fraction Fraction of observations to sample. Default is 1 for sampling with replacement and 0.632 for sampling without replacement. For classification, this can be a vector of class-specific values. 
##' @param case.weights Weights for sampling of training observations. Observations with larger weights will be selected with higher probability in the bootstrap (or subsampled) samples for the trees.
##' @param class.weights Weights for the outcome classes (in order of the factor levels) in the splitting rule (cost sensitive learning). Classification and probability prediction only. For classification the weights are also applied in the majority vote in terminal nodes.
##' @param splitrule Splitting rule. For classification and probability estimation "gini", "extratrees" or "hellinger" with default "gini".
##'   For regression "variance", "extratrees", "maxstat", "beta" or "poisson" with default "variance".
##'   For survival "logrank", "extratrees", "C" or "maxstat" with default "logrank". 
##' @param num.random.splits For "extratrees" splitrule.: Number of random splits to consider for each candidate splitting variable.
##' @param alpha For "maxstat" splitrule: Significance threshold to allow splitting.
##' @param minprop For "maxstat" splitrule: Lower quantile of covariate distribution to be considered for splitting.
##' @param poisson.tau For "poisson" splitrule: The coefficient of variation of the (expected) frequency is \eqn{1/\tau}.
##'   If a terminal node has only 0 responses, the estimate is set to \eqn{\alpha 0 + (1-\alpha) mean(parent)} with \eqn{\alpha = samples(child) mean(parent) / (\tau + samples(child) mean(parent))}.
##' @param split.select.weights Numeric vector with weights between 0 and 1, representing the probability to select variables for splitting. Alternatively, a list of size num.trees, containing split select weight vectors for each tree can be used.  
##' @param always.split.variables Character vector with variable names to be always selected in addition to the \code{mtry} variables tried for splitting.
##' @param respect.unordered.factors Handling of unordered factor covariates. One of 'ignore', 'order' and 'partition'. For the "extratrees" splitrule the default is "partition" for all other splitrules 'ignore'. Alternatively TRUE (='order') or FALSE (='ignore') can be used. See below for details. 
##' @param scale.permutation.importance Scale permutation importance by standard error as in (Breiman 2001). Only applicable if permutation variable importance mode selected.
##' @param regularization.factor Regularization factor (gain penalization), either a vector of length p or one value for all variables.
##' @param regularization.usedepth Consider the depth in regularization. 
##' @param local.importance Calculate and return local importance values as in (Breiman 2001). Only applicable if \code{importance} is set to 'permutation'.
##' @param keep.inbag Save how often observations are in-bag in each tree. 
##' @param inbag Manually set observations per tree. List of size num.trees, containing inbag counts for each observation. Can be used for stratified sampling.
##' @param holdout Hold-out mode. Hold-out all samples with case weight 0 and use these for variable importance and prediction error.
##' @param quantreg Prepare quantile prediction as in quantile regression forests (Meinshausen 2006). Regression only. Set \code{keep.inbag = TRUE} to prepare out-of-bag quantile prediction.
##' @param oob.error Compute OOB prediction error. Set to \code{FALSE} to save computation time, e.g. for large survival forests.
##' @param num.threads Number of threads. Default is number of CPUs available.
##' @param save.memory Use memory saving (but slower) splitting mode. No effect for survival and GWAS data. Warning: This option slows down the tree growing, use only if you encounter memory problems.
##' @param verbose Show computation status and estimated runtime.
##' @param seed Random seed. Default is \code{NULL}, which generates the seed from \code{R}. Set to \code{0} to ignore the \code{R} seed. 
##' @param dependent.variable.name Name of dependent variable, needed if no formula given. For survival forests this is the time variable.
##' @param status.variable.name Name of status variable, only applicable to survival data and needed if no formula given. Use 1 for event and 0 for censoring.
##' @param classification Set to \code{TRUE} to grow a classification forest. Only needed if the data is a matrix or the response numeric. 
##' @param x Predictor data (independent variables), alternative interface to data with formula or dependent.variable.name.
##' @param y Response vector (dependent variable), alternative interface to data with formula or dependent.variable.name. For survival use a \code{Surv()} object or a matrix with time and status.
##' @return Object of class \code{ranger} with elements
##'   \item{\code{forest}}{Saved forest (If write.forest set to TRUE). Note that the variable IDs in the \code{split.varIDs} object do not necessarily represent the column number in R.}
##'   \item{\code{predictions}}{Predicted classes/values, based on out of bag samples (classification and regression only).}
##'   \item{\code{variable.importance}}{Variable importance for each independent variable.}
##'   \item{\code{variable.importance.local}}{Variable importance for each independent variable and each sample, if \code{local.importance} is set to TRUE and \code{importance} is set to 'permutation'.}
##'   \item{\code{prediction.error}}{Overall out of bag prediction error. For classification this is the fraction of missclassified samples, for probability estimation the Brier score, for regression the mean squared error and for survival one minus Harrell's C-index.}
##'   \item{\code{r.squared}}{R squared. Also called explained variance or coefficient of determination (regression only). Computed on out of bag data.}
##'   \item{\code{confusion.matrix}}{Contingency table for classes and predictions based on out of bag samples (classification only).}
##'   \item{\code{unique.death.times}}{Unique death times (survival only).}
##'   \item{\code{chf}}{Estimated cumulative hazard function for each sample (survival only).}
##'   \item{\code{survival}}{Estimated survival function for each sample (survival only).}
##'   \item{\code{call}}{Function call.}
##'   \item{\code{num.trees}}{Number of trees.}
##'   \item{\code{num.independent.variables}}{Number of independent variables.}
##'   \item{\code{mtry}}{Value of mtry used.}
##'   \item{\code{min.node.size}}{Value of minimal node size used.}
##'   \item{\code{treetype}}{Type of forest/tree. classification, regression or survival.}
##'   \item{\code{importance.mode}}{Importance mode used.}
##'   \item{\code{num.samples}}{Number of samples.}
##'   \item{\code{inbag.counts}}{Number of times the observations are in-bag in the trees.}
##' @examples
##' ## Classification forest with default settings
##' ranger(Species ~ ., data = iris)
##'
##' ## Prediction
##' train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
##' iris.train <- iris[train.idx, ]
##' iris.test <- iris[-train.idx, ]
##' rg.iris <- ranger(Species ~ ., data = iris.train)
##' pred.iris <- predict(rg.iris, data = iris.test)
##' table(iris.test$Species, pred.iris$predictions)
##' 
##' ## Quantile regression forest
##' rf <- ranger(mpg ~ ., mtcars[1:26, ], quantreg = TRUE)
##' pred <- predict(rf, mtcars[27:32, ], type = "quantiles")
##' pred$predictions
##'
##' ## Variable importance
##' rg.iris <- ranger(Species ~ ., data = iris, importance = "impurity")
##' rg.iris$variable.importance
##'
##' ## Survival forest
##' require(survival)
##' rg.veteran <- ranger(Surv(time, status) ~ ., data = veteran)
##' plot(rg.veteran$unique.death.times, rg.veteran$survival[1,])
##'
##' ## Alternative interfaces (same results)
##' ranger(dependent.variable.name = "Species", data = iris)
##' ranger(y = iris[, 5], x = iris[, -5])
##' 
##' \dontrun{
##' ## Use GenABEL interface to read Plink data into R and grow a classification forest
##' ## The ped and map files are not included
##' library(GenABEL)
##' convert.snp.ped("data.ped", "data.map", "data.raw")
##' dat.gwaa <- load.gwaa.data("data.pheno", "data.raw")
##' phdata(dat.gwaa)$trait <- factor(phdata(dat.gwaa)$trait)
##' ranger(trait ~ ., data = dat.gwaa)
##' }
##'
##' @author Marvin N. Wright
##' @references
##' \itemize{
##'   \item Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. \url{https://doi.org/10.18637/jss.v077.i01}.
##'   \item Schmid, M., Wright, M. N. & Ziegler, A. (2016). On the use of Harrell's C for clinical risk prediction via random survival forests. Expert Syst Appl 63:450-459. \url{https://doi.org/10.1016/j.eswa.2016.07.018}. 
##'   \item Wright, M. N., Dankowski, T. & Ziegler, A. (2017). Unbiased split variable selection for random survival forests using maximally selected rank statistics. Stat Med 36:1272-1284. \url{https://doi.org/10.1002/sim.7212}.
##'   \item Nembrini, S., Koenig, I. R. & Wright, M. N. (2018). The revival of the Gini Importance? Bioinformatics. \url{https://doi.org/10.1093/bioinformatics/bty373}.
##'   \item Breiman, L. (2001). Random forests. Mach Learn, 45:5-32. \url{https://doi.org/10.1023/A:1010933404324}. 
##'   \item Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random survival forests. Ann Appl Stat 2:841-860. \url{https://doi.org/10.1097/JTO.0b013e318233d835}. 
##'   \item Malley, J. D., Kruppa, J., Dasgupta, A., Malley, K. G., & Ziegler, A. (2012). Probability machines: consistent probability estimation using nonparametric learning machines. Methods Inf Med 51:74-81. \url{https://doi.org/10.3414/ME00-01-0052}.
##'   \item Hastie, T., Tibshirani, R., Friedman, J. (2009). The Elements of Statistical Learning. Springer, New York. 2nd edition.
##'   \item Geurts, P., Ernst, D., Wehenkel, L. (2006). Extremely randomized trees. Mach Learn 63:3-42. \url{https://doi.org/10.1007/s10994-006-6226-1}.
##'   \item Meinshausen (2006). Quantile Regression Forests. J Mach Learn Res 7:983-999. \url{http://www.jmlr.org/papers/v7/meinshausen06a.html}.  
##'   \item Sandri, M. & Zuccolotto, P. (2008). A bias correction algorithm for the Gini variable importance measure in classification trees. J Comput Graph Stat, 17:611-628. \url{https://doi.org/10.1198/106186008X344522}.
##'   \item Coppersmith D., Hong S. J., Hosking J. R. (1999). Partitioning nominal attributes in decision trees. Data Min Knowl Discov 3:197-217. \url{https://doi.org/10.1023/A:1009869804967}.
##'   \item Deng & Runger (2012). Feature selection via regularized trees. The 2012 International Joint Conference on Neural Networks (IJCNN), Brisbane, Australia. \url{https://doi.org/10.1109/IJCNN.2012.6252640}.
##'   }
##' @seealso \code{\link{predict.ranger}}
##' @useDynLib ranger, .registration = TRUE
##' @importFrom Rcpp evalCpp
##' @import stats 
##' @import utils
##' @importFrom Matrix Matrix
##' @export
ranger <- function(formula = NULL, data = NULL, num.trees = 500, mtry = NULL,
                   importance = "none", write.forest = TRUE, probability = FALSE,
                   min.node.size = NULL, max.depth = NULL, replace = TRUE, 
                   sample.fraction = ifelse(replace, 1, 0.632), 
                   case.weights = NULL, class.weights = NULL, splitrule = NULL, 
                   num.random.splits = 1, alpha = 0.5, minprop = 0.1,
                   poisson.tau = 1,
                   split.select.weights = NULL, always.split.variables = NULL,
                   respect.unordered.factors = NULL,
                   scale.permutation.importance = FALSE,
                   local.importance = FALSE, 
                   regularization.factor = 1, regularization.usedepth = FALSE,
                   keep.inbag = FALSE, inbag = NULL, holdout = FALSE,
                   quantreg = FALSE, oob.error = TRUE,
                   num.threads = NULL, save.memory = FALSE,
                   verbose = TRUE, seed = NULL, 
                   dependent.variable.name = NULL, status.variable.name = NULL, 
                   classification = NULL, x = NULL, y = NULL) {
  
  ## By default not in GWAS mode
  snp.data <- as.matrix(0)
  gwa.mode <- FALSE
  
  if (is.null(data)) {
    ## x/y interface
    if (is.null(x) | is.null(y)) {
      stop("Error: Either data or x and y is required.")
    }
  }  else {
    ## GenABEL GWA data
    if (inherits(data, "gwaa.data" )) {
      snp.names <- data@gtdata@snpnames
      snp.data <- data@gtdata@gtps@.Data
      data <- data@phdata
      if ("id" %in% names(data)) {
        data$"id" <- NULL
      }
      gwa.mode <- TRUE
      save.memory <- FALSE
    } 
    
    ## Formula interface. Use whole data frame if no formula provided and depvarname given
    if (is.null(formula)) {
      if (is.null(dependent.variable.name)) {
        if (is.null(y) | is.null(x)) {
          stop("Error: Please give formula, dependent variable name or x/y.")
        } 
      } else {
        if (is.null(status.variable.name)) {
          y <- data[, dependent.variable.name, drop = TRUE]
          x <- data[, !(colnames(data) %in% dependent.variable.name), drop = FALSE]
        } else {
          y <- survival::Surv(data[, dependent.variable.name], data[, status.variable.name]) 
          x <- data[, !(colnames(data) %in% c(dependent.variable.name, status.variable.name)), drop = FALSE]
        }
      }
    } else {
      formula <- formula(formula)
      if (!inherits(formula, "formula")) {
        stop("Error: Invalid formula.")
      }
      data.selected <- parse.formula(formula, data, env = parent.frame())
      y <- data.selected[, 1]
      x <- data.selected[, -1, drop = FALSE]
    }
  }
  
  ## Sparse matrix data
  if (inherits(x, "Matrix")) {
    if (!inherits(x, "dgCMatrix")) {
      stop("Error: Currently only sparse data of class 'dgCMatrix' supported.")
    } 
    if (!is.null(formula)) {
      stop("Error: Sparse matrices only supported with alternative interface. Use dependent.variable.name or x/y instead of formula.")
    }
  }
  
  ## Check missing values
  if (any(is.na(x))) {
    offending_columns <- colnames(x)[colSums(is.na(x)) > 0]
    stop("Missing data in columns: ",
         paste0(offending_columns, collapse = ", "), ".", call. = FALSE)
  }
  if (any(is.na(y))) {
    stop("Missing data in dependent variable.", call. = FALSE)
  }
  
  ## Check response levels
  if (is.factor(y)) {
    if (nlevels(y) != nlevels(droplevels(y))) {
      dropped_levels <- setdiff(levels(y), levels(droplevels(y)))
      warning("Dropped unused factor level(s) in dependent variable: ",
              paste0(dropped_levels, collapse = ", "), ".", call. = FALSE)
    }
  }
  
  ## Treetype
  if (is.factor(y) || is.logical(y)) {
    if (probability) {
      treetype <- 9
    } else {
      treetype <- 1
    }
  } else if (is.numeric(y) && (is.null(ncol(y)) || ncol(y) == 1)) {
    if (!is.null(classification) && classification && !probability) {
      treetype <- 1
    } else if (probability) {
      treetype <- 9
    } else {
      treetype <- 3
    }
  } else if (inherits(y, "Surv") || is.data.frame(y) || is.matrix(y)) {
    treetype <- 5
  } else {
    stop("Error: Unsupported type of dependent variable.")
  }
  
  ## Quantile prediction only for regression
  if (quantreg && treetype != 3) {
    stop("Error: Quantile prediction implemented only for regression outcomes.")
  }

  independent.variable.names <- colnames(x)
  
  ## respect.unordered.factors
  if (is.null(respect.unordered.factors)) {
    if (!is.null(splitrule) && splitrule == "extratrees") {
      respect.unordered.factors <- "partition"
    } else {
      respect.unordered.factors <- "ignore"
    }
  }

  ## Old version of respect.unordered.factors
  if (respect.unordered.factors == TRUE) {
    respect.unordered.factors <- "order"
  } else if (respect.unordered.factors == FALSE) {
    respect.unordered.factors <- "ignore"
  }
  
  ## Recode characters as factors and recode factors if 'order' mode
  if (!is.matrix(x) && !inherits(x, "Matrix") && ncol(x) > 0) {
    character.idx <- sapply(x, is.character)
    
    if (respect.unordered.factors == "order") {
      ## Recode characters and unordered factors
      ordered.idx <- sapply(x, is.ordered)
      factor.idx <- sapply(x, is.factor)
      recode.idx <- character.idx | (factor.idx & !ordered.idx)

      if (any(recode.idx) & (importance == "impurity_corrected" || importance == "impurity_unbiased")) {
        warning("Corrected impurity importance may not be unbiased for re-ordered factor levels. Consider setting respect.unordered.factors to 'ignore' or 'partition' or manually compute corrected importance.")
      }
      
      ## Numeric response
      if (is.factor(y)) {
        num.y <- as.numeric(y)
      } else {
        num.y <- y
      }

      ## Recode each column
      x[recode.idx] <- lapply(x[recode.idx], function(xx) {
        if (!is.factor(xx)) {
          xx <- as.factor(xx)
        } 
        
        if (length(levels(xx)) == 1) {
          ## Don't order if only one level
          levels.ordered <- levels(xx)
        } else if (inherits(y, "Surv")) {
          ## Use median survival if available or largest quantile available in all strata if median not available
          levels.ordered <- largest.quantile(y ~ xx)
          
          ## Get all levels not in node
          levels.missing <- setdiff(levels(xx), levels.ordered)
          levels.ordered <- c(levels.missing, levels.ordered)
        } else if (is.factor(y) & nlevels(y) > 2) {
          levels.ordered <- pca.order(y = y, x = xx)
        } else {
          ## Order factor levels by mean response
          means <- sapply(levels(xx), function(y) {
            mean(num.y[xx == y])
          })
          levels.ordered <- as.character(levels(xx)[order(means)])
        }
        
        ## Return reordered factor
        factor(xx, levels = levels.ordered, ordered = TRUE, exclude = NULL)
      })
      
      ## Save levels
      covariate.levels <- lapply(x, levels)
    } else {
      ## Recode characters only
      x[character.idx] <- lapply(x[character.idx], factor)
    }
  }
  
  ## If gwa mode, add snp variable names
  if (gwa.mode) {
    all.independent.variable.names <- c(independent.variable.names, snp.names)
  } else {
    all.independent.variable.names <- independent.variable.names
  }
  
  ## Error if no covariates
  if (length(all.independent.variable.names) < 1) {
    stop("Error: No covariates found.")
  }
  
  ## Number of trees
  if (!is.numeric(num.trees) || num.trees < 1) {
    stop("Error: Invalid value for num.trees.")
  }
  
  ## mtry as a function
  if (is.function(mtry)) { 
    nv <- length(all.independent.variable.names)
    
    if (length(formals(mtry)) > 1){
      stop("Error: Given mtry function requires single argument (the number of independent variables in the model).")
    }
    
    # Evaluate function
    mtry <- try(mtry(nv), silent = TRUE)
    
    if (inherits(mtry, "try-error")) {
      message("The mtry function produced the error: ", mtry)
      stop("Error: mtry function evaluation resulted in an error.")
    }
    
    ## Check for a single numeric
    if (!is.numeric(mtry) || length(mtry) != 1) {
      stop("Error: Given mtry function should return a single integer or numeric.")
    } else {
      mtry <- as.integer(mtry)
    }
    
    ## Check for limits
    if (mtry < 1 || mtry > nv) {
      stop("Error: Given mtry function should evaluate to a value not less than 1 and not greater than the number of independent variables ( = ", nv, " )")
    }
  }
  
  if (is.null(mtry)) {
    mtry <- 0
  } else if (!is.numeric(mtry) || mtry < 0) {
    stop("Error: Invalid value for mtry")
  }
  
  ## Seed
  if (is.null(seed)) {
    seed <- runif(1 , 0, .Machine$integer.max)
  }
  
  ## Keep inbag
  if (!is.logical(keep.inbag)) {
    stop("Error: Invalid value for keep.inbag")
  }
  
  ## Num threads
  ## Default 0 -> detect from system in C++.
  if (is.null(num.threads)) {
    num.threads = 0
  } else if (!is.numeric(num.threads) || num.threads < 0) {
    stop("Error: Invalid value for num.threads")
  }
  
  ## Minumum node size
  if (is.null(min.node.size)) {
    min.node.size <- 0
  } else if (!is.numeric(min.node.size) || min.node.size < 0) {
    stop("Error: Invalid value for min.node.size")
  }
  
  ## Tree depth
  if (is.null(max.depth)) {
    max.depth <- 0
  } else if (!is.numeric(max.depth) || max.depth < 0) {
    stop("Error: Invalid value for max.depth. Please give a positive integer.")
  }
  
  ## Sample fraction
  if (!is.numeric(sample.fraction)) {
    stop("Error: Invalid value for sample.fraction. Please give a value in (0,1] or a vector of values in [0,1].")
  }
  if (length(sample.fraction) > 1) {
    if (!(treetype %in% c(1, 9))) {
      stop("Error: Invalid value for sample.fraction. Vector values only valid for classification forests.")
    }
    if (any(sample.fraction < 0) || any(sample.fraction > 1)) {
      stop("Error: Invalid value for sample.fraction. Please give a value in (0,1] or a vector of values in [0,1].")
    }
    if (sum(sample.fraction) <= 0) {
      stop("Error: Invalid value for sample.fraction. Sum of values must be >0.")
    }
    if (length(sample.fraction) != nlevels(y)) {
      stop("Error: Invalid value for sample.fraction. Expecting ", nlevels(y), " values, provided ", length(sample.fraction), ".")
    }
    if (!replace & any(sample.fraction * length(y) > table(y))) {
      idx <- which(sample.fraction * length(y) > table(y))[1]
      stop("Error: Not enough samples in class ", names(idx), 
           "; available: ", table(y)[idx], 
           ", requested: ", (sample.fraction * length(y))[idx], ".")
    }
    if (!is.null(case.weights)) {
      stop("Error: Combination of case.weights and class-wise sampling not supported.")
    }
    # Fix order (C++ needs sample.fraction in order as classes appear in data)
    sample.fraction <- sample.fraction[as.numeric(unique(y))]
  } else {
    if (sample.fraction <= 0 || sample.fraction > 1) {
      stop("Error: Invalid value for sample.fraction. Please give a value in (0,1] or a vector of values in [0,1].")
    }
  }
  
  # Regularization
  if (all(regularization.factor == 1)) {
    regularization.factor <-  c(0, 0)
    use.regularization.factor <- FALSE
  } else {
    # Deactivation of paralellization
    if (num.threads != 1) {
      num.threads <- 1
      warning("Paralellization deactivated (regularization used).")
    }
    use.regularization.factor <- TRUE
  } 
  
  if (use.regularization.factor) {
    # A few checkings on the regularization coefficients
    if (max(regularization.factor) > 1) {
      stop("The regularization coefficients cannot be greater than 1.")
    }
    if (max(regularization.factor) <= 0) {
      stop("The regularization coefficients cannot be smaller than 0.")
    }
    p <- length(all.independent.variable.names)
    if (length(regularization.factor) != 1 && length(regularization.factor) != p) {
      stop("You must use 1 or p (the number of predictor variables)
      regularization coefficients.")
    }
    if (length(regularization.factor) == 1) {
      regularization.factor = rep(regularization.factor, p)
    }
  }
  
  ## Importance mode
  if (is.null(importance) || importance == "none") {
    importance.mode <- 0
  } else if (importance == "impurity") {
    importance.mode <- 1
  } else if (importance == "impurity_corrected" || importance == "impurity_unbiased") {
    importance.mode <- 5
  } else if (importance == "permutation") {
    if (local.importance) {
      importance.mode <- 6
    } else if (scale.permutation.importance) {
      importance.mode <- 2
    } else {
      importance.mode <- 3
    }
  } else {
    stop("Error: Unknown importance mode.")
  }
  
  ## Case weights: NULL for no weights or all weights equal
  if (is.null(case.weights) || length(unique(case.weights)) == 1) {
    case.weights <- c(0,0)
    use.case.weights <- FALSE
    if (holdout) {
      stop("Error: Case weights required to use holdout mode.")
    }
  } else {
    use.case.weights <- TRUE
    
    ## Sample from non-zero weights in holdout mode
    if (holdout) {
      sample.fraction <- sample.fraction * mean(case.weights > 0)
    }
    
    if (!replace && sum(case.weights > 0) < sample.fraction * nrow(x)) {
      stop("Error: Fewer non-zero case weights than observations to sample.")
    }
  }
  
  ## Manual inbag selection
  if (is.null(inbag)) {
    inbag <- list(c(0,0))
    use.inbag <- FALSE
  } else if (is.list(inbag)) {
    use.inbag <- TRUE
    if (use.case.weights) {
      stop("Error: Combination of case.weights and inbag not supported.")
    }
    if (length(sample.fraction) > 1) {
      stop("Error: Combination of class-wise sampling and inbag not supported.")
    }
    if (length(inbag) != num.trees) {
      stop("Error: Size of inbag list not equal to number of trees.")
    }
  } else {
    stop("Error: Invalid inbag, expects list of vectors of size num.trees.")
  }
  
  ## Class weights: NULL for no weights (all 1)
  if (is.null(class.weights)) {
    class.weights <- rep(1, nlevels(y))
  } else {
    if (!(treetype %in% c(1, 9))) {
      stop("Error: Argument class.weights only valid for classification forests.")
    }
    if (!is.numeric(class.weights) || any(class.weights < 0)) {
      stop("Error: Invalid value for class.weights. Please give a vector of non-negative values.")
    }
    if (length(class.weights) != nlevels(y)) {
      stop("Error: Number of class weights not equal to number of classes.")
    }

    ## Reorder (C++ expects order as appearing in the data)
    class.weights <- class.weights[unique(as.numeric(y))]
  }
  
  ## Split select weights: NULL for no weights
  if (is.null(split.select.weights)) {
    split.select.weights <- list(c(0,0))
    use.split.select.weights <- FALSE
  } else if (is.numeric(split.select.weights)) {
    if (length(split.select.weights) != length(all.independent.variable.names)) {
      stop("Error: Number of split select weights not equal to number of independent variables.")
    }
    split.select.weights <- list(split.select.weights)
    use.split.select.weights <- TRUE
  } else if (is.list(split.select.weights)) {
    if (length(split.select.weights) != num.trees) {
      stop("Error: Size of split select weights list not equal to number of trees.")
    }
    use.split.select.weights <- TRUE
  } else {
    stop("Error: Invalid split select weights.")
  }
  
  ## Always split variables: NULL for no variables
  if (is.null(always.split.variables)) {
    always.split.variables <- c("0", "0")
    use.always.split.variables <- FALSE
  } else {
    use.always.split.variables <- TRUE
  }
  
  if (use.split.select.weights && use.always.split.variables) {
    stop("Error: Please use only one option of split.select.weights and always.split.variables.")
  }
  
  ## Splitting rule
  if (is.null(splitrule)) {
    if (treetype == 5) {
      splitrule <- "logrank"
    } else if (treetype == 3) {
      splitrule <- "variance"
    } else if (treetype %in% c(1, 9)) {
      splitrule <- "gini"
    }
    splitrule.num <- 1
  } else if (splitrule == "logrank") {
    if (treetype == 5) {
      splitrule.num <- 1
    } else {
      stop("Error: logrank splitrule applicable to survival data only.")
    }
  } else if (splitrule == "gini") {
    if (treetype %in% c(1, 9)) {
      splitrule.num <- 1
    } else {
      stop("Error: Gini splitrule applicable to classification data only.")
    }
  } else if (splitrule == "variance") {
    if (treetype == 3) {
      splitrule.num <- 1
    } else {
      stop("Error: variance splitrule applicable to regression data only.")
    }
  } else if (splitrule == "auc" || splitrule == "C") {
    if (treetype == 5) {
      splitrule.num <- 2
    } else {
      stop("Error: C index splitrule applicable to survival data only.")
    }
  } else if (splitrule == "auc_ignore_ties" || splitrule == "C_ignore_ties") {
    if (treetype == 5) {
      splitrule.num <- 3
    } else {
      stop("Error: C index splitrule applicable to survival data only.")
    }
  } else if (splitrule == "maxstat") {
    if (treetype == 5 || treetype == 3) {
      splitrule.num <- 4
    } else {
      stop("Error: maxstat splitrule applicable to regression or survival data only.")
    }
  } else if (splitrule == "extratrees") {
    splitrule.num <- 5
  } else if (splitrule == "beta") {
    if (treetype == 3) {
      splitrule.num <- 6
    } else {
      stop("Error: beta splitrule applicable to regression data only.")
    }
    
    ## Check for 0..1 outcome
    if (min(y) < 0 || max(y) > 1) {
      stop("Error: beta splitrule applicable to regression data with outcome between 0 and 1 only.")
    }
  } else if (splitrule == "hellinger") {
    if (treetype %in% c(1, 9)) {
      splitrule.num <- 7
    } else {
      stop("Error: Hellinger splitrule only implemented for binary classification.")
    }
    if ((is.factor(y) && nlevels(y) > 2) || (length(unique(y)) > 2)) {
      stop("Error: Hellinger splitrule only implemented for binary classification.")
    }  
  } else if (splitrule == "poisson") {
    if (treetype == 3) {
      splitrule.num <- 8
    } else {
      stop("Error: poisson splitrule applicable to regression data only.")
    }
    
    ## Check for valid responses
    if (min(y) < 0 || sum(y) <= 0) {
      stop("Error: poisson splitrule applicable to regression data with non-positive outcome (y>=0 and sum(y)>0) only.")
    }
  } else {
    stop("Error: Unknown splitrule.")
  }
  
  ## Maxstat splitting
  if (alpha < 0 || alpha > 1) {
    stop("Error: Invalid value for alpha, please give a value between 0 and 1.")
  }
  if (minprop < 0 || minprop > 0.5) {
    stop("Error: Invalid value for minprop, please give a value between 0 and 0.5.")
  }
  if (splitrule == "maxstat" & use.regularization.factor) {
    stop("Error: Regularization cannot be used with 'maxstat' splitrule.")
  }

  ## Extra trees
  if (!is.numeric(num.random.splits) || num.random.splits < 1) {
    stop("Error: Invalid value for num.random.splits, please give a positive integer.")
  }
  if (splitrule.num == 5 && save.memory && respect.unordered.factors == "partition") {
    stop("Error: save.memory option not possible in extraTrees mode with unordered predictors.")
  }
  if (num.random.splits > 1 && splitrule.num != 5) {
    warning("Argument 'num.random.splits' ignored if splitrule is not 'extratrees'.")
  }
  
  if (!is.numeric(poisson.tau) || poisson.tau <= 0) {
    stop("Error: Invalid value for poisson.tau, please give a positive number.")
  }

  ## Unordered factors  
  if (respect.unordered.factors == "partition") {
    ordered.idx <- sapply(x, is.ordered)
    factor.idx <- sapply(x, is.factor)
    unordered.factor.variables <- independent.variable.names[factor.idx & !ordered.idx]
    
    if (length(unordered.factor.variables) > 0) {
      use.unordered.factor.variables <- TRUE
      ## Check level count
      num.levels <- sapply(x[, factor.idx & !ordered.idx, drop = FALSE], nlevels)
      max.level.count <- .Machine$double.digits
      if (max(num.levels) > max.level.count) {
        stop(paste("Too many levels in unordered categorical variable ", unordered.factor.variables[which.max(num.levels)], 
                   ". Only ", max.level.count, " levels allowed on this system. Consider using the 'order' option.", sep = ""))
      } 
    } else {
      unordered.factor.variables <- c("0", "0")
      use.unordered.factor.variables <- FALSE
    } 
  } else if (respect.unordered.factors == "ignore" || respect.unordered.factors == "order") {
    ## Ordering for "order" is handled above
    unordered.factor.variables <- c("0", "0")
    use.unordered.factor.variables <- FALSE
  } else {
    stop("Error: Invalid value for respect.unordered.factors, please use 'order', 'partition' or 'ignore'.")
  }

  ## Unordered maxstat splitting not possible
  if (use.unordered.factor.variables && !is.null(splitrule)) {
    if (splitrule == "maxstat") {
      stop("Error: Unordered factor splitting not implemented for 'maxstat' splitting rule.")
    } else if (splitrule %in% c("C", "auc", "C_ignore_ties", "auc_ignore_ties")) {
      stop("Error: Unordered factor splitting not implemented for 'C' splitting rule.")
    } else if (splitrule == "beta") {
      stop("Error: Unordered factor splitting not implemented for 'beta' splitting rule.")
    } else if (splitrule == "poisson") {
      stop("Error: Unordered factor splitting not implemented for 'poisson' splitting rule.")
    }
  }
  
  ## Warning for experimental 'order' splitting 
  if (respect.unordered.factors == "order") {
    if (treetype == 3 && splitrule == "maxstat") {
      warning("Warning: The 'order' mode for unordered factor handling with the 'maxstat' splitrule is experimental.")
    }
    if (gwa.mode & ((treetype %in% c(1,9) & nlevels(y) > 2) | treetype == 5)) {
      stop("Error: Ordering of SNPs currently only implemented for regression and binary outcomes.")
    }
  }
  
  ## Prediction mode always false. Use predict.ranger() method.
  prediction.mode <- FALSE
  predict.all <- FALSE
  prediction.type <- 1
  
  ## No loaded forest object
  loaded.forest <- list()
  
  ## Use sparse matrix
  if (inherits(x, "dgCMatrix")) {
    sparse.x <- x
    x <- matrix(c(0, 0))
    use.sparse.data <- TRUE
  } else {
    sparse.x <- Matrix(matrix(c(0, 0)))
    use.sparse.data <- FALSE
    if (is.data.frame(x)) {
      x <- data.matrix(x)
    }
  }
  
  if (treetype == 5) {
    y.mat <- as.matrix(y)
  } else {
    y.mat <- as.matrix(as.numeric(y))
  }
  
  if (respect.unordered.factors == "order"){
    order.snps <- TRUE
  } else {
    order.snps <- FALSE
  }
  
  ## No competing risks check
  if (treetype == 5) {
    if (!all(y.mat[, 2] %in% 0:1)) {
      stop("Error: Competing risks not supported yet. Use status=1 for events and status=0 for censoring.")
    }
  }
  
  ## Call Ranger
  result <- rangerCpp(treetype, x, y.mat, independent.variable.names, mtry,
                      num.trees, verbose, seed, num.threads, write.forest, importance.mode,
                      min.node.size, split.select.weights, use.split.select.weights,
                      always.split.variables, use.always.split.variables,
                      prediction.mode, loaded.forest, snp.data,
                      replace, probability, unordered.factor.variables, use.unordered.factor.variables, 
                      save.memory, splitrule.num, case.weights, use.case.weights, class.weights, 
                      predict.all, keep.inbag, sample.fraction, alpha, minprop, poisson.tau,
                      holdout, prediction.type, num.random.splits, sparse.x, use.sparse.data,
                      order.snps, oob.error, max.depth, inbag, use.inbag, 
                      regularization.factor, use.regularization.factor, regularization.usedepth)
  
  if (length(result) == 0) {
    stop("User interrupt or internal error.")
  }
  
  ## Prepare results
  if (importance.mode != 0) {
    names(result$variable.importance) <- all.independent.variable.names
    
    if (importance.mode == 6) {
      # process casewise vimp
      result$variable.importance.local <-
        matrix(
          result$variable.importance.local,
          byrow = FALSE,
          ncol = length(all.independent.variable.names),
          dimnames = list(
            rownames(data),
            all.independent.variable.names
          )
        )
    }
  }

  ## Set predictions
  if (treetype == 1 && oob.error) {
    if (is.factor(y)) {
      result$predictions <- integer.to.factor(result$predictions,
                                              levels(y))
    } 
    result$confusion.matrix <- table(y, result$predictions, 
                                     dnn = c("true", "predicted"), useNA = "ifany")
  } else if (treetype == 5 && oob.error) {
    if (is.list(result$predictions)) {
      result$predictions <- do.call(rbind, result$predictions)
    } 
    if (is.vector(result$predictions)) {
      result$predictions <- matrix(result$predictions, nrow = 1)
    }
    result$chf <- result$predictions
    result$predictions <- NULL
    result$survival <- exp(-result$chf)
  } else if (treetype == 9 && oob.error) {
    if (is.list(result$predictions)) {
      result$predictions <- do.call(rbind, result$predictions)
    } 
    if (is.vector(result$predictions)) {
      result$predictions <- matrix(result$predictions, nrow = 1)
    }
    
    ## Set colnames and sort by levels
    colnames(result$predictions) <- unique(y)
    if (is.factor(y)) {
      result$predictions <- result$predictions[, levels(droplevels(y)), drop = FALSE]
    }
  }
  
  ## Splitrule
  result$splitrule <- splitrule
  if (splitrule == "extratrees") {
    result$num.random.splits <- num.random.splits
  }
  
  ## Set treetype
  if (treetype == 1) {
    result$treetype <- "Classification"
  } else if (treetype == 3) {
    result$treetype <- "Regression"
  } else if (treetype == 5) {
    result$treetype <- "Survival"
  } else if (treetype == 9) {
    result$treetype <- "Probability estimation"
  }
  if (treetype == 3) {
    result$r.squared <- 1 - result$prediction.error / var(y)
  }
  result$call <- sys.call()
  result$importance.mode <- importance
  if (use.sparse.data) {
    result$num.samples <- nrow(sparse.x)
  } else {
    result$num.samples <- nrow(x)
  }
  result$replace <- replace
  
  ## Write forest object
  if (write.forest) {
    if (is.factor(y)) {
      result$forest$levels <- levels(y)
    }
    result$forest$independent.variable.names <- independent.variable.names
    result$forest$treetype <- result$treetype
    class(result$forest) <- "ranger.forest"
    
    ## In 'ordered' mode, save covariate levels
    if (respect.unordered.factors == "order" && ncol(x) > 0) {
      result$forest$covariate.levels <- covariate.levels
    }
  }
  
  class(result) <- "ranger"
  
  ## Prepare quantile prediction
  if (quantreg) {
    terminal.nodes <- predict(result, x, type = "terminalNodes")$predictions + 1
    n <- result$num.samples
    result$random.node.values <- matrix(nrow = max(terminal.nodes), ncol = num.trees)
    
    ## Select one random obs per node and tree
    for (tree in 1:num.trees){
      idx <- sample(1:n, n)
      result$random.node.values[terminal.nodes[idx, tree], tree] <- y[idx]
    }
    
    ## Prepare out-of-bag quantile regression
    if(!is.null(result$inbag.counts)) {
      inbag.counts <- simplify2array(result$inbag.counts)
      random.node.values.oob <- randomObsNode(terminal.nodes, y, inbag.counts)
      
      ## Check num.trees
      minoob <- min(rowSums(inbag.counts == 0))
      if (minoob < 10) {
        stop("Error: Too few trees for out-of-bag quantile regression.")
      }
      
      ## Use the same number of values for all obs, select randomly
      result$random.node.values.oob <- t(apply(random.node.values.oob, 1, function(x) {
        sample(x[!is.na(x)], minoob)
      }))
    }
  }
  
  return(result)
}

