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

# Convert integer to factor
integer.to.factor <- function(x, labels) {
  factor(x, levels = seq_along(labels), labels = labels)
}

# Save version of sample() for length(x) == 1
# See help(sample)
save.sample <- function(x, ...) {
  x[sample.int(length(x), ...)]
}

# Order factor levels with PCA approach 
# Reference: Coppersmith, D., Hong, S.J. & Hosking, J.R. (1999) Partitioning Nominal Attributes in Decision Trees. Data Min Knowl Discov 3:197. \url{https://doi.org/10.1023/A:1009869804967}.
pca.order <- function(y, x) {
  x <- droplevels(x)
  
  ## Create contingency table of the nominal outcome with the nominal covariate
  N <- table(droplevels(y), x)
  
  ## PCA of class probabilites
  P <- N/rowSums(N)
  pc1 <- prcomp(P, rank. = 1)$rotation
  
  ## Return ordered factor levels
  as.character(levels(x)[order(pc1)])
}

# Compute median survival if available or largest quantile available in all strata if median not available.
largest.quantile <- function(formula) {
  ## Fit survival model
  fit <- survfit(formula)
  smry <- summary(fit)
  
  ## Use median survival if available or largest quantile available in all strata if median not available
  max_quant <- max(aggregate(smry$surv ~ smry$strata, FUN = min)[, "smry$surv"])
  quantiles <- quantile(fit, conf.int = FALSE, prob = min(0.5, 1 - max_quant))[, 1]
  names(quantiles) <- gsub(".+=", "", names(quantiles))
  
  ## Return ordered levels
  names(sort(quantiles))
}
