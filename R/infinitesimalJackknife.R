
# Compute variance of estimate from a ranger model
# 
# Computes variances for a prediction from a ranger model, using the infinitesimal jackknife procedure
# 
# This function is a ranger-adaptation of the package \pkg{randomForestCI} of Wager et al. (2014). Their original can be found on github: \url{ https://github.com/swager/randomForestCI/}. 
#
# @param pred A nrow(newdata) by no. of trees matrix which contains numeric predictions
#        from a random forest trained with trees grown on bootstrap samples of the training data
# @param inbag A number of obs. in the training data by no. of trees matrix giving the
#        number of times the ith observation in the training data appeared in the bootstrap sample for the jth tree.
# @param calibrate whether to apply calibration to mitigate Monte Carlo noise
#        warning: if calibrate = FALSE, some variance estimates may be negative
#                 due to Monte Carlo effects if the number of trees in rf is too small
# @param used.trees set of trees to use for variance estimation; uses all tress if NULL
#
# @return A two-column matrix is returned, with predictions in the first column and estimates of prediction variance in the second. 
# @author Stefan Wager
rInfJack = function(pred, inbag, calibrate = TRUE, used.trees = NULL) {
  # original: https://github.com/swager/randomForestCI/blob/master/R/infinitesimalJackknife.R
  
  if (is.null(used.trees)) {
    used.trees = 1:ncol(inbag)
  }
  pred = pred[, used.trees, drop=FALSE]
  
  # Check if sampling without replacement
  no.replacement = (max(inbag) == 1)
  
  # Extract tree-wise predictions and variable counts from random forest
  B = length(used.trees)
  n = nrow(inbag)
  s = sum(inbag) / ncol(inbag)
  
  y.hat = rowMeans(pred)
  pred.centered = pred - rowMeans(pred)
  
  N = Matrix::Matrix(inbag[, used.trees], sparse = TRUE)
  N.avg = Matrix::rowMeans(N)
  
  # Compute raw infinitesimal jackknife
  if (as.numeric(B)^2 > as.numeric(n) * as.numeric(nrow(pred))) {
    
    C = Matrix::tcrossprod(N, pred.centered) -
      Matrix::Matrix(N.avg, nrow(N), 1) %*%
      Matrix::Matrix(rowSums(pred.centered), 1, nrow(pred.centered))
    raw.IJ = Matrix::colSums(C^2) / B^2
    
  } else {
    # Faster implementation when n is large. Uses the fact that
    # colSums((A - B)^2) = T1 - 2 * T2 + T3,
    # where T1 = diag(A'A), T2 = diag(B'A), and T3 = diag(B'B)
    
    NTN = Matrix::crossprod(N, N)
    NTNPT_T = Matrix::tcrossprod(pred.centered, NTN)
    T1 = Matrix::rowSums(pred.centered * NTNPT_T)
    
    RS = rowSums(pred.centered)
    NbarTN = Matrix::crossprod(N.avg, N)
    T2 = RS * Matrix::tcrossprod(NbarTN, pred.centered)
    
    T3 = sum(N.avg^2) * RS^2
    raw.IJ = as.numeric(T1 - 2 * T2 + T3) / B^2
    
  }
  
  # Apply Monte Carlo bias correction
  N.var = mean(Matrix::rowMeans(N^2) - Matrix::rowMeans(N)^2)
  boot.var = rowSums(pred.centered^2) / B
  bias.correction = n * N.var * boot.var / B
  vars = raw.IJ - bias.correction
  
  # Finite sample correction
  if (no.replacement) {
    variance.inflation = 1 / (1 - mean(inbag))^2
    vars = variance.inflation * vars
  }
  
  results = data.frame(y.hat=y.hat, var.hat=vars)
  
  if (isTRUE(calibrate) && nrow(results) <= 20) {
    calibrate = FALSE
    warning("Sample size <=20, no calibration performed.")
  }
  
  # If appropriate, calibrate variance estimates; this step in particular
  # ensures that all variance estimates wil be positive.
  
  if (calibrate) {
    # Compute variance estimates using half the trees
    calibration.ratio = 2
    n.sample = ceiling(B / calibration.ratio)
    results.ss = rInfJack(pred, inbag, calibrate = FALSE, used.trees = sample(used.trees, n.sample))
    
    # Use this second set of variance estimates to estimate scale of Monte Carlo noise
    sigma2.ss = mean((results.ss$var.hat - results$var.hat)^2)
    delta = n.sample / B
    sigma2 = (delta^2 + (1 - delta)^2) / (2 * (1 - delta)^2) * sigma2.ss
    
    # Use Monte Carlo noise scale estimate for empirical Bayes calibration
    results = tryCatch(
      expr = {
        vars.calibrated = calibrateEB(vars, sigma2)
        results$var.hat = vars.calibrated
        results
      }, 
      error = function(e) {
        warning(sprintf("Calibration failed with error:\n%sFalling back to non-calibrated variance estimates.", e))
        results = rInfJack(pred, inbag, calibrate = FALSE, used.trees = used.trees)
        return(results)
      }
    )
  }
  
  return(results)
}

# Fit an empirical Bayes prior in the hierarchical model
#     mu ~ G, X ~ N(mu, sigma^2)
#
# @param X a vector of observations
# @param sigma noise estimate
# @param p tuning parameter -- number of parameters used to fit G
# @param nbin tuning parameter -- number of bins used for discrete approximation
# @param unif.fraction tuning parameter -- fraction of G modeled as "slab"
#
# @return posterior density estimate g
#
# @section References:
# For more details about "g-estimation", see: B Efron. Two modeling strategies for
# empirical Bayes estimation. Stat. Sci., 29: 285-301, 2014.
# @author Stefan Wager
gfit = function(X, sigma, p = 2, nbin = 1000, unif.fraction = 0.1) {
  
  xvals = seq(min(min(X) - 2 * sd(X), 0), max(max(X) + 2 * sd(X), sd(X)), length.out = nbin)
  binw = xvals[2] - xvals[1]
  
  zero.idx = max(which(xvals <= 0))
  noise.kernel = dnorm(xvals / sigma) * binw / sigma
  
  if (zero.idx > 1) {
    noise.rotate = noise.kernel[c(zero.idx:length(xvals), 1:(zero.idx - 1))]
  } else {
    noise.rotate = noise.kernel
  }
  
  XX = sapply(1:p, function(j) xvals^j * as.numeric(xvals >= 0))
  neg.loglik = function(eta) {
    g.eta.raw = exp(XX %*% eta) * as.numeric(xvals >= 0)
    if ((sum(g.eta.raw) == Inf) | (sum(g.eta.raw) <= 100 * .Machine$double.eps)) {
      return (1000 * (length(X) + sum(eta^2)))
    }
    g.eta.main = g.eta.raw / sum(g.eta.raw)
    g.eta = (1 - unif.fraction) * g.eta.main +
      unif.fraction * as.numeric(xvals >= 0) / sum(xvals >= 0)
    f.eta = convolve(g.eta, noise.rotate)
    sum(approx(xvals, -log(pmax(f.eta, 0.0000001)), X)$y)
  }
  
  eta.hat = nlm(neg.loglik, rep(-1, p))$estimate
  g.eta.raw = exp(XX %*% eta.hat) * as.numeric(xvals >= 0)
  g.eta.main = g.eta.raw / sum(g.eta.raw)
  g.eta = (1 - unif.fraction) * g.eta.main +
    unif.fraction * as.numeric(xvals >= 0) / sum(xvals >= 0)
  
  return(data.frame(x=xvals, g=g.eta))
}

# Bayes posterior estimation with Gaussian noise
#
# @param x0 an obsevation
# @param g.est a prior density, as returned by gfit
# @param sigma noise estimate
#
# @return posterior estimate E[mu | x0]
# @author Stefan Wager
gbayes = function(x0, g.est, sigma) {
  Kx = dnorm((g.est$x - x0) / sigma)
  post = Kx * g.est$g
  post = post / sum(post)
  sum(post * g.est$x)
}

# Empirical Bayes calibration of noisy variance estimates
#
# @param vars list of variance estimates
# @param sigma2 estimate of the Monte Carlo noise in vars
#
# @return calibrated variance estimates
# @author Stefan Wager
calibrateEB = function(vars, sigma2) {
  
  if(sigma2 <= 0 | min(vars) == max(vars)) {
    return(pmax(vars, 0))
  }
  
  sigma = sqrt(sigma2)
  eb.prior = gfit(vars, sigma)
  
  if (length(vars >= 200)) {
    # If there are many test points, use interpolation to speed up computations
    calib.x = unique(quantile(vars, q = seq(0, 1, by = 0.02)))
    calib.y = sapply(calib.x, function(xx) gbayes(xx, eb.prior, sigma))
    calib.all = approx(x=calib.x, y=calib.y, xout=vars)$y
  } else {
    calib.all = sapply(vars, function(xx) gbayes(xx, eb.prior, sigma))
  }
  return(calib.all)
}

