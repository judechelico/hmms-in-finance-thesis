library(MASS)
library(ellipse)
library(mclust)

genDataGMM = function(n, numsets){
  xs = list()
  ys = list()
  
  centers = list()
  sigs = list()
  
  for (i in 1:numsets){
    mu = rnorm(2, 0, 7)
    sig = crossprod(
      matrix(rnorm(4), 2)
    ) # positive definite 
    # covariance matrix
    
    centers[[i]] = mu
    sigs[[i]] = sig
    
    sample = mvrnorm(n[i], mu, sig)
    xs[[i]] = sample[, 1]
    ys[[i]] = sample[, 2]
  }
  
  comb = data.frame(X = unlist(xs), Y = unlist(ys))
  return(list('data' = comb, 'true mus' = centers,
              'true covs' = sigs, 'true pis' = n / sum(n)))
}

# Log-likelihood: X is N x d data.frame or matrix
logLikelihood <- function(X, pis, mus, covs) {
  Xmat <- as.matrix(X)
  N <- nrow(Xmat)
  k <- length(pis)
  ll <- 0
  for (i in 1:N) {
    xi <- Xmat[i, ]
    probs <- numeric(k)
    for (j in 1:k) {
      probs[j] <- pis[j] * dmvnorm(
        t(xi), mean = t(mus[j, ]), sigma = covs[[j]]
      )
    }
    # avoid log(0)
    s <- sum(probs)
    ll <- ll + log(s + .Machine$double.eps)
  }
  return(ll)
}

# E-step: compute responsibilities (N x k)
eStep <- function(X, pis, mus, covs) {
  Xmat <- as.matrix(X)
  N <- nrow(Xmat)
  k <- length(pis)
  gamma <- matrix(0, nrow = N, ncol = k)
  for (i in 1:N) {
    xi <- Xmat[i, ]
    numerators <- numeric(k)
    for (j in 1:k) {
      numerators[j] <- pis[j] * dmvnorm(
        t(xi), mean = t(mus[j, ]), sigma = covs[[j]]
        )
    }
    denom <- sum(numerators)
    # normalize; add tiny eps to denom to avoid 
    # division-by-zero
    denom <- denom + .Machine$double.eps
    gamma[i, ] <- numerators / denom
  }
  return(gamma)  # returns N x k matrix
}

# M-step: compute parameters
mStep = function(Xmat, gamma, pis, mus, covs, N, k, d){
  reg = 1e-6
  Nks <- colSums(gamma)  # vector length k
  # update pis
  pis <- Nks / N
  # update mus
  for (j in 1:k) {
    # weighted mean
    mus[j, ] <- colSums(
      gamma[, j] * Xmat) / (Nks[j] + .Machine$double.eps)
  }
  # update covariances
  for (j in 1:k) {
    S <- matrix(0, nrow = d, ncol = d)
    mu_j <- mus[j, ]
    for (i in 1:N) {
      diff <- matrix(Xmat[i, ] - mu_j, ncol = 1)
      S <- S + gamma[i, j] * (diff %*% t(diff))
    }
    cov_j <- S / (Nks[j] + .Machine$double.eps)
    # regularize to ensure positive-definite
    covs[[j]] <- cov_j + reg * diag(1, d)
  }
  result = list("pis" = pis, "mus" = mus, "covs" = covs)
  return(result)
}

# EM algorithm
emAlg <- function(data, k, trueMus = NA, trueCovs = NA, truePis = NA,
                  maxIter = 200, tol = 1e-4, reg = 1e-6,
                  truths = FALSE) {
  Xmat <- as.matrix(data)
  N <- nrow(Xmat)
  d <- ncol(Xmat)
  
  # initialization uses k-means clustering
  km <- kmeans(Xmat, centers = k, nstart = 5)
  mus <- as.matrix(km$centers)
  pis <- as.numeric(table(km$cluster) / N)
  covs <- vector("list", length = k)
  for (i in 1:k) covs[[i]] <- diag(1, d)
  
  logLikeOld <- -Inf
  logLikeNew <- logLikelihood(data, pis, mus, covs)
  logLikes <- numeric()
  iter <- 1
  
  while (
    abs(logLikeNew - logLikeOld) > tol && iter <= maxIter
  ) {
    logLikeOld <- logLikeNew
    # E-step
    gamma <- eStep(data, pis, mus, covs)  # N x k
    
    # M-step
    maxed = mStep(Xmat, gamma, pis, mus, covs, N, k, d)
    pis = maxed[["pis"]]
    mus = maxed[["mus"]]
    covs = maxed[["covs"]]
    
    # compute log-likelihood
    logLikeNew <- logLikelihood(data, pis, mus, covs)
    logLikes <- c(logLikes, logLikeNew)
    iter <- iter + 1
  }
  
  if (iter > maxIter) message("Reached max iterations.")
  
  par(mfrow = c(1, 2))
  plot(Xmat, main = 'Model Fit',
       xlab = "Simulated x-values",
       ylab = "Simulated y-values")
  points(mus, pch = 20, col = 'purple')
  if (truths){
      for (i in 1:length(trueMus)){
      points(
        trueMus[[i]][2] ~ trueMus[[i]][1], pch = 18, col = 'red'
      )
      }
    }
  for (i in 1:k){
    points(ellipse(covs[[i]], center = mus[i,]), 
           type = 'l', col = 'purple', lwd = 2)
    if (truths){
      points(ellipse(trueCovs[[i]], center = trueMus[[i]]), 
             type = 'l', col = 'red', lty = 2, lwd = 2)
    }
  }
  plot(logLikes, type = 'l', 
       main = 'Log Likelihood over Iterations', 
       xlab = 'Iteration', ylab = 'Log-Likelihood')
  
  # return a clean list
  return(list(
    muhats = mus, covhats = covs, pihats = pis, 
    mus = trueMus, covs = trueCovs, pis = truePis,
    logLikes = logLikes
  ))
}

gmmWrapper = function(seed, k){
  set.seed(seed)
  # number of data points per cluster is sampled randomly
  truth = genDataGMM(sample(15:90, k), k)
  data = truth[['data']]
  mus = truth[['true mus']]
  covs = truth[['true covs']]
  pis = truth[['true pis']]
  emAlg(data, k, mus, covs, pis, truths = TRUE)
}

# model's prediction are the solid purple ellipses,
# truths are the dotted red ellipses
# gmmWrapper(8, 5)
