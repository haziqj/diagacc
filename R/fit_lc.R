#' @export
fit_lc <- function(X, method = c("EM", "MCMC"), n.sample = 2000, n.chains = 1,
                   n.thin = 1,  n.burnin = 800, n.adapt = 200, raw = FALSE,
                   runjags.method = "rjags", silent = FALSE, calcSE = TRUE) {
  method <- match.arg(method, c("EM", "MCMC"))

  res <- NULL
  if (method == "EM") {
    res <- fit_lc_randomLCA(X = X, calcSE = calcSE)
  }
  if (method == "MCMC") {
    res <- fit_lc_mcmc(X = X, n.chains = n.chains, n.sample = n.sample,
                       n.thin = n.thin, n.burnin = n.burnin, n.adapt = n.adapt,
                       runjags.method = runjags.method, silent = silent)
  }

  res$diagaccmod <- "LC"
  if (isTRUE(raw)) {
    return(res)
  } else {
    convert_mod_diagacc(res, X)
  }
}


fit_lc_randomLCA <- function(X, calcSE = TRUE) {
  randomLCA::randomLCA(X, nclass = 2, probit = TRUE, calcSE = calcSE)
}

fit_lc_mcmc <- function(X, n.sample = 2000, n.chains = 1, n.thin = 1,
                        n.burnin = 800, n.adapt = 200, runjags.method = "rjags",
                        silent = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  X <- as.matrix(X)

  inits <- list(
    tau = 0.1,
    d = rep(0, n),
    pi = matrix(c(0.9, 0.1), nrow = p, ncol = 2, byrow = TRUE)
  )

  mod.jags <- "model{
    for (i in 1:n) {
      d[i] ~ dbern(tau)
      for (j in 1:p) {
        pi.x[i, j] <- d[i] * pi[j, 1] + (1 - d[i]) * pi[j, 2]
        X[i, j] ~ dbern(pi.x[i, j])
      }
    }

    # Priors
    tau ~ dbeta(1, 1)
    for (j in 1:p) {
      pi[j, 1] ~ dbeta(1, 1)
      pi[j, 2] ~ dbeta(1, 1)
    }

    # Sensitivities and specificities
    for (j  in 1:p) {
      sens[j] <- pi[j, 1]
      spec[j] <- 1 - pi[j, 2]
    }
  }

  #data# X, n, p
  #monitor# tau, sens, spec
  "

  if (isTRUE(silent)) {
    runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE,
                             modules = "lecuyer")
  } else {
    runjags::runjags.options(silent.jags = FALSE, silent.runjags = FALSE,
                             modules = "lecuyer")
  }

  runjags::run.jags(mod.jags, n.chains = n.chains, sample = n.sample,
                    thin = n.thin, inits = inits, burnin = n.burnin,
                    adapt = n.adapt, method = runjags.method)
}
