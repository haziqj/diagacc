#' @export
fit_lcre <- function(X, method = c("EM", "MCMC"), n.sample = 2000, n.chains = 1,
                     n.thin = 1,  n.burnin = 800, n.adapt = 200, raw = FALSE,
                     runjags.method = "rjags", silent = FALSE, quad.points = 21,
                     calcSE = TRUE) {
  method <- match.arg(method, c("EM", "MCMC"))

  res <- NULL
  if (method == "EM") {
    res <- fit_lcre_randomLCA(X = X, quad.points = quad.points, calcSE = calcSE)
  }
  if (method == "MCMC") {
    res <- fit_lcre_mcmc(X = X, n.chains = n.chains, n.sample = n.sample,
                         n.thin = n.thin, n.burnin = n.burnin, n.adapt = n.adapt,
                         runjags.method = runjags.method, silent = silent)
  }

  res$diagaccmod <- "LCRE"
  if (isTRUE(raw)) {
    return(res)
  } else {
    convert_mod_diagacc(res, X)
  }
}

fit_lcre_randomLCA <- function(X, quad.points = 21, calcSE = TRUE) {
  randomLCA::randomLCA(X, nclass = 2, probit = TRUE, random = TRUE,
                       calcSE = calcSE, constload = TRUE, byclass = TRUE,
                       quadpoints = quad.points)
}

fit_lcre_mcmc <- function(X, n.sample = 2000, n.chains = 1, n.thin = 1,
                          n.burnin = 800, n.adapt = 200,
                          runjags.method = "rjags", silent = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  X <- as.matrix(X)

  inits <- list(
    tau = 0.1,
    d = rep(0, n),
    beta = matrix(c(1, -1), nrow = p, ncol = 2, byrow = TRUE)
  )

  mod.jags <- "model{
    for (i in 1:n) {
      d[i] ~ dbern(tau)
      for (j in 1:p) {
        u1[i, j] ~ dnorm(beta[j, 1], psi[1])
        u0[i, j] ~ dnorm(beta[j, 2], psi[2])
        pi.x[i, j] <- d[i] * phi(u1[i, j]) + (1 - d[i]) * phi(u0[i, j])
        X[i, j] ~ dbern(pi.x[i, j])
      }
    }

    # Priors
    tau ~ dbeta(1,1)
    for (j in 1:p) {
      beta[j, 1] ~ dnorm(0, 1)
      beta[j, 2] ~ dnorm(0, 1)
    }
    for (k in 1:2) {
      psi[k] ~ dgamma(1, 1)
      sigma[k] <- 1 / sqrt(psi[k])
    }

    # Sensitivities and specificities
    for (j  in 1:p) {
      sens[j] <- phi(beta[j, 1] / sqrt(1 + pow(sigma[1], 2)))
      spec[j] <- 1 - phi(beta[j, 2] / sqrt(1 + pow(sigma[2], 2)))
    }
  }

  #data# X, n, p
  #monitor# tau, sens, spec, beta
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
