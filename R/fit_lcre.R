#' Fit a latent class with random effects model
#'
#' Fit a latent class with random effects model using MCMC.
#'
#' Priors for probabilities are Unif(0, 1). For the betas, we found that a
#' standard normal prior works best, because assigned a wider range for the
#' betas a priori disrupts the sensitivities and specificities calculations
#' (sticky chains at values close 1). Initial value for the prevalence is set
#' at 0.1, the disease indicators to zero for all units.
#'
#' Note that when \code{gold.std} is \code{TRUE}, then the last column in
#' \code{X} is assumed to be the gold standard item responses. Thus, the
#' sensitivities and specificities attached to this item is fixed to 1.
#'
#' @inheritParams fit_lc
#' @param quad.points (numeric, positive) Number of quadrature points for
#'   randomLCA fit. Check randomLCA documentation.
#'
#' @export
fit_lcre <- function(X, n.sample = 2000, n.chains = 2, n.thin = 1,
                     n.burnin = 800, n.adapt = 200, raw = FALSE,
                     runjags.method = "rjags", silent = FALSE, quad.points = 21,
                     calcSE = TRUE, gold.std = FALSE, method = c("MCMC", "EM")) {
  # method <- match.arg(method, c("MCMC", "EM"))
  method <- "MCMC"  # Deprecate the method option
  if (all(is.na(X[, ncol(X)]))) {
    gold.std <- FALSE
    X <- X[, -ncol(X)]
  }

  res <- NULL
  if (method == "EM") {
    # res <- fit_lcre_randomLCA(X = X, quad.points = quad.points, calcSE = calcSE)
    stop("Removed functionality from package. Use MCMC.")
  }
  if (method == "MCMC") {
    res <- fit_lcre_mcmc(X = X, n.chains = n.chains, n.sample = n.sample,
                         n.thin = n.thin, n.burnin = n.burnin, n.adapt = n.adapt,
                         runjags.method = runjags.method, silent = silent,
                         gold.std = gold.std)
  }

  res$diagaccmod <- "LCRE"
  if (isTRUE(raw)) {
    return(res)
  } else {
    convert_mod_diagacc(res, X)
  }
}

fit_lcre_randomLCA <- function(X, quad.points = 21, calcSE = TRUE) {
  # DEPRECATED. Helper function to fit LCRE models using EM algorithm. This uses
  # the randomLCA package.
  #
  # Args: X data matrix, quad.points determine the number of quadrature points
  # for the adaptive quadrature in the E-step, calcSE logical to calculate SE or
  # not.
  #
  # Returns: A randomLCA fit object.
  randomLCA::randomLCA(X, nclass = 2, probit = TRUE, random = TRUE,
                       calcSE = calcSE, constload = TRUE, byclass = TRUE,
                       quadpoints = quad.points)
}

fit_lcre_mcmc <- function(X, n.sample = 2000, n.chains = 1, n.thin = 1,
                          n.burnin = 800, n.adapt = 200,
                          runjags.method = "rjags", silent = FALSE,
                          gold.std = FALSE) {
  # Helper function to fit LCRE models using MCMC. This uses JAGS. There are two
  # versions of the MCMC model, one where the last item of X is the gold
  # standard (and therefore the sensitivities and specificities are fixed to 1),
  # and the other is where there is no gold standard available.
  #
  # Args: X data matrix, gold.std logical, and the rest are standard rjags
  # options.
  #
  # Returns: A runjags object.
  n <- nrow(X)
  p <- ncol(X)
  X <- as.matrix(X)

  inits <- list(
    tau = 0.1,
    d = rep(0, n)
  )
  
  if (isTRUE(gold.std)) {
    # This is the model for gold standard at the final column ------------------
    inits$beta <- matrix(c(rep(c(1, -1), p - 1), 100, -100), nrow = p,
                         ncol = 2, byrow = TRUE)
    mod.jags.lcre <- "model{
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
      for (j in 1:(p-1)) {
        beta[j, 1] ~ dnorm(0, 0.01)
        beta[j, 2] ~ dnorm(0, 0.01)
      }
      beta[p, 1] ~ dnorm(100, 100)   # This fixes sens and spec to 1.
      beta[p, 2] ~ dnorm(-100, 100)  #
      for (k in 1:2) {
        psi[k] ~ dgamma(0.1, 0.1)
        sigma[k] <- 1 / sqrt(psi[k])
      }

      # Sensitivities and specificities
      for (j  in 1:p) {
        sens[j] <- phi(beta[j, 1] / sqrt(1 + pow(sigma[1], 2)))
        spec[j] <- 1 - phi(beta[j, 2] / sqrt(1 + pow(sigma[2], 2)))
      }
    }

    #data# X, n, p
    #monitor# tau, sens, spec, beta, sigma, deviance
    #inits# tau, d, beta
    "
  } else {
    # This is the model for NO gold standard -----------------------------------
    inits$beta <- matrix(c(1, -1), nrow = p, ncol = 2, byrow = TRUE)
    mod.jags.lcre <- "model{
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
        beta[j, 1] ~ dnorm(0, 0.01)
        beta[j, 2] ~ dnorm(0, 0.01)
      }
      for (k in 1:2) {
        psi[k] ~ dgamma(0.1, 0.1)
        sigma[k] <- 1 / sqrt(psi[k])
      }

      # Sensitivities and specificities
      for (j  in 1:p) {
        sens[j] <- phi(beta[j, 1] / sqrt(1 + pow(sigma[1], 2)))
        spec[j] <- 1 - phi(beta[j, 2] / sqrt(1 + pow(sigma[2], 2)))
      }
    }

    #data# X, n, p
    #monitor# tau, sens, spec, beta, sigma, deviance
    #inits# tau, d, beta
    "
  }

  if (isTRUE(silent)) {
    runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE,
                             modules = "lecuyer")
  } else {
    runjags::runjags.options(silent.jags = FALSE, silent.runjags = FALSE,
                             modules = "lecuyer")
  }

  runjags::run.jags(mod.jags.lcre, n.chains = n.chains, sample = n.sample,
                    thin = n.thin, burnin = n.burnin, inits = inits,
                    adapt = n.adapt, method = runjags.method)
}
