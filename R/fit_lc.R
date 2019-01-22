#' Fit a latent class model
#'
#' Fit a latent class model using MCMC.
#'
#' Uninformative priors are used, e.g. Unif(0, 1) for probabilities. Initial
#' value for the prevalence is set at 0.1, the disease indicators to zero for
#' all units, and the sensitivities and specificities to 0.9.
#'
#' Note that when \code{gold.std} is \code{TRUE}, then the last column in
#' \code{X} is assumed to be the gold standard item responses. Thus, the
#' sensitivities and specificities attached to this item is fixed to 1.
#'
#' @param method (DEPRECATED--Use MCMC only) EM algorithm or MCMC.
#' @param calcSE (logical) Calculate standard error of estimates for randomLCA
#'   fit.
#' @param X (Matrix) Item responses data set.
#' @param n.sample Number of MCMC samples.
#' @param n.chains Number of chains.
#' @param n.thin Thinning value.
#' @param n.burnin Number of burn-in.
#' @param n.adapt Number of adaptation samples.
#' @param raw (Logical) Return the randomLCA or runjags object?
#' @param runjags.method Parallel or normal method. See runjags documentation.
#' @param silent (Logical) Suppress output.
#' @param gold.std (Logical) Is the last item/column in X the gold standard?
#'
#' @export
fit_lc <- function(X, n.sample = 2000, n.chains = 2, n.thin = 1, n.burnin = 800,
                   n.adapt = 200, raw = FALSE, runjags.method = "rjags",
                   silent = FALSE, calcSE = TRUE, gold.std = FALSE,
                   method = c("MCMC", "EM")) {
  # method <- match.arg(method, c("MCMC", "EM"))
  method <- "MCMC"  # Deprecate the method option
  if (all(is.na(X[, ncol(X)]))) {
    gold.std <- FALSE
    X <- X[, -ncol(X)]
  }

  res <- NULL
  if (method == "EM") {
    # res <- fit_lc_randomLCA(X = X, calcSE = calcSE)
    stop("Removed functionality from package. Use MCMC.")
  }
  if (method == "MCMC") {
    res <- fit_lc_mcmc(X = X, n.chains = n.chains, n.sample = n.sample,
                       n.thin = n.thin, n.burnin = n.burnin, n.adapt = n.adapt,
                       runjags.method = runjags.method, silent = silent,
                       gold.std = gold.std)
  }

  res$diagaccmod <- "LC"
  if (isTRUE(raw)) {
    return(res)
  } else {
    convert_mod_diagacc(res, X)
  }
}

fit_lc_randomLCA <- function(X, calcSE = TRUE) {
  # DEPRECATED. Helper function to fit LC models using EM algorithm. This uses
  # the randomLCA package.
  #
  # Args: X data matrix, calcSE logical to calculate SE or not.
  #
  # Returns: A randomLCA fit object.
  randomLCA::randomLCA(X, nclass = 2, probit = TRUE, calcSE = calcSE)
}

fit_lc_mcmc <- function(X, n.sample = 2000, n.chains = 1, n.thin = 1,
                        n.burnin = 800, n.adapt = 200, runjags.method = "rjags",
                        silent = FALSE, gold.std = FALSE) {
  # Helper function to fit LC models using MCMC. This uses JAGS. There are two
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
    inits$pi <- matrix(c(rep(c(0.9, 0.1), p - 1), 0.999999, 0.000001), nrow = p,
                       ncol = 2, byrow = TRUE)
    mod.jags.lc <- "model{
      for (i in 1:n) {
        d[i] ~ dbern(tau)
        for (j in 1:p) {
          pi.x[i, j] <- d[i] * pi[j, 1] + (1 - d[i]) * pi[j, 2]
          X[i, j] ~ dbern(pi.x[i, j])
        }
      }

      # Priors (uninformative)
      tau ~ dbeta(1, 1)
      for (j in 1:(p-1)) {
        pi[j, 1] ~ dbeta(1, 1)
        pi[j, 2] ~ dbeta(1, 1)
      }
      pi[p, 1] ~ dunif(0.99999, 1.00)  # This fixes sens and spec to 1. Initial
      pi[p, 2] ~ dunif(0.00, 0.00001)  # values must be within this range.

      # Sensitivities and specificities
      for (j  in 1:p) {
        sens[j] <- pi[j, 1]
        spec[j] <- 1 - pi[j, 2]
      }
    }

    #data# X, n, p
    #monitor# tau, sens, spec, deviance
    "
  } else {
    # This is the model for NO gold standard -----------------------------------
    inits$pi <- matrix(c(0.9, 0.1), nrow = p, ncol = 2, byrow = TRUE)
    mod.jags.lc <- "model{
      for (i in 1:n) {
        d[i] ~ dbern(tau)
        for (j in 1:p) {
          pi.x[i, j] <- d[i] * pi[j, 1] + (1 - d[i]) * pi[j, 2]
          X[i, j] ~ dbern(pi.x[i, j])
        }
      }
      # Priors (uninformative)
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
    #monitor# tau, sens, spec, deviance
    "
  }

  if (isTRUE(silent)) {
    runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE,
                             modules = "lecuyer")
  } else {
    runjags::runjags.options(silent.jags = FALSE, silent.runjags = FALSE,
                             modules = "lecuyer")
  }

  runjags::run.jags(mod.jags.lc, n.chains = n.chains, sample = n.sample,
                    thin = n.thin, inits = inits, burnin = n.burnin,
                    adapt = n.adapt, method = runjags.method)
}
