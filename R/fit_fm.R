#' Fit a finite mixture model
#'
#' Fit a finite mixture model using MCMC.
#'
#' Uninformative priors are used, e.g. Unif(0, 1) for probabilities. Initial
#' value for the prevalence is set at 0.1, the disease indicators to zero for
#' all units, probabilities of correctly diagnosing patients (eta) to 0.1, and
#' probabilities of the tests correctly diagnosing patients when the patient
#' truly has or does not have the diseas as 0.9 and 0.7 respectively.
#'
#' Note that when \code{gold.std} is \code{TRUE}, then the last column in
#' \code{X} is assumed to be the gold standard item responses. Thus, the
#' sensitivities and specificities attached to this item is fixed to 1.
#'
#' @param X (matrix) Data set.
#' @param n.sample Number of MCMC samples.
#' @param n.chains Number of chains.
#' @param n.thin Thinning value.
#' @param n.burnin Number of burn-in.
#' @param n.adapt Number of adaptation samples.
#' @param raw (logical) Return the randomLCA or runjags object?
#' @param runjags.method Parallel or normal method. See runjags documentation.
#' @param silent (logical) Suppress output.
#'
#' @export
fit_fm <- function(X, n.sample = 2000, n.chains = 2, n.thin = 1, n.burnin = 800,
                   n.adapt = 200, raw = FALSE, runjags.method = "rjags",
                   silent = FALSE, gold.std = FALSE) {
  if (all(is.na(X[, ncol(X)]))) {
    gold.std <- FALSE
    X <- X[, -ncol(X)]
  }
  Y <- as.matrix(X)
  I <- nrow(Y)
  J <- ncol(Y)
  K <- 2

  # Initial values
    l    = rep(0, I),
    tau  = 0.1,
    eta  = c(0.1, 0.1)
  

  if (isTRUE(gold.std)) {
    # This is the model for gold standard at the final column ------------------
    w <- matrix(c(rep(c(0.9, 0.7), J - 1), 0.999999, 0.999999), nrow = J,
                      ncol = 2, byrow = TRUE)
    mod.jags.fm <- "model{
      # NOTE: class 1 = diseased, class 2 = healthy
      for (i in 1:I) {
        for (j in 1:J) {
          Y[i,j] ~ dbern(A[i,j])
          A[i,j] <- pow(w[j, 2] * d[i] + (1 - w[j,1]) * (1 - d[i]), 1 - l[i]) - l[i] * (1 - d[i])
        }
        l[i] ~ dbern(pi.l[i])
        pi.l[i] <- (1 - d[i]) * eta[1] + d[i] * eta[2]
        d[i] ~ dbern(tau)
      }

      # Priors and sens/spec
      tau ~ dbeta(1,1)
      for (k in 1:K) {
        eta[k] ~ dbeta(1,1)
        w[J,k] ~ dunif(0.99999, 1.0)
      }
      for (j in 1:(J-1)) {
        w[j,1] ~ dbeta(1,1)
        w[j,2] ~ dbeta(1,1)
        sens[j] <- eta[2] + (1 - eta[2]) * w[j,2]
        spec[j] <- eta[1] + (1 - eta[1]) * w[j,1]
      }
      sens[J] <- eta[2] + (1 - eta[2]) * w[J,2]
      spec[J] <- eta[1] + (1 - eta[1]) * w[J,1]
    }

    #data# I, J, K, Y
    #monitor# tau, sens, spec, w, eta, deviance
    #inits# l, tau, eta, w
    "
  } else {
    # This is the model for NO gold standard -----------------------------------
    w <- matrix(c(0.9, 0.7), nrow = J, ncol = 2, byrow = TRUE)
    mod.jags.fm <- "model{
      # NOTE: class 1 = diseased, class 2 = healthy
      for (i in 1:I) {
        for (j in 1:J) {
          Y[i,j] ~ dbern(A[i,j])
          A[i,j] <- pow(w[j, 2] * d[i] + (1 - w[j,1]) * (1 - d[i]), 1 - l[i]) - l[i] * (1 - d[i])
        }
        l[i] ~ dbern(pi.l[i])
        pi.l[i] <- (1 - d[i]) * eta[1] + d[i] * eta[2]
        d[i] ~ dbern(tau)
      }

      # Priors and sens/spec
      tau ~ dbeta(1,1)
      for (k in 1:K) {
        eta[k] ~ dbeta(1,1)
      }
      for (j in 1:J) {
        w[j,1] ~ dbeta(1,1)
        w[j,2] ~ dbeta(1,1)
        sens[j] <- eta[2] + (1 - eta[2]) * w[j,2]
        spec[j] <- eta[1] + (1 - eta[1]) * w[j,1]
      }
    }

    #data# I, J, K, Y
    #monitor# tau, sens, spec, w, eta, deviance
    #inits# l, tau, eta, w
    "
  }

  if (isTRUE(silent)) {
    runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE,
                             modules = "lecuyer")
  } else {
    runjags::runjags.options(silent.jags = FALSE, silent.runjags = FALSE,
                             modules = "lecuyer")
  }

  mod <- runjags::run.jags(mod.jags.fm, n.chains = n.chains, sample = n.sample,
                           thin = n.thin, burnin = n.burnin,
                           adapt = n.adapt, method = runjags.method)

  mod$diagaccmod <- "FM"

  if (isTRUE(raw)) {
    return(mod)
  } else {
    convert_mod_diagacc(mod, X)
  }
}
