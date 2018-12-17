#' Latent class model data generation
#'
#' Generate a data set using the simple latent class model as the underlying
#' data generating mechanism.
#'
#' Drops the last item if miss.prop = 1.
#'
#' @param n Sample size.
#' @param tau Prevalence of the disease, i.e. \eqn{\text{P}(\delta=1)}.
#' @param miss.prop Proportion of missing values in the gold standard item.
#' @param seed Random seed.
#' @param name.items Logical. Use item names as set in the options?
#'
#' @return Data frame.
#' @export
#'
#' @examples
#' gen_lc(n = 10)
gen_lc <- function(n = 250, tau = 0.08, miss.prop = 0.2, seed = NULL,
                   name.items = TRUE) {
  # Initialise and check -------------------------------------------------------
  read_diagacc_opt(environment())  # obtains sens, spec, item.names and p = length(sens)
  if (tau > 1 | tau < 0) stop("tau needs to be between 0 and 1")
  if (miss.prop > 1 | miss.prop < 0)
    stop("miss.prop needs to be between 0 and 1")

  # Generate data set ----------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  X <- data.frame(matrix(NA, nrow = n, ncol = p))
  d <- rbinom(n, size = 1, prob = tau)
  for (i in 1:n) {
    if (d[i] == 1) {
      X[i, ] <- rbinom(p, size = 1, prob = sens)
    } else {
      X[i, ] <- rbinom(p, size = 1, prob = 1 - spec)
    }
  }
  if (isTRUE(name.items)) colnames(X) <- item.names

  # Generate missing values in the gold standard -------------------------------
  pos.of.gs <- getOption("diagacc.gold")
  X[sample(seq(1, n, by = 1), n * miss.prop), pos.of.gs] <- NA

  X
}

#' Latent class with random effects model data generation
#'
#' Generate a data set using the latent class with random effects model as the
#' underlying data generating mechanism.
#'
#' Drops the last item if miss.prop = 1.
#'
#' @param n Sample size.
#' @param tau Prevalence of the disease, i.e. \eqn{\text{P}(\delta=1)}.
#' @param miss.prop Proportion of missing values in the gold standard item.
#' @param seed Random seed.
#' @param name.items Logical. Use item names as set in the options?
#' @param sigma The standard deviation of the random effects.
#'
#' @return Data frame.
#' @export
#'
#' @examples
#' gen_lcre(n = 10)
gen_lcre <- function(n = 250, tau = 0.08, miss.prop = 0.2, seed = NULL,
                     sigma = c(1.5, 1.5), name.items = TRUE) {
  # Initialise and check -------------------------------------------------------
  read_diagacc_opt(environment())  # obtains sens, spec, item.names and p = length(sens)
  if (tau > 1 | tau < 0) stop("tau needs to be between 0 and 1")
  if (miss.prop > 1 | miss.prop < 0)
    stop("miss.prop needs to be between 0 and 1")
  if (length(sigma) == 1) sigma <- c(sigma, sigma)

  beta <- matrix(NA, ncol = 2, nrow = p)
  for (j in 1:p) {
    beta[j, 1] <- qnorm(sens[j]) * sqrt(1 + sigma[1] ^ 2)
    beta[j, 2] <- qnorm(1 - spec[j]) * sqrt(1 + sigma[2] ^ 2)
  }

  # Generate data set ----------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  X <- data.frame(matrix(NA, nrow = n, ncol = p))
  d <- rbinom(n, size = 1, prob = tau)
  for (i in 1:n) {
    if (d[i] == 1) {
      X[i, ] <- rbinom(p, size = 1, prob = pnorm(beta[, 1] + sigma[1] * rnorm(1)))
    } else {
      X[i, ] <- rbinom(p, size = 1, prob = pnorm(beta[, 2] + sigma[2] * rnorm(1)))
    }
  }
  if (isTRUE(name.items)) colnames(X) <- item.names

  # Generate missing values in the gold standard -------------------------------
  pos.of.gs <- getOption("diagacc.gold")
  X[sample(seq(1, n, by = 1), n * miss.prop), pos.of.gs] <- NA

  X
}

#' Finite mixture model data generation
#'
#' Generate a data set using the finite mixture model as the
#' underlying data generating mechanism.
#'
#' Drops the last item if miss.prop = 1.
#'
#' @param n Sample size.
#' @param tau Prevalence of the disease, i.e. \eqn{\text{P}(\delta=1)}.
#' @param miss.prop Proportion of missing values in the gold standard item.
#' @param seed Random seed.
#' @param name.items Logical. Use item names as set in the options?
#' @param eta Eta.
#'
#' @return Data frame.
#' @export
#'
#' @author Elena Erosheva
#' @examples
#' gen_fm(n = 10)
#'
gen_fm <- function(n = 250, tau = 0.08, miss.prop = 0.2, seed = NULL,
                     eta = c(0.5, 0.2), name.items = TRUE) {
  # Initialise and check -------------------------------------------------------
  read_diagacc_opt(environment())  # obtains sens, spec, item.names and p = length(sens)
  if (tau > 1 | tau < 0) stop("tau needs to be between 0 and 1")
  if (miss.prop > 1 | miss.prop < 0)
    stop("miss.prop needs to be between 0 and 1")
  if (length(eta) == 1) eta <- c(eta, eta)

  N <- n  # sample size
  J <- p  # number of items/tests
  tau <- c(1 - tau, tau) # prevalence

  w <- matrix(rep(0, J * 2), ncol = 2)
  w[, 2] <- (sens - eta[1]) / (1 - eta[1])
  w[, 1] <- (spec - eta[2]) / (1 - eta[2])

  # Generate data set ----------------------------------------------------------
  resp <- data.frame(matrix(NA, nrow = n, ncol = p))
  if (!is.null(seed)) set.seed(seed)
  d <- rbinom(N, 1, tau[2])
  l0 <- rbinom(N, 1, eta[2])
  l1 <- rbinom(N, 1, eta[1])

  for (i in 1:N) {
    if (d[i] == 1) {
      if (l1[i] == 1) {
        for (j in 1:J) {
          resp[i, j] <- 1
        }
      } else {
        for (j in 1:J) {
          resp[i,j] <- rbinom(1, 1, w[j,2])
        }
      }
    }
    if (d[i] == 0) {
      if (l0[i] == 1) {
        for (j in 1:J) {
          resp[i, j] <- 0
        }
      } else {
        for (j in 1:J) {
          resp[i,j] <- rbinom(1, 1, 1 - w[j,1])
        }
      }
    }
  }
  X <- resp
  if (isTRUE(name.items)) colnames(X) <- item.names

  # Generate missing values in the gold standard -------------------------------
  pos.of.gs <- getOption("diagacc.gold")
  X[sample(seq(1, n, by = 1), n * miss.prop), pos.of.gs] <- NA

  X
}
