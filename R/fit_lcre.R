#' @export
fit_lcre <- function(X, method = c("EM", "MCMC"), n.chains = 1, n.sample = 10000,
                     n.thin = 1,  n.burnin = 800, n.adapt = 200, raw = FALSE,
                     runjags.method = "rjags", silent = FALSE, quad.points = 21) {
  method <- match.arg(method, c("EM", "MCMC"))

  res <- NULL
  if (method == "EM") {
    res <- fit_lcre_randomLCA(X = X, raw = raw, quad.points = quad.points)
  }
  if (method == "MCMC") {
    res <- fit_lcre_mcmc(X = X, n.chains = n.chains, n.sample = n.sample,
                         n.thin = n.thin, n.burnin = n.burnin, n.adapt = n.adapt,
                         raw = raw, runjags.method = runjags.method,
                         silent = silent)
  }
  res
}


fit_lcre_randomLCA <- function(X, raw = FALSE, quad.points = 21) {
  mod <- randomLCA::randomLCA(X, nclass = 2, probit = TRUE, random = TRUE,
                              calcSE = FALSE, constload = TRUE, byclass = TRUE,
                              quadpoints = quad.points)

  if (isTRUE(raw)) {
    return(mod)
  } else {
    class.probs <- randomLCA::classProbs(mod)  # class probs P[delta]
    probs <- mod$outcomep  # P[X = 1 | delta]

    delta <- which_class_diseased(probs)

    sens.fit <- probs[delta, ]  # sensitivity
    spec.fit <- 1 - probs[-delta, ]  # specificity
    sens.and.spec <- data.frame(sens.fit, spec.fit)
    colnames(sens.and.spec) <- c("Sensitivity", "Specificity")
    rownames(sens.and.spec) <- getOption("diagacc.item.names")

    # Remove gold standard -----------------------------------------------------
    pos.of.gs <- getOption("diagacc.gold")
    if (!is.na(pos.of.gs)) sens.and.spec <- sens.and.spec[-pos.of.gs, ]

    list(prevalence = class.probs[delta], sens.and.spec = sens.and.spec)
  }
}

fit_lcre_mcmc <- function(X, n.chains = 1, n.sample = 10000, n.thin = 1,
                        n.burnin = 800, n.adapt = 200, raw = FALSE,
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

  mod <- runjags::run.jags(mod.jags, n.chains = n.chains, sample = n.sample,
                           thin = n.thin, inits = inits, burnin = n.burnin,
                           adapt = n.adapt, method = runjags.method)

  if (isTRUE(raw)) {
    return(mod)
  } else {
    res <- summary(mod)[, "Mean"]
    class.probs <- res[grep("tau", names(res))]
    sens.fit <- res[grep("sens", names(res))]  # sensitivity
    spec.fit <- res[grep("spec", names(res))]  # specificity
    sens.and.spec <- data.frame(sens.fit, spec.fit)
    colnames(sens.and.spec) <- c("Sensitivity", "Specificity")
    rownames(sens.and.spec) <- getOption("diagacc.item.names")

    # Remove gold standard -----------------------------------------------------
    pos.of.gs <- getOption("diagacc.gold")
    if (!is.na(pos.of.gs)) sens.and.spec <- sens.and.spec[-pos.of.gs, ]

    list(prevalence = class.probs, sens.and.spec = sens.and.spec)
  }
}
