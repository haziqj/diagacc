#' @export
fit_lc <- function(X, method = c("EM", "MCMC"), n.sample = 10000, n.chains = 1,
                   n.thin = 1,  n.burnin = 800, n.adapt = 200, raw = FALSE,
                   runjags.method = "rjags", silent = FALSE, calcSE = TRUE) {
  method <- match.arg(method, c("EM", "MCMC"))

  res <- NULL
  if (method == "EM") {
    res <- fit_lc_randomLCA(X = X, raw = raw, calcSE = calcSE)
  }
  if (method == "MCMC") {
    res <- fit_lc_mcmc(X = X, n.chains = n.chains, n.sample = n.sample,
                       n.thin = n.thin, n.burnin = n.burnin, n.adapt = n.adapt,
                       raw = raw, runjags.method = runjags.method,
                       silent = silent)
  }
  res
}


fit_lc_randomLCA <- function(X, raw = FALSE, calcSE = TRUE) {
  mod <- randomLCA::randomLCA(X, nclass = 2, probit = TRUE, calcSE = calcSE)

  if (isTRUE(raw)) {
    return(mod)
  } else {
    class.probs <- randomLCA::classProbs(mod)  # class probs P[delta]
    probs <- mod$outcomep  # P[X = 1 | delta]

    delta <- which_class_diseased(probs)

    sens.fit <- probs[delta, ]  # sensitivity
    spec.fit <- 1 - probs[-delta, ]  # specificity
    sens.and.spec <- data.frame(sens.fit, spec.fit)

    # Obtain standard errors ---------------------------------------------------
    se.prev <- mod$se[1]
    se.sens.and.spec <- as.data.frame(t(matrix(mod$se[-1], nrow = 2)))

    # Remove gold standard -----------------------------------------------------
    pos.of.gs <- getOption("diagacc.gold")
    if (!is.na(pos.of.gs)) {
      sens.and.spec <- sens.and.spec[-pos.of.gs, ]
      se.sens.and.spec <- se.sens.and.spec[-pos.of.gs, ]
    }

    colnames(sens.and.spec) <- colnames(se.sens.and.spec) <-
      c("Sensitivity", "Specificity")
    rownames(sens.and.spec) <- rownames(se.sens.and.spec) <-
      getOption("diagacc.item.names")[seq_len(getOption("diagacc.p"))]

    list(prevalence = class.probs[delta], sens.and.spec = sens.and.spec,
         se.prev = se.prev, se.sens.and.spec = se.sens.and.spec)
  }
}

fit_lc_mcmc <- function(X, n.sample = 100, n.chains = 1, n.thin = 1,
                        n.burnin = 800, n.adapt = 200, raw = FALSE,
                        runjags.method = "rjags", silent = FALSE) {
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

  mod <- runjags::run.jags(mod.jags, n.chains = n.chains, sample = n.sample,
                           thin = n.thin, inits = inits, burnin = n.burnin,
                           adapt = n.adapt, method = runjags.method)

  if (isTRUE(raw)) {
    return(mod)
  } else {
    res <- summary(mod)[, "Mean"]
    class.probs <- res[grep("tau", names(res))]
    sens.fit <- res[grep("sens", names(res))]  # sensitivity
    spec.fit <- res[grep("spec", names(res))]   # specificity
    sens.and.spec <- data.frame(sens.fit, spec.fit)

    # Obtain post. stand. deviation --------------------------------------------
    res <- summary(mod)[, "SD"]
    se.prev <- res[grep("tau", names(res))]
    se.sens <- res[grep("sens", names(res))]  # sensitivity
    se.spec <- res[grep("spec", names(res))]   # specificity
    se.sens.and.spec <- data.frame(se.sens, se.spec)

    colnames(sens.and.spec) <- colnames(se.sens.and.spec) <-
      c("Sensitivity", "Specificity")
    rownames(sens.and.spec) <- rownames(se.sens.and.spec) <-
      getOption("diagacc.item.names")

    # Remove gold standard -----------------------------------------------------
    pos.of.gs <- getOption("diagacc.gold")
    if (!is.na(pos.of.gs)) {
      sens.and.spec <- sens.and.spec[-pos.of.gs, ]
      se.sens.and.spec <- se.sens.and.spec[-pos.of.gs, ]
    }


    list(prevalence = class.probs, sens.and.spec = sens.and.spec,
         se.prev = se.prev, se.sens.and.spec = se.sens.and.spec)
  }
}
