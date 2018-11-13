#' @export
fit_fm <- function(X, n.chains = 1, n.sample = 1000, n.thin = 1,
                   n.burnin = 800, n.adapt = 200, raw = FALSE,
                   runjags.method = "rjags", silent = FALSE) {
  Y <- as.matrix(X)
  I <- nrow(Y)
  J <- ncol(Y)
  K <- 2

  inits <- list(
    l    = rep(0, I),
    pi.d = 0.1,
    eta  = c(0.1, 0.1),
    w    = matrix(c(0.9, 0.7), nrow = J, ncol = 2, byrow = TRUE)
  )

  # class 1 = diseased, class 2 = healthy
  mod.jags <- "model{
    for (i in 1:I) {
      for (j in 1:J) {
        Y[i,j] ~ dbern(A[i,j])
        A[i,j] <- pow(w[j, 2] * d[i] + (1 - w[j,1]) * (1 - d[i]), 1 - l[i]) - l[i] * (1 - d[i])
      }
      l[i] ~ dbern(pi.l[i])
      pi.l[i] <- (1 - d[i]) * eta[1] + d[i] * eta[2]
      d[i] ~ dbern(pi.d)
    }

    # Priors and sens/spec
    pi.d ~ dbeta(1,1)
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
  #monitor# pi.d, sens, spec, eta
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
    class.probs <- res[grep("pi.d", names(res))]
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
