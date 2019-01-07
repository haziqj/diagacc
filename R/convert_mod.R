convert_mod_diagacc <- function(mod, X) {
  item.names <- colnames(X)

  if (inherits(mod, "randomLCA")) {
    class.probs <- randomLCA::classProbs(mod)  # class probs P[delta]
    probs <- mod$outcomep  # P[X = 1 | delta]

    delta <- which_class_diseased(probs)
    class.probs <- class.probs[delta]  # get prevalence

    sens.fit <- probs[delta, ]  # sensitivity
    spec.fit <- 1 - probs[-delta, ]  # specificity
    sens.and.spec <- data.frame(sens.fit, spec.fit)

    # Obtain standard errors ---------------------------------------------------
    se.prev <- mod$se[1]
    se.sens.and.spec <- as.data.frame(
      t(matrix(mod$se[-1][seq_len(2 * ncol(X))], nrow = 2))
    )
  }

  if (inherits(mod, "runjags")) {
    res <- summary(mod)[, c("Mean", "SD", "Lower95", "Upper95")]
    class.probs <- res[grep("tau", rownames(res)), "Mean"]
    sens.fit <- res[grep("sens", rownames(res)), "Mean"]  # sensitivity
    spec.fit <- res[grep("spec", rownames(res)), "Mean"]  # specificity
    sens.and.spec <- data.frame(sens.fit, spec.fit)

    # Obtain post. stand. deviation --------------------------------------------
    se.prev <- res[grep("tau", rownames(res)), "SD"]
    se.sens <- res[grep("sens", rownames(res)), "SD"]  # sensitivity
    se.spec <- res[grep("spec", rownames(res)), "SD"]   # specificity
    se.sens.and.spec <- data.frame(se.sens, se.spec)

    # Prepare table of results for MCMC run ------------------------------------
    res <- as.data.frame(res)
    res <- res[c(grep("tau", rownames(res)),
                 grep("sens", rownames(res)),
                 grep("spec", rownames(res))), ]
    colnames(res) <- c("Mean", "SE", "2.5%", "97.5%")
    rownames(res) <- c("Prevalence", paste0("Sens.", item.names),
                       paste0("Spec.", item.names))

    # Deviance -----------------------------------------------------------------
    the.dev <- mod$deviance.sum  # first entry is deviance, second entry is penalty
    names(the.dev) <- NULL
  } else {
    res <- NULL
  }

  # # Remove gold standard -----------------------------------------------------
  # pos.of.gs <- getOption("diagacc.gold")
  # if (!is.na(pos.of.gs)) {
  #   sens.and.spec <- sens.and.spec[-pos.of.gs, ]
  #   se.sens.and.spec <- se.sens.and.spec[-pos.of.gs, ]
  # }

  colnames(sens.and.spec) <- colnames(se.sens.and.spec) <-
    c("Sensitivity", "Specificity")
  rownames(sens.and.spec) <- rownames(se.sens.and.spec) <-
    item.names

  res <- list(prevalence = class.probs, sens.and.spec = sens.and.spec,
              se.prev = se.prev, se.sens.and.spec = se.sens.and.spec,
              X = X, diagaccmod = mod$diagaccmod, MCMC.res = res,
              deviance = the.dev[1], DIC = the.dev[1] + the.dev[2])
  class(res) <- "diagaccMod"
  return(res)
}
