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
    res <- summary(mod)[, "Mean"]
    class.probs <- res[grep("tau", names(res))]
    sens.fit <- res[grep("sens", names(res))]  # sensitivity
    spec.fit <- res[grep("spec", names(res))]  # specificity
    sens.and.spec <- data.frame(sens.fit, spec.fit)

    # Obtain post. stand. deviation --------------------------------------------
    res <- summary(mod)[, "SD"]
    se.prev <- res[grep("tau", names(res))]
    se.sens <- res[grep("sens", names(res))]  # sensitivity
    se.spec <- res[grep("spec", names(res))]   # specificity
    se.sens.and.spec <- data.frame(se.sens, se.spec)
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
              X = X, diagaccmod = mod$diagaccmod)
  class(res) <- "diagaccMod"
  return(res)
}
