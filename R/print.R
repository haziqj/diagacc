sim_res <- function(res, part.of.sim = TRUE) {
  if (isTRUE(part.of.sim)) {
    p <- getOption("diagacc.p")
    item.names <- getOption("diagacc.item.names")[1:p]
  } else {
    item.names <- rownames(res[[1]]$sens.and.spec)
    p <- length(item.names)
  }

  prev.vec <- sapply(res, function(x) x$prevalence)
  mean.prev <- mean(prev.vec)
  sd.prev <- sd(prev.vec)

  se.prev.vec <- sapply(res, function(x) x$se.prev)
  se.mean.prev <- mean(se.prev.vec)
  se.sd.prev <- sd(se.prev.vec)

  sens.and.spec.array <- array(
    unlist(lapply(res, function(x) x$sens.and.spec[1:p, ])), dim = c(p, 2, 2)
  )
  mean.sens.and.spec <- apply(sens.and.spec.array, c(1, 2), mean)
  sd.sens.and.spec <- apply(sens.and.spec.array, c(1, 2), sd)

  se.sens.and.spec.array <- array(
    unlist(lapply(res, function(x) x$se.sens.and.spec)), dim = c(p, 2, 2)
  )
  se.mean.sens.and.spec <- apply(se.sens.and.spec.array, c(1, 2), mean)
  se.sd.sens.and.spec <- apply(se.sens.and.spec.array, c(1, 2), sd)

  tab.est <- data.frame(cbind(
    Mean = c(mean.prev, c(mean.sens.and.spec)),
    SE   = c(sd.prev, c(sd.sens.and.spec))
  ))
  tab.est <- cbind(Est = tab.est[, 1], "2.5%" = tab.est$Mean - 1.96 * tab.est$SE,
                   "97.5%" = tab.est$Mean + 1.96 * tab.est$SE)

  tab.se <- data.frame(cbind(
    Mean = c(se.mean.prev, c(se.mean.sens.and.spec)),
    SE   = c(se.sd.prev, c(se.sd.sens.and.spec))
  ))
  tab.se <- cbind(SE = tab.se[, 1], "2.5%" = tab.se$Mean - 1.96 * tab.se$SE,
                  "97.5%" = tab.se$Mean + 1.96 * tab.se$SE)

  if (isTRUE(part.of.sim)) {
    tab <- cbind(tab.est, tab.se)
  } else {
    tab <- data.frame(cbind(
      Mean = c(mean.prev, c(mean.sens.and.spec)),
      SE   = c(se.mean.prev, c(se.mean.sens.and.spec))
    ))
    tab$`2.5%` <- tab$Mean - 1.96 * tab$SE
    tab$`97.5%` <- tab$Mean + 1.96 * tab$SE
  }

  rownames(tab) <- c("Prevalence", paste0("Sens.", item.names),
                     paste0("Spec.", item.names))
  tab

}

sim_tab_hat <- function(res, est.or.sd = c("est", "sd")) {
  # The function to tabulate parameter replicates of run_sim() per $LC, $LCRE
  # and $FM

  est.or.sd <- match.arg(est.or.sd, c("est", "sd"))
  item.names <- rownames(res[[1]]$sens.and.spec)
  p <- length(item.names)

  if (est.or.sd == "est") {
    prev.hat <- sapply(res, function(x) x$prevalence)
    sens.and.spec <- lapply(res, function(x) c(x$sens.and.spec))
  }
  if (est.or.sd == "sd") {
    prev.hat <- sapply(res, function(x) x$se.prev)
    sens.and.spec <- lapply(res, function(x) c(x$se.sens.and.spec))
  }

  sens.hat <- matrix(
    unlist(lapply(sens.and.spec, function(y) y$Sensitivity)), nrow = p
  )
  spec.hat <- matrix(
    unlist(lapply(sens.and.spec, function(y) y$Specificity)), nrow = p
  )
  tab.hat <- rbind(prev.hat, sens.hat, spec.hat)
  rownames(tab.hat) <- c(
    "Prevalence", paste0("Sens.", item.names), paste0("Spec.", item.names)
  )

  tab.hat
}

sim_res2 <- function(x, raw = FALSE) {
  # Helper function to summarise results from diagaccSim1 object.
  if (extract_miss.prop(x) == 1) {
    truth <- true_vals(x)
  } else {
    truth <- true_vals_with_gs(x)
  }
  res.est <- list(
    LC = sim_tab_hat(x$LC, "est"),
    LCRE = sim_tab_hat(x$LCRE, "est"),
    FM = sim_tab_hat(x$FM, "est")
  )
  res.sd <- list(
    LC = sim_tab_hat(x$LC, "sd"),
    LCRE = sim_tab_hat(x$LCRE, "sd"),
    FM = sim_tab_hat(x$FM, "sd")
  )

  if (isTRUE(raw)) {
    return(list(EST = res.est, SD = res.sd))
  } else {
    # Get the mean, MSE, SD, MSE for SD?
    EST <- lapply(res.est, function(x) apply(x, 1, mean))
    BIAS <- lapply(res.est, function(x) apply(x - truth, 1, mean))
    MSE <- lapply(res.est, function(x) apply((x - truth) ^ 2, 1, mean))
    SD <- lapply(res.sd, function(x) apply(x, 1, mean))
    res <- list(
      LC   = t(mapply(EST$LC, BIAS$LC, MSE$LC, SD$LC, FUN = cbind)),
      LCRE = t(mapply(EST$LCRE, BIAS$LCRE, MSE$LCRE, SD$LCRE, FUN = cbind)),
      FM   = t(mapply(EST$FM, BIAS$FM, MSE$FM, SD$FM, FUN = cbind))
    )
    colnames(res$LC) <- colnames(res$LCRE) <- colnames(res$FM) <- c(
      "EST", "BIAS", "MSE", "SD"
    )
    return(res)
  }
}

#' @export
print.diagaccSim1 <- function(x, raw = FALSE, ...) {
  if (isTRUE(raw)) {
    sim_res2(x, raw = TRUE)
  } else {
    res <- sim_res2(x, raw = FALSE)
    cat("LC model fit\n")
    print(res$LC)
    cat("\nLCRE model fit\n")
    print(res$LCRE)
    cat("\nFM model fit\n")
    print(res$FM)
  }
}

# print.diagaccSim1 <- function(x, ...) {
#   cat("LC model fit\n")
#   tmp <- cbind(truth = true_vals(x), sim_res(x$LC))
#   tmp <- paste0(utils::capture.output(iprior::dec_plac(tmp, 3)), collapse = "\n")
#   tmp <- gsub('"', " ", tmp)
#   tmp <- gsub("truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
#               "  truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
#               tmp)
#   cat(gsub('"', " ", tmp))
#
#   cat("\n\nLCRE model fit\n")
#   tmp <- cbind(truth = true_vals(x), sim_res(x$LCRE))
#   tmp <- paste0(utils::capture.output(iprior::dec_plac(tmp, 3)), collapse = "\n")
#   tmp <- gsub('"', " ", tmp)
#   tmp <- gsub("truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
#               "  truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
#               tmp)
#   cat(gsub('"', " ", tmp))
#
#   cat("\n\nFM model fit\n")
#   tmp <- cbind(truth = true_vals(x), sim_res(x$FM))
#   tmp <- paste0(utils::capture.output(iprior::dec_plac(tmp, 3)), collapse = "\n")
#   tmp <- gsub('"', " ", tmp)
#   tmp <- gsub("truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
#               "  truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
#               tmp)
#   cat(gsub('"', " ", tmp))
#   cat("\n")
# }

#' @export
print.diagaccSim2 <- function(x, sim.key, ...) {
  if (missing(sim.key)) {
    cat("Use print(object, sim.key = <number>), where <number> is from the following:\n")
    cat(paste0(utils::capture.output(x$sim.key), collapse = "\n"))
  } else {
    print(x[[sim.key]])
  }
}
