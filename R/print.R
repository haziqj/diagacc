sim_tab_hat <- function(res, est.or.sd = c("est", "sd")) {
  # Helper function to tabulate parameter replicates of run_sim() per $LC, $LCRE
  # and $FM.
  #
  # Args: A diagaccSim1 object. Also choose between parameter estimates or their
  # post. SD.
  #
  # Returns: A (2p + 1) x B matrix, i.e. estimates along the rows, and
  # replicates along the columns.

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


#' Helper function to provide summary of simulation results
#'
#' This is an internal helper function to provide summary of simulation results.
#'
#' @param x A \code{diagaccSim1} object.
#' @param type Choose between \code{"est"}, \code{"bias"} or \code{"sd"}.
#' @param raw (Logical) Provide the raw matrix of replicates?
#'
#' @export
sim_res <- function(x, type = c("est", "bias", "sd"), raw = FALSE) {
  if (extract_miss.prop(x) == 1) {
    truth <- true_vals(x)
  } else {
    truth <- true_vals_with_gs(x)
  }
  type <- match.arg(type, c("est", "bias", "sd"))
  res <- list(
    LC = sim_tab_hat(x$LC, ifelse(type == "sd", "sd", "est")),
    LCRE = sim_tab_hat(x$LCRE, ifelse(type == "sd", "sd", "est")),
    FM = sim_tab_hat(x$FM, ifelse(type == "sd", "sd", "est"))
  )
  if (type == "bias") {
    res <- lapply(res, function(z) z - truth)
  }

  if (isTRUE(raw)) {
    return(res)
  } else {
    MEAN <- lapply(res, function(x) apply(x, 1, mean))
    SE <-  lapply(res, function(x) apply(x, 1, sd))
    res <- list(
      LC   = t(mapply(MEAN$LC, SE$LC, FUN = cbind)),
      LCRE = t(mapply(MEAN$LCRE, SE$LCRE, FUN = cbind)),
      FM   = t(mapply(MEAN$FM, SE$FM, FUN = cbind))
    )
    UL <- lapply(res, function(z) z[, 1] + qnorm(0.975) * z[, 2])
    LL <- lapply(res, function(z) z[, 1] - qnorm(0.975) * z[, 2])
    res <- mapply(res, LL, UL, FUN = cbind, SIMPLIFY = FALSE)

    colnames(res$LC) <- colnames(res$LCRE) <- colnames(res$FM) <- c(
      toupper(type), "SE", "2.5%", "97.5%"
    )
    res
  }
}

#' @export
print.diagaccSim1 <- function(x, type = c("est", "bias", "sd"), digits = 5, ...) {
  cat("Summary from", x$sim.settings$B, "replications (")
  sim.msg <- x$sim.settings$sim.msg
  cat(gsub("\\n", ")\n", sim.msg))
  cat("Statistic: ")
  type <- match.arg(type, c("est", "bias", "sd"))
  if (type == "est") type.print <- "Parameter estimates"
  if (type == "bias") type.print <- "Parameter bias"
  if (type == "sd") type.print <- "Posterior SD of parameters"
  cat(type.print, "\n\n")

  res <- sim_res(x, type)
  cat("LC model fit\n")
  print(round(res$LC, digits))
  cat("\nLCRE model fit\n")
  print(round(res$LCRE, digits))
  cat("\nFM model fit\n")
  print(round(res$FM, digits))
}

#' @export
print.diagaccSim2 <- function(x, sim.key, ...) {
  if (missing(sim.key)) {
    cat("Use print(object, sim.key = <number>), where <number> is from the following:\n")
    cat(paste0(utils::capture.output(x$sim.key), collapse = "\n"))
  } else {
    print(x[[sim.key]])
  }
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
