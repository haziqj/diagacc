sim_res <- function(res) {
  p <- getOption("diagacc.p")
  item.names <- getOption("diagacc.item.names")[1:p]

  prev.vec <- sapply(res, function(x) x$prevalence)
  mean.prev <- mean(prev.vec)
  sd.prev <- sd(prev.vec)

  se.prev.vec <- sapply(res, function(x) x$se.prev)
  se.mean.prev <- mean(se.prev.vec)
  se.sd.prev <- sd(se.prev.vec)

  sens.and.spec.array <- array(
    unlist(lapply(res, function(x) x$sens.and.spec)), dim = c(p, 2, 2)
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

  tab <- cbind(tab.est, tab.se)
  rownames(tab) <- c("Prevalence", paste0("Sens.", item.names),
                     paste0("Spec.", item.names))

  tab
}



#' @export
print.diagaccSim1 <- function(x, ...) {
  cat("LC model fit\n")
  tmp <- cbind(truth = true_vals(x), sim_res(x$LC))
  tmp <- paste0(utils::capture.output(iprior::dec_plac(tmp, 3)), collapse = "\n")
  tmp <- gsub('"', " ", tmp)
  tmp <- gsub("truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
              "  truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
              tmp)
  cat(gsub('"', " ", tmp))

  cat("\n\nLCRE model fit\n")
  tmp <- cbind(truth = true_vals(x), sim_res(x$LCRE))
  tmp <- paste0(utils::capture.output(iprior::dec_plac(tmp, 3)), collapse = "\n")
  tmp <- gsub('"', " ", tmp)
  tmp <- gsub("truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
              "  truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
              tmp)
  cat(gsub('"', " ", tmp))

  cat("\n\nFM model fit\n")
  tmp <- cbind(truth = true_vals(x), sim_res(x$FM))
  tmp <- paste0(utils::capture.output(iprior::dec_plac(tmp, 3)), collapse = "\n")
  tmp <- gsub('"', " ", tmp)
  tmp <- gsub("truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
              "  truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
              tmp)
  cat(gsub('"', " ", tmp))
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
