#' Combine results from two separate simulation studies
#'
#' Combine results from two separate simulation studies
#'
#' This can be either a diagaccSim1 or diagaccSim2 object.
#'
#' @param res1 Simulation object 1
#' @param res2 Simulation object 2
#'
#' @export
combine_res <- function(res1, res2) {
  # Checks ---------------------------------------------------------------------
  combine.sim <- combine.study <- FALSE
  if (!is.diagaccSim1(res1) &
      !is.diagaccSim1(res2) &
      !is.diagaccSim2(res1) &
      !is.diagaccSim2(res2)) {
    stop("Input only diagaccSim objects")
  }
  if (all(is.diagaccSim1(res1), is.diagaccSim1(res2), TRUE)) {
    settings1 <- res1$sim.settings
    settings1$B <- NULL
    settings2 <- res2$sim.settings
    settings2$B <- NULL
    combine.sim <- TRUE
  } else if (all(is.diagaccSim2(res1), is.diagaccSim2(res2), TRUE)) {
    settings1 <- res1$study.settings
    settings1$B <- NULL
    settings2 <- res2$study.settings
    settings2$B <- NULL
    combine.study <- TRUE
  } else {
    stop("Cannot combine sim results with study results")
  }
  if (!all.equal(settings1, settings2)) {
    stop("The two objects have different simulation settings")
  }

  # Combine sim results --------------------------------------------------------
  if (isTRUE(combine.sim)) {
    res <- combine_sim(res1, res2)
  }
  if (isTRUE(combine.study)) {
    res <- mapply(FUN = combine_sim, res1[seq_len(nrow(res1$sim.key))],
                  res2[seq_len(nrow(res2$sim.key))], SIMPLIFY = FALSE)
    res$study.settings <- res1$study.settings
    res$study.settings$B <- extract_B(res1) + extract_B(res2)
    res$sim.key <- res1$sim.key
    class(res) <- "diagaccSim2"
  }
  res
}

combine_sim <- function(res1, res2) {
  res1$LC <- c(res1$LC, res2$LC)
  res1$LCRE <- c(res1$LCRE, res2$LCRE)
  res1$FM <- c(res1$FM, res2$FM)
  res1$sim.settings$B <- extract_B(res1) + extract_B(res2)
  res1
}
