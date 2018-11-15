#' @export
run_study <- function(object = NULL, B = 3, n = c(25, 100), tau = c(0.08, 0.4),
                      miss.prop = c(0.2, 0.5, 1.0)) {
  # res
  # |_ LC
  #     |_ N
  #     |_ TAU  (... 12 scenarios ...)
  #     |_ MISSPROP
  #           |_ 3 x results
  # |_ LCRE
  #     ...
  # |_ FM
  #     ...
  res <- list()
  no.scen <- length(n) * length(tau) * length(miss.prop) * 3  # data.gen

  pb <- progress::progress_bar$new(
    format = "  running simulations [:bar] :percent in :elapsed",
    total = 3 * B * no.scen, clear = FALSE, width = 90, show_after = 0
  )

  i <- 0
  for (dg in c("lc", "lcre", "fm")) {
    for (N in n) {
      for (TAU in tau) {
        for (MISSPROP in miss.prop) {
          i <- i + 1
          cat("\nN = ", N, " TAU = ", TAU, " MISSPROP = ", MISSPROP,
              " DATA GEN = ", dg, "\n")
          res[[i]] <- run_sim(object, B, N, TAU, MISSPROP, dg, pb)
          # names(res)[length(res)] <- extract_sim.msg(res[[i]])
        }
      }
    }
  }

  res
}

#' @export
run_study_par <- function(object = NULL, B = 4, n = c(25, 100),
                          tau = c(0.08, 0.4), miss.prop = c(0.2, 0.5, 1.0)) {
  # Initialise -----------------------------------------------------------------
  if (!is.null(object) & is.diagaccSim2(object)) {  # Add additional simulations
    n <- extract_n(object)
    tau <- extract_tau(object)
    miss.prop <- extract_miss.prop(object)
  }

  # Create sim key -------------------------------------------------------------
  sim.key <- expand.grid(
    "missing gold" = paste0(miss.prop * 100, "%"),
    "prevalence" = tau,
    "n" = n,
    "data gen. mech." = c("LC", "LCRE", "FM")
  )
  sim.key <- sim.key[, c(3, 2, 1, 4)]

  if (!is.null(object)) {
    cat(paste0("Running a total of ", nrow(sim.key), " simulation scenarios, ",
               extract_B(object), " + ", B, " replications each.\n"))

  } else {
    cat(paste0("Running a total of ", nrow(sim.key), " simulation scenarios, ",
               B, " replications each.\n"))
  }
  # print(n)

  res <- list()
  i <- 0
  for (dg in c("lc", "lcre", "fm")) {
    for (N in n) {
      for (TAU in tau) {
        for (MISSPROP in miss.prop) {
          i <- i + 1
          cat(paste0("[", i, "] "))
          res[[i]] <- run_sim_par(object[[i]], B, N, TAU, MISSPROP, dg)
          sim.msg <- extract_sim.msg(res[[i]])
          names(res)[length(res)] <- gsub("\n", "", sim.msg)
        }
      }
    }
  }

  if (!is.null(object)) B <- extract_B(object) + B

  study.settings <- list(B = B, n = n, tau = tau, miss.prop = miss.prop)
  res$study.settings <- study.settings
  res$sim.key <- sim.key
  class(res) <- "diagaccSim2"
  res
}

#' @export
print.diagaccSim2 <- function(x, sim.key) {
  if (missing(sim.key)) {
    cat("Use print(object, sim.key = <number>), where <number> is from the following:\n")
    cat(paste0(capture.output(x$sim.key), collapse = "\n"))
  } else {
    print(x[[sim.key]])
  }
}

#' @export
plot.diagaccSim2 <- function(x, sim.key, type = c("est", "se")) {
  if (missing(sim.key)) {
    cat("Use plot(object, sim.key = <number>, ...), where <number> is from the following:\n")
    cat(paste0(capture.output(x$sim.key), collapse = "\n"))
  } else {
    plot(x[[sim.key]], type = type)
  }
}

