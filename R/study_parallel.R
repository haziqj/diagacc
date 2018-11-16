#' @export
run_study_par <- function(object = NULL, B = 4, n = c(25, 100),
                          tau = c(0.08, 0.4), miss.prop = c(0.2, 0.5, 1.0),
                          no.cores = parallel::detectCores()) {
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

  res <- list()
  i <- 0
  for (dg in c("lc", "lcre", "fm")) {
    for (N in n) {
      for (TAU in tau) {
        for (MISSPROP in miss.prop) {
          i <- i + 1
          cat(paste0("[", i, "] "))
          res[[i]] <- run_sim_par(object[[i]], B, N, TAU, MISSPROP, dg, no.cores)
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