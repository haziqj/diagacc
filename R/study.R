#' @export
run_study <- function(object = NULL, B = 4, n = c(25, 100), tau = c(0.08, 0.4),
                      miss.prop = c(0.2, 0.5, 1.0), lc.method = c("EM", "MCMC"),
                      lcre.method = c("EM", "MCMC")) {
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
  # Initialise -----------------------------------------------------------------
  if (!is.null(object) & is.diagaccSim2(object)) {  # Add additional simulations
    n <- extract_n(object)
    tau <- extract_tau(object)
    miss.prop <- extract_miss.prop(object)
    lc.method <- extract_lc.method(object)
    lcre.method <- extract_lcre.method(object)
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
               extract_B(object), " + ", B, " replications each."))

  } else {
    cat(paste0("Running a total of ", nrow(sim.key), " simulation scenarios, ",
               B, " replications each."))
  }

  pb <- progress::progress_bar$new(
    format = "  running simulations [:bar] :percent in :elapsed",
    total = 3 * B * nrow(sim.key), clear = FALSE, width = 90, show_after = 0
  )

  res <- list()
  i <- 0
  for (dg in c("lc", "lcre", "fm")) {
    for (N in n) {
      for (TAU in tau) {
        for (MISSPROP in miss.prop) {
          i <- i + 1
          cat(paste0("\n[", i,"] n = ", N, " prev. = ", TAU, " missing gold = ",
                     MISSPROP, " data gen. mech. = ", toupper(dg), "\n"))
          res[[i]] <- run_sim(object[[i]], B, N, TAU, MISSPROP, dg, pb,
                              lc.method = lc.method, lcre.method = lcre.method)
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
