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

run_study_par <- function(object = NULL, B = 4, n = c(25, 100), tau = c(0.08, 0.4),
                      miss.prop = c(0.2, 0.5, 1.0)) {

  # Create sim key -------------------------------------------------------------
  sim.key <- expand.grid(
    "missing gold" = miss.prop,
    "prevalence" = tau,
    "n" = n,
    "data gen. mech." = c("LC", "LCRE", "FM")
  )
  sim.key <- sim.key[, c(3, 2, 1, 4)]
  cat(paste0("Running a total of ", nrow(sim.key), " simulation scenarios.\n"))

  res <- list()

  i <- 0
  for (dg in c("lc", "lcre", "fm")) {
    for (N in n) {
      for (TAU in tau) {
        for (MISSPROP in miss.prop) {
          i <- i + 1
          cat(paste0("[", i, "] "))
          res[[i]] <- run_sim_par(object, B, N, TAU, MISSPROP, dg)
          sim.msg <- extract_sim.msg(res[[i]])
          names(res)[length(res)] <- gsub("\n", "", sim.msg)
        }
      }
    }
  }

  class(res) <- "diagaccSim2"
  res
}



# 1) Function to combine
# 2) Print
# 3) Plot
