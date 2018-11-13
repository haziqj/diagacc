is.try_error <- function(x) inherits(x, "try-error")

run_sim <- function(B = 3, n = 250, tau = 0.08, miss.prop = 0.2,
                    data.gen = c("lc", "lcre", "fm")) {
  # Initialise -----------------------------------------------------------------
  data.gen <- match.arg(data.gen, c("lc", "lcre", "fm"))
  if (data.gen == "lc") gen_data <- gen_lc
  if (data.gen == "lcre") gen_data <- gen_lcre
  if (data.gen == "fm") gen_data <- gen_fm

  res.lc <- res.lcre <- res.fm <- list()

  if (!exists("pb")) {
    cat("Running simulation with these settings:\n")
    sim.msg <- paste0("n = ", n, ", tau = ", tau, ", miss prop = ", miss.prop,
                      ", data gen. mech. = ", toupper(data.gen), "\n")
    cat(sim.msg)

    pb <- progress::progress_bar$new(
      format = "  running simulations [:bar] :percent in :elapsed", total = B * 3,
      clear = FALSE, width = 90, show_after = 0
    )
    pb$tick(0)
  }

  i <- 1
  while (i <= B) {
    X <- gen_data(n = n, tau = tau, miss.prop = miss.prop)

    # Latent class model fit ---------------------------------------------------
    res.lc[[i]] <- try(fit_lc(X))
    pb$tick()

    # Latent class with random effects model fit -------------------------------
    res.lcre[[i]] <- try(fit_lcre(X, quad.points = 189), silent = TRUE)
    pb$tick()

    if (is.try_error(res.lc[[i]]) | is.try_error(res.lcre[[i]])) {
      pb$tick(-2)
    } else {
      # If no errors in randomLCA fit, proceed to FM fit (MCMC) ----------------
      suppressWarnings(res.fm[[i]] <- fit_fm(X, n.sample = 10000, silent = TRUE))
      i <- i + 1
      pb$tick()
    }
  }

  res <- list(LC = res.lc, LCRE = res.lcre, FM = res.fm, sim.msg = sim.msg)
  class(res) <- "diagaccSim1"
  res
}

sim_res <- function(res) {
  p <- getOption("diagacc.p")
  item.names <- getOption("diagacc.item.names")[1:p]

  prev.vec <- sapply(res, function(x) x[[1]])
  mean.prev <- mean(prev.vec)
  sd.prev <- sd(prev.vec)

  sens.and.spec.array <- array(unlist(lapply(res, function(x) x[[2]])), dim = c(5, 2, 2))
  mean.sens.and.spec <- apply(sens.and.spec.array, c(1, 2), mean)
  sd.sens.and.spec <- apply(sens.and.spec.array, c(1, 2), sd)

  tab <- data.frame(cbind(
    Mean = c(mean.prev, c(mean.sens.and.spec)),
    SE   = c(sd.prev, c(sd.sens.and.spec))
  ))
  tab <- cbind(tab, "2.5%" = tab$Mean - 1.96 * tab$SE,
               "97.5%" = tab$Mean + 1.96 * tab$SE)
  rownames(tab) <- c("Prevalence", paste0("Sens.", item.names),
                     paste0("Spec.", item.names))

  tab
}

print.diagaccSim1 <- function(x) {
  # Get the true values --------------------------------------------------------
  p <- getOption("diagacc.p")
  sens <- getOption("diagacc.sens")[1:p]
  spec <- getOption("diagacc.spec")[1:p]
  tau <- as.numeric(
    strsplit(strsplit(tmp$sim.msg, "tau = ")[[1]][2], ", miss prop")[[1]][1]
  )
  truth <- c(tau, sens, spec)

  # Output ---------------------------------------------------------------------
  cat("LC data generation\n")
  tmp <- cbind(sim_res(x$LC), truth)
  cat(paste0(capture.output(iprior::dec_plac(tmp, 3)), collapse = "\n"))

  cat("\n\nLCRE data generation\n")
  tmp <- cbind(sim_res(x$LCRE), truth)
  cat(paste0(capture.output(iprior::dec_plac(tmp, 3)), collapse = "\n"))

  cat("\n\nFM data generation\n")
  tmp <- cbind(sim_res(x$FM), truth)
  cat(paste0(capture.output(iprior::dec_plac(tmp, 3)), collapse = "\n"))
}
