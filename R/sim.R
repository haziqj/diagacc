is.try_error <- function(x) inherits(x, "try-error")

#' @export
run_sim <- function(object = NULL, B = 3, n = 250, tau = 0.08, miss.prop = 0.2,
                    data.gen = c("lc", "lcre", "fm"), pb) {
  # Initialise -----------------------------------------------------------------
  if (!is.null(object)) {  # Add additional simulations
    n <- extract_n(object)
    tau <- extract_tau(object)
    miss.prop <- extract_miss.prop(object)
    data.gen <- extract_data.gen(object)
  }
  data.gen <- match.arg(data.gen, c("lc", "lcre", "fm"))
  if (data.gen == "lc") gen_data <- gen_lc
  if (data.gen == "lcre") gen_data <- gen_lcre
  if (data.gen == "fm") gen_data <- gen_fm

  res.lc <- res.lcre <- res.fm <- list()

  sim.msg <- paste0("n = ", n, ", prev. = ", tau, ", missing gold = ", miss.prop * 100,
                    "%, data gen. mech. = ", toupper(data.gen), "\n")

  if (missing(pb)) {
    cat(paste0("Running ", B, ifelse(
      is.null(object),
      " replications\n",
      paste0(" additional replications (on top of ", extract_B(object), " runs)\n"
    ))))
    cat(sim.msg)

    pb <- progress::progress_bar$new(
      format = "  running simulations [:bar] :percent in :elapsed", total = 3 * B,
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
      suppressWarnings(res.fm[[i]] <- fit_fm(X, n.sample = 100, silent = TRUE))
      i <- i + 1
      pb$tick()
    }
  }

  if (!is.null(object)) {
    res.lc <- c(object$LC, res.lc)
    res.lcre <- c(object$LCRE, res.lcre)
    res.fm <- c(object$FM, res.fm)
    B <- extract_B(object) + B
  }

  sim.settings <- list(B = B, n = n, tau = tau, miss.prop = miss.prop,
                       data.gen = data.gen, sim.msg = sim.msg)

  res <- list(LC = res.lc, LCRE = res.lcre, FM = res.fm,
              sim.settings = sim.settings)
  class(res) <- "diagaccSim1"
  res
}

