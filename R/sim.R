is.try_error <- function(x) inherits(x, "try-error")

#' Run a simulation
#'
#' Repeatedly fit the LC, LCRE and FM model onto a specific scenaria.
#'
#' The framework for the simulation is as follows: \enumerate{\item Randomly
#' generate a specific data set according to the \code{data.gen} mechanism and
#' the settings as per \code{n}, \code{tau}, and \code{miss.prop}. \item Fit a
#' LC, LCRE and FM model onto the data set. \item Repeat steps 1--2 a total of
#' \code{B} number of times.}
#'
#' Note that it is possible to continue on a simulation run by calling the saved
#' \code{object} in the argument. See the example section.
#'
#' @param object A diagaccSim1 object
#' @param B Number of replications
#' @param n (numeric) Sample size
#' @param tau (numeric, between 0 and 1) Prevalence
#' @param miss.prop (numeric, between 0 and 1) Missing proportion
#' @param data.gen The data generating mechanism
#' @param pb Not used
#' @param lc.method (DEPRECATED) The latent class model estimation method
#' @param lcre.method (DEPRECATED) The latent class model with random effects
#'   estimation method
#'
#' @return A diagaccSim1 object. It is a list of length four containing the
#'   following items: \itemize{\item \code{LC}: A list of length \code{B}
#'   containing the LC model fit using the \code{fit_lc()} function. \item
#'   \code{LCRE}: A list of length \code{B} containing the LCRE model fit using
#'   the \code{fit_lcre()} function. \item \code{FM}: A list of length \code{B}
#'   containing the LC model fit using the \code{fit_fm()} function. \item \code{sim.settings}: A list of length \code{B}
#'   containing the LC model fit using the \code{fit_lc()} function.}
#' @export
run_sim <- function(B = 4, n = 250, tau = 0.08, miss.prop = 0.2,
                    data.gen = c("lc", "lcre", "fm"), lc.method = "MCMC",
                    lcre.method = "MCMC", pb, object = NULL) {
  # Initialise -----------------------------------------------------------------
  if (!is.null(object)) {  # Add additional simulations
    n <- extract_n(object)
    tau <- extract_tau(object)
    miss.prop <- extract_miss.prop(object)
    data.gen <- extract_data.gen(object)
    lc.method <- extract_lc.method(object)
    lcre.method <- extract_lcre.method(object)
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
    res.lc[[i]] <- try(suppressWarnings(
      fit_lc(X, method = lc.method, silent = TRUE, gold.std = TRUE)
    ), silent = TRUE)
    pb$tick()

    # Latent class with random effects model fit -------------------------------
    res.lcre[[i]] <- try(suppressWarnings(
      fit_lcre(X, quad.points = 189, method = lcre.method, silent = TRUE,
               gold.std = TRUE)
    ), silent = TRUE)
    pb$tick()

    if (is.try_error(res.lc[[i]]) | is.try_error(res.lcre[[i]])) {
      pb$tick(-2)
    } else {
      # If no errors in randomLCA fit, proceed to FM fit (MCMC) ----------------
      suppressWarnings(res.fm[[i]] <- fit_fm(X, n.sample = 2000, silent = TRUE,
                                             gold.std = TRUE))
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
                       data.gen = data.gen, sim.msg = sim.msg,
                       lc.method = lc.method, lcre.method = lcre.method)

  res <- list(LC = res.lc, LCRE = res.lcre, FM = res.fm,
              sim.settings = sim.settings)
  class(res) <- "diagaccSim1"
  res
}
