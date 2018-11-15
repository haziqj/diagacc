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
print.diagaccSim1 <- function(x) {
  cat("LC model fit\n")
  tmp <- cbind(truth = true_vals(x), sim_res(x$LC))
  tmp <- paste0(capture.output(iprior::dec_plac(tmp, 3)), collapse = "\n")
  tmp <- gsub('"', " ", tmp)
  tmp <- gsub("truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
              "  truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
              tmp)
  cat(gsub('"', " ", tmp))

  cat("\n\nLCRE model fit\n")
  tmp <- cbind(truth = true_vals(x), sim_res(x$LCRE))
  tmp <- paste0(capture.output(iprior::dec_plac(tmp, 3)), collapse = "\n")
  tmp <- gsub('"', " ", tmp)
  tmp <- gsub("truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
              "  truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
              tmp)
  cat(gsub('"', " ", tmp))

  cat("\n\nFM model fit\n")
  tmp <- cbind(truth = true_vals(x), sim_res(x$FM))
  tmp <- paste0(capture.output(iprior::dec_plac(tmp, 3)), collapse = "\n")
  tmp <- gsub('"', " ", tmp)
  tmp <- gsub("truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
              "  truth    Est      2.5%     97.5%    SE       2.5%     97.5%",
              tmp)
  cat(gsub('"', " ", tmp))
}

#' @export
plot.diagaccSim1 <- function(x, type = c("est", "se")) {
  type <- match.arg(type, c("est", "se"))
  p <- getOption("diagacc.p")
  item.names <- getOption("diagacc.item.names")

  plot.df <- rbind(sim_res(x$LC), sim_res(x$LCRE), sim_res(x$FM))
  if (type == "est") {
    plot.df <- as.data.frame(plot.df[, 1:3])
    plot.df <- cbind(plot.df,
                     x = c("Prevalence", rep(item.names[1:p], 2)),
                     model = rep(c("LC", "LCRE", "FM"), each = 2 * p + 1))
    plot.df$Est <- plot.df$Est - true_vals(x)
    plot.df$`2.5%` <- plot.df$`2.5%` - true_vals(x)
    plot.df$`97.5%` <- plot.df$`97.5%` - true_vals(x)
    plot.df$model <- factor(plot.df$model, levels = c("LC", "LCRE", "FM"))
    plot.df$x <- factor(plot.df$x, levels = c(item.names, "Prevalence"))
    rownames(plot.df) <- NULL
    plot.df$Type <- c("", rep("Sensitivity", 5), rep("Specificity", 5))

    the.title <- extract_sim.msg(x)

    p <- ggplot(plot.df, aes(x = x, y = Est, col = model)) +
      geom_hline(yintercept = 0, linetype = "dashed", col = "grey60") +
      geom_point(position = position_dodge(width = 0.25)) +
      geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.3,
                     position = position_dodge(width = 0.25)) +
      facet_grid(. ~ Type, scales = "free", space = "free") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      labs(x = "", y = "Bias", col = "") +
      ggtitle(the.title) +
      coord_cartesian(ylim = c(-0.5, 0.5))
  }
  if (type == "se") {
    plot.df <- as.data.frame(plot.df[, -(1:3)])
    plot.df <- cbind(plot.df,
                     x = c("Prevalence", rep(item.names[1:p], 2)),
                     model = rep(c("LC", "LCRE", "FM"), each = 2 * p + 1))
    plot.df$model <- factor(plot.df$model, levels = c("LC", "LCRE", "FM"))
    plot.df$x <- factor(plot.df$x, levels = c(item.names, "Prevalence"))
    rownames(plot.df) <- NULL
    plot.df$Type <- c("", rep("Sensitivity", 5), rep("Specificity", 5))

    the.title <- extract_sim.msg(x)

    p <- ggplot(plot.df, aes(x = x, y = SE, col = model)) +
      geom_hline(yintercept = 0, linetype = "dashed", col = "grey60") +
      geom_point(position = position_dodge(width = 0.25)) +
      geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.3,
                     position = position_dodge(width = 0.25)) +
      facet_grid(. ~ Type, scales = "free", space = "free") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      labs(x = "", y = "SE/Post. SD of Estimate", col = "") +
      ggtitle(the.title)
  }

  p
}
