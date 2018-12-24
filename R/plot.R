prep_plot_df <- function(x, type, monitor) {
  type <- match.arg(type, c("est", "bias", "mse", "sd", "sd.bias", "sd.mse"))
  monitor <- match.arg(monitor, c("all", "prev", "sens", "spec"))
  tmp <- sim_res2(x)
  typ.ind <- match(type, c("est", "bias", "mse", "sd", "sd.bias", "sd.mse"))
  plot.df <- data.frame(
    LC   = tmp$LC[, typ.ind],
    LCRE = tmp$LCRE[, typ.ind],
    FM   = tmp$FM[, typ.ind]
  )
  p <- (nrow(plot.df) - 1) / 2

  plot.df$param <- c("Prevalence", rep("Sensitivity", p), rep("Specificity", p))
  plot.df$name <- {
    tmp.nam <- rownames(plot.df)
    tmp.nam <- gsub("Sens.", "", tmp.nam)
    tmp.nam <- gsub("Spec.", "", tmp.nam)
    tmp.nam
  }
  plot.df$name <-  factor(plot.df$name, levels = getOption("diagacc.item.names"))
  rownames(plot.df) <- NULL

  plot.df$n <- extract_n(x)
  plot.df$tau <- extract_tau(x)
  plot.df$miss.prop <- paste0(extract_miss.prop(x) * 100, "%")
  plot.df$data.gen <- extract_data.gen(x)

  if (monitor == "all") {
    mon.ind <- seq_along(plot.df$param)
  } else {
    mon.ind <- grep(monitor, plot.df$param, ignore.case = TRUE)
  }

  plot.df[mon.ind, ]
}

#' @export
plot_paper <- function(x, type = "bias", monitor = "sens", n = 250, tau = 0.08,
                       data.gen = "lc") {
  # Plots for the paper. Takes in a diagaccSim2 object. Plot either bias, mse or
  # sd for 3 x data.gen and 3 x miss.prop
  sk.ind <- which(x$sim.key$n == n & x$sim.key$prevalence == tau)
  plot.df <- NULL
  for (i in sk.ind) {
    plot.df <- rbind(plot.df, prep_plot_df(x[[i]], type = type, monitor = monitor))
  }
  plot.df <- reshape2::melt(plot.df, var = "model", id = c("param", "name",
                                                           "data.gen", "n",
                                                           "tau", "miss.prop"))
  plot.df <- plot.df[plot.df$data.gen == data.gen, ]
  plot.df$miss.prop <- factor(plot.df$miss.prop, levels = c("20%", "50%", "100%"))

  pp <- ggplot(plot.df)
  if (monitor == "sens") {
    the.title <- "Sensitivities"
  }
  if (monitor == "spec") {
    the.title <- "Specificities"
  }
  the.title <- paste0(the.title, " under ", toupper(data.gen), " data gen.")
  if (type == "bias") {
    yaxis.lab <- "Bias"
    pp <- pp + geom_abline(intercept = 0, slope = 0, linetype = "dashed",
                           col = "grey")
  }
  if (type == "mse") {
    yaxis.lab <- "MSE"
  }
  if (type == "sd") {
    yaxis.lab <- "Posterior S.D."
  }
  if (type == "sd.bias") {
    yaxis.lab <- "Bias of Posterior S.D."
  }
  if (type == "sd.mse") {
    yaxis.lab <- "MSE of Posterior S.D."
  }

  pp +
    geom_point(aes(x = name, y = value)) +
    facet_grid(miss.prop ~ model) +
    labs(x = NULL, y = yaxis.lab) +
    ggtitle(the.title) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


#' Plot for diagaccSim1 object
#'
#' @param x A diagaccSim1 object.
#' @param type The plot type.
#' @param ... Not used.
#'
#' @export
plot.diagaccSim1 <- function(x, type = c("est", "bias", "sd"), ...) {

  plot.df <- sim_res(x, type)
  stuff <- rownames(plot.df[[1]])
  stuff <- gsub("Sens.", "", stuff)
  stuff <- gsub("Spec.", "", stuff)
  p <- (length(stuff) - 1) / 2

  plot.df <- as.data.frame(rbind(plot.df$LC, plot.df$LCRE, plot.df$FM))
  plot.df <- cbind(plot.df, name = stuff,
                   model = rep(c("LC", "LCRE", "FM"), each = length(stuff)))
  plot.df$model <- factor(plot.df$model, levels = c("LC", "LCRE", "FM"))
  plot.df$name <- factor(plot.df$name, levels = c(stuff[-seq_len(p + 1)],
                                                  "Prevalence"))
  plot.df$monitor <- c("", rep("Sensitivity", p), rep("Specificity", p))
  rownames(plot.df) <- NULL
  colnames(plot.df)[1] <- "y"

  the.title <- extract_sim.msg(x)

  pp <- ggplot(plot.df, aes(x = name, y = y, col = model))
  if (type == "bias") {
    pp <- pp + geom_hline(yintercept = 0, linetype = "dashed", col = "grey60") +
      labs(x = "", y = "Bias", col = "") +
      coord_cartesian(ylim = c(-0.5, 0.5))
  }
  if (type == "est") {
    pp <- pp + labs(x = "", y = "Estimate", col = "") +
      coord_cartesian(ylim = c(0, 1))
  }
  if (type == "sd") {
    pp <- pp + labs(x = "", y = "Posterior SD", col = "")
  }
  pp + geom_point(position = position_dodge(width = 0.25)) +
    geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.3,
                   position = position_dodge(width = 0.25)) +
    facet_grid(. ~ monitor, scales = "free", space = "free") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(the.title)

}

#' Plot diagaccSim2 study object
#'
#' @param x diagaccSim2 object
#' @param sim.key The study result to display
#' @param type Plot type, estimates or standard error.
#' @param ... Not used.
#'
#' @export
plot.diagaccSim2 <- function(x, sim.key, type = c("est", "se"), ...) {
  if (missing(sim.key)) {
    cat("Use plot(object, sim.key = <number>, ...), where <number> is from the following:\n")
    cat(paste0(utils::capture.output(x$sim.key), collapse = "\n"))
  } else {
    plot(x[[sim.key]], type = type)
  }
}
