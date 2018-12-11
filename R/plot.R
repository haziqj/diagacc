prep_plot_df <- function(x, type, monitor) {
  type <- match.arg(type, c("est", "bias", "mse", "sd"))
  monitor <- match.arg(monitor, c("all", "prev", "sens", "spec"))
  tmp <- sim_res2(x)
  typ.ind <- grep(type, colnames(tmp[[1]]), ignore.case = TRUE)
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

plot_paper <- function(x, type = c("bias"), monitor = c("sens"), n = 250,
                       tau = 0.08, data.gen = "lc") {
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

  ggplot(plot.df) +
    geom_abline(intercept = 0, slope = 0, linetype = "dashed", col = "grey") +
    geom_point(aes(x = name, y = value)) +
    facet_grid(miss.prop ~ model) +
    labs(x = NULL, y = tools::toTitleCase(type)) +
    ggtitle(paste0("Sensitivities under ", toupper(data.gen),
                   " data generating mechanism")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

}


#' Plot for diagaccSim1 object
#'
#' @param x diagaccSim1 object.
#' @param type Plot type, estimates or standard error.
#' @param ... Not used.
#'
#' @export
plot.diagaccSim1 <- function(x, type = c("est", "bias", "mse", "sd"),
                             monitor = c("all", "prev", "sens", "spec"), ...) {

# tmp.nam





#
#   p <- getOption("diagacc.p")
#   item.names <- getOption("diagacc.item.names")
#
#   plot.df <- rbind(sim_res(x$LC), sim_res(x$LCRE), sim_res(x$FM))
#   if (type == "est") {
#     plot.df <- as.data.frame(plot.df[, 1:3])
#     plot.df <- cbind(plot.df,
#                      x = c("Prevalence", rep(item.names[1:p], 2)),
#                      model = rep(c("LC", "LCRE", "FM"), each = 2 * p + 1))
#     plot.df$Est <- plot.df$Est - true_vals(x)
#     plot.df$`2.5%` <- plot.df$`2.5%` - true_vals(x)
#     plot.df$`97.5%` <- plot.df$`97.5%` - true_vals(x)
#     plot.df$model <- factor(plot.df$model, levels = c("LC", "LCRE", "FM"))
#     plot.df$x <- factor(plot.df$x, levels = c(item.names, "Prevalence"))
#     rownames(plot.df) <- NULL
#     plot.df$Type <- c("", rep("Sensitivity", 5), rep("Specificity", 5))
#
#     the.title <- extract_sim.msg(x)
#
#     p <- ggplot(plot.df, aes(x = x, y = Est, col = model)) +
#       geom_hline(yintercept = 0, linetype = "dashed", col = "grey60") +
#       geom_point(position = position_dodge(width = 0.25)) +
#       geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.3,
#                      position = position_dodge(width = 0.25)) +
#       facet_grid(. ~ Type, scales = "free", space = "free") +
#       theme_bw() +
#       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#       labs(x = "", y = "Bias", col = "") +
#       ggtitle(the.title) +
#       coord_cartesian(ylim = c(-0.5, 0.5))
#   }
#   if (type == "se") {
#     plot.df <- as.data.frame(plot.df[, -(1:3)])
#     plot.df <- cbind(plot.df,
#                      x = c("Prevalence", rep(item.names[1:p], 2)),
#                      model = rep(c("LC", "LCRE", "FM"), each = 2 * p + 1))
#     plot.df$model <- factor(plot.df$model, levels = c("LC", "LCRE", "FM"))
#     plot.df$x <- factor(plot.df$x, levels = c(item.names, "Prevalence"))
#     rownames(plot.df) <- NULL
#     plot.df$Type <- c("", rep("Sensitivity", 5), rep("Specificity", 5))
#
#     the.title <- extract_sim.msg(x)
#
#     p <- ggplot(plot.df, aes(x = x, y = SE, col = model)) +
#       geom_hline(yintercept = 0, linetype = "dashed", col = "grey60") +
#       geom_point(position = position_dodge(width = 0.25)) +
#       geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.3,
#                      position = position_dodge(width = 0.25)) +
#       facet_grid(. ~ Type, scales = "free", space = "free") +
#       theme_bw() +
#       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#       labs(x = "", y = "SE/Post. SD of Estimate", col = "") +
#       ggtitle(the.title)
#   }
#
#   p
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

