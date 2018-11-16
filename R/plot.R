
#' @export
plot.diagaccSim1 <- function(x, type = c("est", "se"), ...) {
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



#' @export
plot.diagaccSim2 <- function(x, sim.key, type = c("est", "se"), ...) {
  if (missing(sim.key)) {
    cat("Use plot(object, sim.key = <number>, ...), where <number> is from the following:\n")
    cat(paste0(utils::capture.output(x$sim.key), collapse = "\n"))
  } else {
    plot(x[[sim.key]], type = type)
  }
}

