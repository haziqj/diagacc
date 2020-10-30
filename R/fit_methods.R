#' Print method for diagaccMod
#'
#' @param x diagaccMod object
#' @param ... Not used
#'
#' @return Print diagaccMod object values
#' @export
print.diagaccMod <- function(x, ...) {
  if (is.diagaccLC(x)) cat("Latent class model fit\n")
  if (is.diagaccLCRE(x)) cat("Latent class with random effects model fit\n")
  if (is.diagaccFM(x)) cat("Finite mixture model fit\n")

  if (!is.null(x$MCMC.res)) {
    tmp <- x$MCMC.res
  } else {
    tmp <- sim_res(list(x), part.of.sim = FALSE)
  }
  tmp <- paste0(utils::capture.output(iprior::dec_plac(tmp, 3)), collapse = "\n")
  cat(tmp)
}

#' Obtain fitted values
#'
#' @param x diagaccMod object
#' @param a Item 1
#' @param b Item 2
#' @param ... Not used
#'
#' @return Print fitted values
#' @export
fitted.diagaccMod <- function(x, a, b, ...) {
  prev <- x$prevalence
  sens.fit <- x$sens.and.spec[, "Sensitivity"]
  spec.fit <- x$sens.and.spec[, "Specificity"]
  p <- length(sens.fit)
  X <- x$X

  # class/delta = 1 is diseased, class/delta = 2 is non-diseased

  tmp <- rbind(sens.fit, 1 - spec.fit) # these are P[X=1|delta]

  tab <- expand.grid(X = 1:0, class = 1:2, item = seq_len(p))

  calc_prob <- function(X, class, item) {
    res <- tmp[class, item]
    if (X == 0) return(1 - res)
    if (X == 1) return(res)
  }

  tab$prob <- NA
  for (i in seq_len(nrow(tab))) {
    tab$prob[i] <- calc_prob(tab$X[i], tab$class[i], tab$item[i])
  }
  tab

  expd.list <- list(0:1)
  for (i in seq_len(p-1)) {
    expd.list <- c(expd.list, list(0:1))
  }
  names(expd.list) <- paste0("X", p:1)
  fit.tab <- expand.grid(expd.list)[, p:1]
  fit.tab$pattern <- apply(fit.tab, 1, function(x) paste0(x, collapse = ""))

  obs.counts <- table(apply(X, 1, function(x) paste0(x, collapse = "")))
  tmpp <- rep(0, nrow(fit.tab))
  tmpp[fit.tab$pattern %in% names(obs.counts)] <- obs.counts
  fit.tab$obs <- tmpp

  calc_prob2 <- function(x.pat, n = nrow(X)) {
    # p.class1
    tab.class1 <- as.data.frame(cbind(X = as.numeric(x.pat)[1:p], class = 1, item = 1:p))
    tab.class1$prob <- NA
    for (i in seq_len(nrow(tab.class1))) {
      tab.class1$prob[i] <- calc_prob(tab.class1$X[i], 1, tab.class1$item[i])
    }

    # p.class2
    tab.class2 <- as.data.frame(cbind(X = as.numeric(x.pat)[1:p], class = 2, item = 1:p))
    tab.class2$prob <- NA
    for (i in seq_len(nrow(tab.class2))) {
      tab.class2$prob[i] <- calc_prob(tab.class2$X[i], 2, tab.class2$item[i])
    }

    n * (
      prev * prod(tab.class1$prob) + (1 - prev) * prod(tab.class2$prob)
    )
  }

  fit.tab$exp <- NA
  fit.tab$exp <- apply(fit.tab[, 1:p], 1, calc_prob2)

  colnames(fit.tab) <- c(colnames(X), "Pattern", "Obs.", "Exp.")
  for (j in seq_len(p)) {
    fit.tab[, j] <- as.character(fit.tab[, j])
  }

  if (!missing(a) & !missing(b)) {
    ind1 <- which(fit.tab[, a] == 0 & fit.tab[, b] == 0)
    ind2 <- which(fit.tab[, a] == 0 & fit.tab[, b] == 1)
    ind3 <- which(fit.tab[, a] == 1 & fit.tab[, b] == 0)
    ind4 <- which(fit.tab[, a] == 1 & fit.tab[, b] == 1)
    tmp.tab <- fit.tab[c(ind1[1], ind2[1], ind3[1], ind4[1]), -(p + 1)]

    tmp.tab[1, (p + 1):(p + 2)] <- apply(fit.tab[ind1, -seq_len(p + 1)], 2, sum)
    tmp.tab[2, (p + 1):(p + 2)] <- apply(fit.tab[ind2, -seq_len(p + 1)], 2, sum)
    tmp.tab[3, (p + 1):(p + 2)] <- apply(fit.tab[ind3, -seq_len(p + 1)], 2, sum)
    tmp.tab[4, (p + 1):(p + 2)] <- apply(fit.tab[ind4, -seq_len(p + 1)], 2, sum)
    fit.tab <- tmp.tab[, c(a, b, p + 1, p + 2)]
  }

  # Adjusted table (observed cell counts > 0)
  # fit.tab <- fit.tab[fit.tab$Obs. > 0, ]
  rownames(fit.tab) <- NULL

  res.chisq <- sum((fit.tab$Obs. - fit.tab$Exp.) ^ 2 / fit.tab$Exp.)

  list(tab = fit.tab, chisq = res.chisq)
}
