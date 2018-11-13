diagacc_default_options <- list(
  diagacc.sens = c(0.60, 0.73, 0.90, 0.90, 0.95, 1.00),
  diagacc.spec = c(0.99, 0.45, 0.87, 0.50, 0.90, 1.00),
  diagacc.item.names = c("Microscopy", "Dipsticks", "CAA", "Antibody", "LAMP", "Gold std."),
  diagacc.gold = 6,
  diagacc.p = 5,
  diagacc.pwg = 6
)

restore_diagacc_options <- function() {
  options(diagacc_default_options)
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(diagacc_default_options) %in% names(op))
  if (any(toset)) options(diagacc_default_options[toset])

  invisible()
}

.onAttach <- function(libname, pkgname) {

  tab.msg1 <- "Using the following sensitivities and specificities for the tests.\nUse diagacc_opt() to change these settings.\n"

  tab <- cbind(
    getOption("diagacc.sens"),
    getOption("diagacc.spec")
  )
  rownames(tab) <- getOption("diagacc.item.names")
  colnames(tab) <- c("Sensitivity", "Specificity")
  tab.msg2 <- paste0(capture.output(tab), collapse = "\n")

  packageStartupMessage(paste0(tab.msg1, tab.msg2))
}


