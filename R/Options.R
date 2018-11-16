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

  tab <- get_sens_spec_tab()
  tab.msg2 <- paste0(capture.output(tab), collapse = "\n")

  packageStartupMessage(paste0(tab.msg1, tab.msg2))
}

get_sens_spec_tab <- function() {
  tab <- cbind(
    getOption("diagacc.sens"),
    getOption("diagacc.spec")
  )
  rownames(tab) <- getOption("diagacc.item.names")
  colnames(tab) <- c("Sensitivity", "Specificity")
  tab
}

#' @export
diagacc_opt <- function(sens, spec, item.names, pos.gold.std, default = FALSE) {
  if (isTRUE(default)) {
    restore_diagacc_options()
  } else {
    if (!missing(sens)) options(diagacc.sens = sens)
    if (!missing(spec)) options(diagacc.spec = spec)
    if (!missing(item.names)) options(diagacc.item.names = item.names)
    if (!missing(pos.gold.std)) options(diagacc.gold = pos.gold.std)
    check_diagacc_opt()
  }
  print(get_sens_spec_tab())
}

check_diagacc_opt <- function() {
  sens <- getOption("diagacc.sens")
  spec <- getOption("diagacc.spec")
  if (length(sens) != length(spec)) {
    restore_diagacc_options()
    stop("Length of sensitivities and specificities do not match.")
  } else {
    pos.gold.std <- getOption("diagacc.gold")
    p <- length(sens) - length(pos.gold.std)
    options(diagacc.p = p)
    options(diagacc.pwg = p + length(pos.gold.std))
  }
  if (any(sens > 1) | any(sens < 0)) {
    restore_diagacc_options()
    stop("Sensitivities must be between 0 and 1.")
  }
  if (any(spec > 1) | any(spec < 0)) {
    restore_diagacc_options()
    stop("Specificities must be between 0 and 1.")
  }
}
