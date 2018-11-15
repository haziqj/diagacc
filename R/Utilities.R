read_diagacc_opt <- function(env) {
  sens <- getOption("diagacc.sens")
  spec <- getOption("diagacc.spec")
  item.names <- getOption("diagacc.item.names")
  p <- getOption("diagacc.pwg")

  list2env(list(sens = sens, spec = spec, item.names = item.names, p = p),
           env)
}

is.diagaccSim1 <- function(x) inherits(x, "diagaccSim1")

is.diagaccSim2 <- function(x) inherits(x, "diagaccSim2")

extract_B <- function(x) {
  if (is.diagaccSim1(x)) return(x$sim.settings$B)
  if (is.diagaccSim2(x)) return(x$study.settings$B)
}

extract_n <- function(x) {
  if (is.diagaccSim1(x)) return(x$sim.settings$n)
  if (is.diagaccSim2(x)) return(x$study.settings$n)
}

extract_tau <- function(x) {
  if (is.diagaccSim1(x)) return(x$sim.settings$tau)
  if (is.diagaccSim2(x)) return(x$study.settings$tau)
}

extract_miss.prop <- function(x) {
  if (is.diagaccSim1(x)) return(x$sim.settings$miss.prop)
  if (is.diagaccSim2(x)) return(x$study.settings$miss.prop)
}

extract_data.gen <- function(x) {
  if (is.diagaccSim1(x)) return(x$sim.settings$data.gen)
  if (is.diagaccSim2(x)) return(x$study.settings$data.gen)
}

extract_sim.msg <- function(x) {
  if (is.diagaccSim1(x)) return(x$sim.settings$sim.msg)
  if (is.diagaccSim2(x)) return(x$study.settings$sim.msg)
}

true_vals <- function(x) {
  p <- getOption("diagacc.p")
  sens <- getOption("diagacc.sens")[1:p]
  spec <- getOption("diagacc.spec")[1:p]
  tau <- extract_tau(x)
  c(tau, sens, spec)
}
