read_diagacc_opt <- function(env) {
  sens <- getOption("diagacc.sens")
  spec <- getOption("diagacc.spec")
  item.names <- getOption("diagacc.item.names")
  p <- getOption("diagacc.pwg")

  list2env(list(sens = sens, spec = spec, item.names = item.names, p = p),
           env)
}

is.try_error <- function(x) inherits(x, "try-error")

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

which_class_diseased <- function(x) {
  # Helper function for the randomLCA method to determine which class is
  # diseased. Compares the probabilities P[X = 1 | row1] vs P[X = 1 | row2] and
  # determines that if the probabilities are larger in row1 then row2 is the
  # diseased class, and vice versa. The logic is that the larger probabilities
  # should imply the sensitivities of the tests, which should be high.
  #
  # Args: x is a (2 by no. of items) matrix. Each row is the probabilities P[X =
  # 1 | delta]. This is obtained from the fit of randomLCA.
  #
  # Returns: Numeric. The row number for the diseased class (1 or 2).
  tmp <- apply(x, 2, function(y) which(y == max(y)))
  as.numeric(names(sort(table(tmp)))[1])
}


