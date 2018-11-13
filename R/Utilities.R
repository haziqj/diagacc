read_diagacc_opt <- function(env) {
  sens <- getOption("diagacc.sens")
  spec <- getOption("diagacc.spec")
  item.names <- getOption("diagacc.item.names")
  p <- getOption("diagacc.pwg")

  list2env(list(sens = sens, spec = spec, item.names = item.names, p = p),
           env)
}
