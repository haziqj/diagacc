context("LC model fit")

test_that("EM algorithm", {
  X <- gen_lc(n = 30)
  mod <- fit_lc(X, method = "EM")
  expect_is(mod, "diagaccMod")
  expect_true(is.diagaccLC(mod))

  mod <- fit_lc(X, method = "EM", raw = TRUE)
  expect_is(mod, "randomLCA")
})

test_that("MCMC", {
  X <- gen_lc(n = 30)
  suppressWarnings(
    mod <- fit_lc(X, method = "MCMC", n.adapt = 10, n.sample = 100, n.burnin = 10,
                  n.chains = 1, silent = TRUE)
  )
  expect_is(mod, "diagaccMod")
  expect_true(is.diagaccLC(mod))

  suppressWarnings(
    mod <- fit_lc(X, method = "MCMC", n.adapt = 10, n.sample = 100, n.burnin = 10,
                  n.chains = 2, raw = TRUE, silent = TRUE)
  )
  expect_is(mod, "runjags")
})

context("LCRE model fit")

test_that("EM algorithm", {
  X <- gen_lcre(n = 30, seed = 123)
  mod <- fit_lcre(X, method = "EM", quad.points = 21)
  expect_is(mod, "diagaccMod")
  expect_true(is.diagaccLCRE(mod))

  mod <- fit_lcre(X, method = "EM", raw = TRUE, quad.points = 21)
  expect_is(mod, "randomLCA")
})

test_that("MCMC", {
  X <- gen_lcre(n = 30)
  suppressWarnings(
    mod <- fit_lcre(X, method = "MCMC", n.adapt = 10, n.sample = 100, n.burnin = 10,
                    n.chains = 1, silent = TRUE)
  )
  expect_is(mod, "diagaccMod")
  expect_true(is.diagaccLCRE(mod))

  suppressWarnings(
    mod <- fit_lcre(X, method = "MCMC", n.adapt = 10, n.sample = 100, n.burnin = 10,
                    n.chains = 2, raw = TRUE, silent = TRUE)
  )
  expect_is(mod, "runjags")
})

context("FM model fit")

test_that("MCMC", {
  X <- gen_lcre(n = 30)
  suppressWarnings(
    mod <- fit_fm(X, n.adapt = 10, n.sample = 100, n.burnin = 10,
                  n.chains = 1, silent = TRUE)
  )
  expect_is(mod, "diagaccMod")
  expect_true(is.diagaccFM(mod))

  suppressWarnings(
    mod <- fit_fm(X, n.adapt = 10, n.sample = 100, n.burnin = 10,
                  n.chains = 2, raw = TRUE, silent = TRUE)
  )
  expect_is(mod, "runjags")
})

context("Item names")

test_that("Item names", {
  X <- gen_lc(n = 30)
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  mod <- fit_lc(X, method = "EM")
  expect_identical(rownames(mod$sens.and.spec), colnames(X))
})
