context("LC model fit")

test_that("MCMC", {
  X <- gen_lc(n = 10)
  suppressWarnings(
    mod <- fit_lc(X, n.adapt = 1, n.sample = 100, n.burnin = 10, n.chains = 1,
                  silent = TRUE)
  )
  expect_is(mod, "diagaccMod")
  expect_true(is.diagaccLC(mod))

  suppressWarnings(
    mod <- fit_lc(X, n.adapt = 1, n.sample = 100, n.burnin = 0, n.chains = 2,
                  raw = TRUE, silent = TRUE)
  )
  expect_is(mod, "runjags")
})

context("LCRE model fit")

test_that("MCMC", {
  X <- gen_lcre(n = 10)
  suppressWarnings(
    mod <- fit_lcre(X, n.adapt = 1, n.sample = 100, n.burnin = 0, n.chains = 1,
                    silent = TRUE)
  )
  expect_is(mod, "diagaccMod")
  expect_true(is.diagaccLCRE(mod))

  suppressWarnings(
    mod <- fit_lcre(X, n.adapt = 1, n.sample = 100, n.burnin = 0, n.chains = 2,
                    raw = TRUE, silent = TRUE)
  )
  expect_is(mod, "runjags")
})

context("FM model fit")

test_that("MCMC", {
  X <- gen_lcre(n = 10)
  suppressWarnings(
    mod <- fit_fm(X, n.adapt = 1, n.sample = 100, n.burnin = 0, n.chains = 1,
                  silent = TRUE)
  )
  expect_is(mod, "diagaccMod")
  expect_true(is.diagaccFM(mod))

  suppressWarnings(
    mod <- fit_fm(X, n.adapt = 1, n.sample = 100, n.burnin = 0, n.chains = 2,
                  raw = TRUE, silent = TRUE)
  )
  expect_is(mod, "runjags")
})

context("Item names")

test_that("Item names", {
  X <- gen_lc(n = 10, name.items = FALSE)
  suppressWarnings(
    mod <- fit_fm(X, n.adapt = 1, n.sample = 100, n.burnin = 0, n.chains = 1,
                  silent = TRUE)
  )
  expect_identical(rownames(mod$sens.and.spec), colnames(X))
})
