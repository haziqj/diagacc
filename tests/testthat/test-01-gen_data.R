context("Generate data")

test_that("LC data", {
  X <- gen_lc(n = 5, miss.prop = 1, seed = 123)
  expect_true(all(is.na(X[, ncol(X)])))

  XX <- gen_lc(seed = 123, drop.gs = TRUE)
  expect_equal(ncol(XX), ncol(X) - 1)
})

test_that("LCRE data", {
  X <- gen_lcre(n = 5, miss.prop = 1, seed = 123)
  expect_true(all(is.na(X[, ncol(X)])))

  XX <- gen_lcre(seed = 123, drop.gs = TRUE)
  expect_equal(ncol(XX), ncol(X) - 1)
})

test_that("FM data", {
  X <- gen_fm(n = 5, miss.prop = 1, seed = 123)
  expect_true(all(is.na(X[, ncol(X)])))

  XX <- gen_fm(seed = 123, drop.gs = TRUE)
  expect_equal(ncol(XX), ncol(X) - 1)
})

test_that("Item names in data gen", {
  X <- gen_lc(n = 5, seed = 123, name.items = FALSE)
  expect_equal(colnames(X), paste0("X", seq_len(ncol(X))))
})
