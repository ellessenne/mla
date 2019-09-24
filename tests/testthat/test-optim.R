# Comparison with 'optim'

## Test 1
f1 <- function(b) {
  return(4 * (b[1] - 5)^2 + (b[2] - 6)^2)
}
gr <- function(b) {
  return(c(8 * (b[1] - 5), 2 * (b[2] - 6)))
}
b <- c(8, 9)

testthat::test_that("Test 1, comparison with optim", {
  t1.mla <- mla::mla(par = b, fn = f1)
  t1.optim.nm <- optim(par = b, fn = f1)
  t1.optim.bfgs <- optim(par = b, fn = f1, method = "BFGS")
  testthat::expect_equal(object = t1.mla$par, expected = t1.optim.nm$par, tolerance = 1e-4)
  testthat::expect_equal(object = t1.mla$par, expected = t1.optim.bfgs$par, tolerance = 1e-6)
})

testthat::test_that("Test 1 with gradient, comparison with optim", {
  t1.gr.mla <- mla::mla(par = b, fn = f1, gr = gr)
  t1.gr.optim.nm <- optim(par = b, fn = f1, gr = gr)
  t1.gr.optim.bfgs <- optim(par = b, fn = f1, gr = gr, method = "BFGS")
  testthat::expect_equal(object = t1.gr.mla$par, expected = t1.gr.optim.nm$par, tolerance = 1e-4)
  testthat::expect_equal(object = t1.gr.mla$par, expected = t1.gr.optim.bfgs$par, tolerance = 1e-6)
})

## Test 2: Rosenbrock function
fr <- function(x) {
  x1 <- x[1]
  x2 <- x[2]
  100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
  x1 <- x[1]
  x2 <- x[2]
  c(
    -400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
    200 * (x2 - x1 * x1)
  )
}
b <- c(-1.2, 1)

testthat::test_that("Test 2: Rosenbrock function, comparison with optim", {
  t2.mla <- mla::mla(par = b, fn = fr)
  t2.optim.nm <- optim(par = b, fn = fr)
  t2.optim.bfgs <- optim(par = b, fn = fr, method = "BFGS")
  testthat::expect_equal(object = t2.mla$par, expected = t2.optim.nm$par, tolerance = 1e-3)
  testthat::expect_equal(object = t2.mla$par, expected = t2.optim.bfgs$par, tolerance = 1e-3)
})

testthat::test_that("Test 2: Rosenbrock function with gradient, comparison with optim", {
  t2.gr.mla <- mla::mla(par = b, fn = fr, gr = grr)
  t2.gr.optim.nm <- optim(par = b, fn = fr, gr = grr)
  t2.gr.optim.bfgs <- optim(par = b, fn = fr, gr = grr, method = "BFGS")
  testthat::expect_equal(object = t2.gr.mla$par, expected = t2.gr.optim.nm$par, tolerance = 1e-3)
  testthat::expect_equal(object = t2.gr.mla$par, expected = t2.gr.optim.bfgs$par, tolerance = 1e-3)
})

### Test 3: Booth function
fr <- function(par) (par[1] + 2 * par[2] - 7)^2 + (2 * par[1] + par[2] - 5)^2
b <- c(-1.2, 1)

testthat::test_that("Test: Booth function", {
  t3.mla <- mla::mla(par = b, fn = fr)
  t3.optim.nm <- optim(par = b, fn = fr)
  t3.optim.bfgs <- optim(par = b, fn = fr, method = "BFGS")
  testthat::expect_equal(object = t3.mla$par, expected = t3.optim.nm$par, tolerance = 1e-3)
  testthat::expect_equal(object = t3.mla$par, expected = t3.optim.bfgs$par, tolerance = 1e-3)
})
