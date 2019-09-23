### First test
fit.tol <- 1e-6

f1 <- function(b) {
  return(4 * (b[1] - 5)^2 + (b[2] - 6)^2)
}
gr <- function(b) {
  return(c(8 * (b[1] - 5), 2 * (b[2] - 6)))
}
hes <- function(b) {
  return(c(-8, 0, -2))
}

testthat::test_that("Test #1.1", {
  t1.1 <- mla::mla(m = 2, maxiter = 100, epsa = fit.tol, epsb = fit.tol, epsd = fit.tol, fn = f1)
  testthat::expect_equal(object = t1.1$b, expected = c(5, 6), tolerance = 1e-6)
})

testthat::test_that("Test #1.2", {
  t1.2 <- mla::mla(b = c(8, 9), maxiter = 100, epsa = fit.tol, epsb = fit.tol, epsd = fit.tol, fn = f1)
  testthat::expect_equal(object = t1.2$b, expected = c(5, 6), tolerance = 1e-6)
})

testthat::test_that("Test #1.3", {
  t1.3 <- mla::mla(b = c(8, 9), m = 2, maxiter = 100, epsa = fit.tol, epsb = fit.tol, epsd = fit.tol, fn = f1)
  testthat::expect_equal(object = t1.3$b, expected = c(5, 6), tolerance = 1e-6)
})

testthat::test_that("Test #1.4 (with gradient)", {
  t1.4 <- mla::mla(b = c(8, 9), m = 2, maxiter = 100, epsa = fit.tol, epsb = fit.tol, epsd = fit.tol, fn = f1, gr = gr)
  testthat::expect_equal(object = t1.4$b, expected = c(5, 6), tolerance = 1e-6)
})

testthat::test_that("Test #1.5 (with hessian)", {
  t1.5 <- mla::mla(b = c(8, 9), m = 2, maxiter = 100, epsa = fit.tol, epsb = fit.tol, epsd = fit.tol, fn = f1, hess = hes)
  testthat::expect_equal(object = t1.5$b, expected = c(5, 6), tolerance = 1e-6)
})

testthat::test_that("Test #1.6 (with gradient and hessian)", {
  t1.6 <- mla::mla(b = c(8, 9), m = 2, maxiter = 100, epsa = fit.tol, epsb = fit.tol, epsd = fit.tol, fn = f1, gr = gr, hess = hes)
  testthat::expect_equal(object = t1.6$b, expected = c(5, 6), tolerance = 1e-5)
})


### Second test
f2 <- function(b) {
  (b[1] + 10 * b[2])^2 + 5 * (b[3] - b[4])^2 + (b[2] - 2 * b[3])^4 + 10 * (b[1] - b[4])^4
}

testthat::test_that("Test #2.1", {
  t2.1 <- mla::mla(m = 4, maxiter = 100, epsa = fit.tol, epsb = fit.tol, epsd = fit.tol, fn = f2)
  testthat::expect_equal(object = length(t2.1$b), expected = 4)
  testthat::expect_equal(object = t2.1$b, expected = rep(0, 4), tolerance = 1e-3)
})

testthat::test_that("Test #2.2", {
  t2.2 <- mla::mla(b = c(3, -1, 0, 1), maxiter = 100, epsa = fit.tol, epsb = fit.tol, epsd = fit.tol, fn = f2)
  testthat::expect_equal(object = length(t2.2$b), expected = 4)
  testthat::expect_equal(object = t2.2$b, expected = rep(0, 4), tolerance = 1e-3)
})

testthat::test_that("Test #2.2", {
  t2.3 <- mla::mla(b = c(3, -1, 0, 1), m = 4, maxiter = 100, epsa = fit.tol, epsb = fit.tol, epsd = fit.tol, fn = f2)
  testthat::expect_equal(object = length(t2.3$b), expected = 4)
  testthat::expect_equal(object = t2.3$b, expected = rep(0, 4), tolerance = 1e-3)
})


### Rosenbrock Banana function
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

testthat::test_that("Test: Rosenbrock Banana function", {
  t <- mla::mla(m = 2, epsa = fit.tol, epsb = fit.tol, epsd = fit.tol, fn = fr)
  testthat::expect_equal(object = length(t$b), expected = 2)
  testthat::expect_equal(object = t$b, expected = c(1, 1), tolerance = 1e-4)

  t <- mla::mla(b = c(-1.2, 1), epsa = fit.tol, epsb = fit.tol, epsd = fit.tol, fn = fr)
  testthat::expect_equal(object = length(t$b), expected = 2)
  testthat::expect_equal(object = t$b, expected = c(1, 1), tolerance = 1e-4)

  t <- mla::mla(m = 2, epsa = fit.tol, epsb = fit.tol, epsd = fit.tol, fn = fr, gr = grr)
  testthat::expect_equal(object = length(t$b), expected = 2)
  testthat::expect_equal(object = t$b, expected = c(1, 1), tolerance = 1e-6)

  t <- mla::mla(b = c(-1.2, 1), epsa = fit.tol, epsb = fit.tol, epsd = fit.tol, fn = fr, gr = grr)
  testthat::expect_equal(object = length(t$b), expected = 2)
  testthat::expect_equal(object = t$b, expected = c(1, 1), tolerance = 1e-6)
})
