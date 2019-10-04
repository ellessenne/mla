### First test
fit.tol <- 1e-8

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
  t1.1 <- mla::mla(par = c(8, 9), control = list(maxiter = 100, epsa = fit.tol, epsb = fit.tol, epsd = fit.tol), fn = f1)
  testthat::expect_equal(object = t1.1$par, expected = c(5, 6), tolerance = 1e-6)
})

testthat::test_that("Test #1.2 (with gradient)", {
  t1.2 <- mla::mla(par = c(8, 9), control = list(maxiter = 100, epsa = fit.tol, epsb = fit.tol, epsd = fit.tol), fn = f1, gr = gr)
  testthat::expect_equal(object = t1.2$par, expected = c(5, 6), tolerance = 1e-6)
})

testthat::test_that("Test #1.3 (with hessian)", {
  t1.3 <- mla::mla(par = c(8, 9), control = list(maxiter = 100, epsa = fit.tol, epsb = fit.tol, epsd = fit.tol), fn = f1, hessian = hes)
  testthat::expect_equal(object = t1.3$par, expected = c(5, 6), tolerance = 1e-6)
})

testthat::test_that("Test #1.4 (with gradient and hessian)", {
  t1.4 <- mla::mla(par = c(8, 9), control = list(maxiter = 100, epsa = fit.tol, epsb = fit.tol, epsd = fit.tol), fn = f1, gr = gr, hessian = hes)
  testthat::expect_equal(object = t1.4$par, expected = c(5, 6), tolerance = 1e-5)
})


### Second test
f2 <- function(b) {
  (b[1] + 10 * b[2])^2 + 5 * (b[3] - b[4])^2 + (b[2] - 2 * b[3])^4 + 10 * (b[1] - b[4])^4
}

testthat::test_that("Test #2.1", {
  t2.2 <- mla::mla(par = c(3, -1, 0, 1), control = list(maxiter = 100, epsa = fit.tol, epsb = fit.tol, epsd = fit.tol), fn = f2)
  testthat::expect_equal(object = length(t2.2$par), expected = 4)
  testthat::expect_equal(object = t2.2$par, expected = rep(0, 4), tolerance = 1e-3)
})

### Rosenbrock function
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

testthat::test_that("Test: Rosenbrock function", {
  t <- mla::mla(par = c(-1.2, 1), control = list(epsa = fit.tol, epsb = fit.tol, epsd = fit.tol), fn = fr)
  testthat::expect_equal(object = length(t$par), expected = 2)
  testthat::expect_equal(object = t$par, expected = c(1, 1), tolerance = 1e-4)

  t <- mla::mla(par = c(-1.2, 1), control = list(epsa = fit.tol, epsb = fit.tol, epsd = fit.tol), fn = fr, gr = grr)
  testthat::expect_equal(object = length(t$par), expected = 2)
  testthat::expect_equal(object = t$par, expected = c(1, 1), tolerance = 1e-6)
})

### Booth function
fun <- function(par) (par[1] + 2 * par[2] - 7)^2 + (2 * par[1] + par[2] - 5)^2

testthat::test_that("Test: Booth function", {
  t <- mla::mla(par = c(-1.2, 1), control = list(epsa = fit.tol, epsb = fit.tol, epsd = fit.tol), fn = fun)
  testthat::expect_equal(object = length(t$par), expected = 2)
  testthat::expect_equal(object = t$par, expected = c(1, 3), tolerance = 1e-6)
})
