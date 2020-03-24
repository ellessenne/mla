# 'lm' estimation
ll <- function(par, data) {
  b0 <- par[1]
  b1 <- par[2]
  s <- par[3]
  n <- length(data$y)
  ll <- -n / 2 * log(2 * pi) - n * log(s) - 1 / (2 * s^2) * sum((data$y - (b0 + b1 * data$x))^2)
  return(-ll)
}
make_data <- function(n = 1e3) {
  x <- runif(n)
  y <- 1 + 0.5 * x + rnorm(n)
  return(data.frame(y, x))
}


testthat::test_that("Check average of estimates from B replication", {
  set.seed(3497657)
  B <- 100
  res <- vector(mode = "list", length = B)
  for (i in seq_along(res)) {
    data.tmp <- make_data()
    res[[i]] <- suppressWarnings(mla::mla(par = runif(3), fn = ll, data = data.tmp)$par)
  }
  res <- do.call(rbind, res)
  res <- apply(res, MARGIN = 2, FUN = mean)
  testthat::expect_equal(object = res, expected = c(1, 0.5, 1), tol = 1e-2)
})

testthat::test_that("Comparison with optim and nlminb", {
  data.tmp <- make_data(n = 1e5)
  b <- runif(3)
  t1.mla <- suppressWarnings(mla::mla(par = b, fn = ll, data = data.tmp))
  t1.optim.nm <- suppressWarnings(optim(par = b, fn = ll, data = data.tmp))
  t1.optim.bfgs <- suppressWarnings(optim(par = b, fn = ll, data = data.tmp, method = "BFGS"))
  t1.nlminb <- suppressWarnings(nlminb(start = b, objective = ll, data = data.tmp))

  testthat::expect_equal(object = t1.mla$par, expected = t1.optim.nm$par, tolerance = 1e-3)
  testthat::expect_equal(object = t1.mla$par, expected = t1.optim.bfgs$par, tolerance = 1e-3)
  testthat::expect_equal(object = t1.mla$par, expected = t1.nlminb$par, tolerance = 1e-3)
})

testthat::test_that("Comparison with lm", {
  data.tmp <- make_data(n = 1e5)
  b <- runif(3)
  fit.mla <- suppressWarnings(mla::mla(par = b, fn = ll, data = data.tmp))
  fit.lm <- suppressWarnings(lm(y ~ x, data = data.tmp))
  testthat::expect_equal(object = fit.mla$par, expected = c(coef(fit.lm), sqrt(var(fit.lm$residuals))), tolerance = 1e-3, check.attributes = FALSE)
})
