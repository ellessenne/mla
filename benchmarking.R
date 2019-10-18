# benchmarking
library(marqLevAlg)
library(microbenchmark)
devtools::load_all()

ll <- function(par, data) {
  b0 <- par[1]
  b1 <- par[2]
  s <- par[3]
  n <- length(data$y)
  ll <- -n / 2 * log(2 * pi) - n * log(s) - 1 / (2 * s^2) * sum((data$y - (b0 + b1 * data$x))^2)
  return(-ll)
}

make_data <- function(n = 1e5) {
  x <- runif(n)
  y <- 1 + 0.5 * x + rnorm(n)
  return(data.frame(y, x))
}

mlaFun <- function() {
  mla::mla(par = rep(0.5, 3), fn = ll, data = make_data())
}

mlaFun()

marqLevAlgFun <- function() {
  data.in <- make_data()
  ll2 <- function(par) {
    b0 <- par[1]
    b1 <- par[2]
    s <- par[3]
    n <- length(data.in$y)
    ll <- -n / 2 * log(2 * pi) - n * log(s) - 1 / (2 * s^2) * sum((data.in$y - (b0 + b1 * data.in$x))^2)
    return(-ll)
  }
  marqLevAlg::marqLevAlg(b = rep(0.5, 3), fn = ll2)
}

bench <- microbenchmark::microbenchmark(
  "mla" = mlaFun(),
  "marqLevAlg" = marqLevAlgFun(),
  times = 100
)

library(ggplot2)
autoplot(bench)

