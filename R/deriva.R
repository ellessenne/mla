#' @title Numerical approch to derivate.
#' @description Function to return the first, second derivate and the information score matrix. The central finite-difference and forward finite-difference will be used.
#' @param b Value of parameters to be optimized over.
#' @param funcpa Function to be minimized (or maximized), with argument the vector of parameters over which minimization is to take place. It should return a scalar result.
#' @return A list with the following elements:
#' * `v`, the information score matrix;
#' * `rl`, log-likelihood or likelihood of the model.
#' @references Donald W. Marquardt (1963). _An algorithm for least-squares estimation of nonlinear parameters_. Journal of the Society for Industrial and Applied Mathematics, 11(2):431--441
#' @author Daniel Commenges
#' @export
#' @examples
#' b <- 0.1
#' f <- function(b) {
#'   return((2 * b[1]**2 + 3 * b[1]))
#' }
#' d <- deriva(b = b, funcpa = f)
deriva <- function(b, funcpa) {
  m <- length(b)
  bh2 <- bh <- rep(0, m)
  v <- rep(0, (m * (m + 3) / 2))
  fcith <- fcith2 <- rep(0, m)
  # function
  rl <- funcpa(b)
  # gradient null
  for (i in 1:m) {
    bh <- bh2 <- b
    th <- max(1E-7, (1E-4 * abs(b[i])))
    bh[i] <- bh[i] + th
    bh2[i] <- bh2[i] - th

    fcith[i] <- funcpa(bh)
    fcith2[i] <- funcpa(bh2)
  }
  k <- 0
  m1 <- m * (m + 1) / 2
  l <- m1
  for (i in 1:m) {
    l <- l + 1
    bh <- b
    thn <- -max(1E-7, (1E-4 * abs(b[i])))
    v[l] <- -(fcith[i] - fcith2[i]) / (2 * thn)
    for (j in 1:i) {
      bh <- b
      k <- k + 1
      thi <- max(1E-7, (1E-4 * abs(b[i])))
      thj <- max(1E-7, (1E-4 * abs(b[j])))
      th <- thi * thj
      bh[i] <- bh[i] + thi
      bh[j] <- bh[j] + thj
      temp <- funcpa(bh)
      v[k] <- -(temp - (fcith[j]) - (fcith[i]) + rl) / th
    }
  }
  result <- list(v = v, rl = rl)
  return(result)
}
