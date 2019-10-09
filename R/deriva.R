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
  g <- numDeriv::grad(func = funcpa, x = b)
  h <- -numDeriv::hessian(func = funcpa, x = b)
  result <- list(v = c(h[upper.tri(h, diag = TRUE)], g), rl = funcpa(b), h = h)
  return(result)
}
