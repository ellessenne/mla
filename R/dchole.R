#' @keywords internal
dchole2 <- function(a, k, nq, idpos) {
  hessian <- matrix(0, nrow = k, ncol = k)
  hessian[upper.tri(hessian, diag = TRUE)] <- a[1:(k * (k + 1) / 2)]
  hessian[lower.tri(hessian, diag = FALSE)] <- t(hessian)[lower.tri(hessian, diag = FALSE)]
  asd <- chol(hessian)
  a[1:(k * (k + 1) / 2)] <- asd[upper.tri(asd, diag = TRUE)]

  # the last elements are the solution of the hessian and gradients:
  # it's the step pk of the BFGS algorithm!
  a[-(1:(k * (k + 1) / 2))] <- solve(hessian, a[-(1:(k * (k + 1) / 2))])

  list(fu = a, k = k, nq = nq, idpos = idpos)
}
