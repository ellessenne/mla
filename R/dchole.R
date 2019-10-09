#' @keywords internal
dchole2 <- function(a, k, nq, idpos) {
  hessian <- matrix(0, nrow = k, ncol = k)
  hessian[upper.tri(hessian, diag = TRUE)] <- a[1:(k * (k + 1) / 2)]
  hessian[lower.tri(hessian, diag = FALSE)] <- t(hessian)[lower.tri(hessian, diag = FALSE)]
  asd <- chol(hessian)
  a[1:(k * (k + 1) / 2)] <- asd[upper.tri(asd, diag = TRUE)]

  # what are the gradients divided for?
  # for 2/3 parameters they are divided by the diag of the hessian...
  a[-(1:(k * (k + 1) / 2))] <- a[-(1:(k * (k + 1) / 2))] / diag(hessian)

  list(fu = a, k = k, nq = nq, idpos = idpos)
}
