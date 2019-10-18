#' @title An Algorithm for Least-squares Curve Fitting
#'
#' @description This algorithm provides a numerical solution to the
#' problem of minimizing a function. This is more efficient than
#' the Gauss-Newton-like algorithm when starting from points very
#' far from the final minimum. A new convergence test is
#' implemented (RDM) in addition to the usual stopping criterion:
#' stopping rule is when the gradients are small enough in the
#' parameters metric (GH-1G).
#'
#' @name mla
#' @docType package
#' @author Alessandro Gasparini (alessandro.gasparini@@ki.se)
#' @references Donald W. Marquardt (1963). _An algorithm for least-squares estimation of nonlinear parameters_. Journal of the Society for Industrial and Applied Mathematics, 11(2):431--441
#' @references Daniel Commenges,  H Jacqmin-Gadda, C. Proust, J. Guedj (2006). _A Newton-like algorithm for likelihood maximization the robust-variance scoring algorithm_. arxiv:math/0610402v2
NULL
