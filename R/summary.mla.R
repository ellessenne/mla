#' @title Summary of optimization.
#' @description A short summary of parameters estimates via the `mla` algorithm.
#' @param object An object of class `mla`.
#' @param digits Number of digits to print in outputs. Default value is 8.
#' @param ... Unused.
#' @seealso mla
#' @seealso print.mla
#' @author Daniel Commenges
#' @author Melanie Prague
#' @author Amadou Diakite
#' @keywords summary
#' @export
#' @examples
#' f1 <- function(b) {
#'   return(4 * (b[1] - 5)^2 + (b[2] - 6)^2)
#' }
#' test.mla <- mla(
#'   b = c(8, 9), m = 2, maxiter = 100, epsa = 0.001, epsb = 0.001,
#'   epsd = 0.001, fn = f1
#' )
#' 
#' summary(test.mla)
summary.mla <- function(object, digits = 8, ...) {
  x <- object
  if (!inherits(x, "mla")) stop("use only with \"marqLevAlg\" objects")

  cat(" \n")
  cat("Values:", "\n")
  cat(" \n")
  id <- 1:length(x$b)
  indice <- rep(id * (id + 1) / 2)
  se <- sqrt(x$v[indice])
  wald <- (x$b / se)**2
  z <- abs(stats::qnorm((1 + .95) / 2))
  tmp <- data.frame("coef" = format(round(x$b, digits)), "SE coef" = format(round(se, digits)), "Wald" = format(wald, digits), "P-value" = format.pval(1 - stats::pchisq(wald, 1), digits = digits, eps = 0.0001))

  print(tmp, row.names = F)
}
