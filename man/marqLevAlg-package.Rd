\name{marqLevAlg-package}
\Rdversion{1.1}
\alias{marqLevAlg-package}
\docType{package}

\title{An algorithm for least-squares curve fitting.}

\description{
This algorithm provides a numerical solution to the problem of minimizing a function. This is more efficient than the Gauss-Newton-like algorithm when starting from points vey far from the final minimum. A new convergence test is implemented (RDM) in addition to the usual stopping criterion : stopping rule is when the gradients are small enough in the parameters metric (GH-1G). 
}

\details{
  \tabular{ll}{
     Package: \tab marqLevAlg\cr
     Type: \tab Package\cr
     Version: \tab 1.1\cr
     Date: \tab 2012-03-09\cr
     License: \tab GPL (>= 2.0)\cr
     LazyLoad: \tab yes\cr
  } 
This algorithm provides a numerical solution to the problem of minimizing a function. This is more efficient than the Gauss-Newton-like algorithm when starting from points vey far from the final minimum. A new convergence test is implemented (RDM) in addition to the usual stopping criterion : stopping rule is when the gradients are small enough in the parameters metric (GH-1G). 
}

\author{D. Commenges <Daniel.Commenges@isped.u-bordeaux2.fr>, M. Prague <Melanie.Prague@isped.u-bordeaux2.fr> and Amadou Diakite <Amadou.Diakite@isped.u-bordeaux2.fr>} 

\references{
 \emph{marqLevAlg Algorithm}
 
Donald W. marquardt An algorithm for Least-Squares Estimation of Nonlinear Parameters. Journal of the Society for Industrial and Applied Mathematics, Vol. 11, No. 2. (Jun, 1963), pp. 431-441.

\emph{Convergence criteria : Relative distance to Minimim}

Commenges D. Jacqmin-Gadda H. Proust C. Guedj J. A Newton-like algorithm for likelihood maximization the robust-variance scoring algorithm arxiv:math/0610402v2 (2006)

}

\examples{

### 1
### initial values
b <- c(8,9)
### your function
f1 <- function(b){	
	return(4*(b[1]-5)^2+(b[2]-6)^2)	
}
## Call
test1 <- marqLevAlg(b=b,fn=f1)

### 2
### initial values
b <- c(3,-1,0,1)
### your function
f2 <- function(b){	
	return((b[1]+10*b[2])^2+5*(b[3]-b[4])^2+(b[2]-2*b[3])^4+10*(b[1]-b[4])^4)	
}

## Call
test2 <- marqLevAlg(b=b,fn=f2)
test2
}

\keyword{marqLevAlg}
\keyword{algorithm}
\keyword{optimization}
\keyword{minimization}
\keyword{maximisation}
\keyword{package}