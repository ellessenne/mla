\name{marqLevAlg-package}
\Rdversion{1.0}
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
     Version: \tab 1.0\cr
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

## first test

f1 <- function(b){	
	return(4*(b[1]-5)^2+(b[2]-6)^2)	
}

## derivate
gr <- function(b){
	return(c(8*(b[1]-5),2*(b[2]-6)))
}

## hessian
hes <- function(b){
	return(c(-8,0,-2))
}

### without parameters
test.marq <- marqLevAlg(m=2,maxiter=100,epsa=0.001,epsb=0.001,epsd=0.001,fn=f1)
test.marq

### without length of parameters
test.marq <- marqLevAlg(b=c(8,9),maxiter=100,epsa=0.001,epsb=0.001,epsd=0.001,fn=f1)
test.marq

### with parameters and length of parameters
test.marq <- marqLevAlg(b=c(8,9),m=2,maxiter=100,epsa=0.001,epsb=0.001,epsd=0.001,fn=f1)

### with gradiant
test.marq <- marqLevAlg(b=c(8,9),m=2,maxiter=100,epsa=0.001,epsb=0.001,epsd=0.001,fn=f1,gr=gr)

### with hessian
test.marq <- marqLevAlg(b=c(8,9),m=2,maxiter=100,epsa=0.001,epsb=0.001,epsd=0.001,fn=f1,gr=gr,hess=hes)
test.marq


## second test

f2 <- function(b){(b[1]+10*b[2])^2+5*(b[3]-b[4])^2+(b[2]-2*b[3])^4+10*(b[1]-b[4])^4}


### without parameters
test.marq2 <- marqLevAlg(m=4,maxiter=100,epsa=0.001,epsb=0.001,epsd=0.001,fn=f2)
test.marq2

### without length of parameters
test.marq2 <- marqLevAlg(b=c(3,-1,0,1),maxiter=100,epsa=0.001,epsb=0.001,epsd=0.001,fn=f2)
test.marq2

### with parameters and length of parameters
test.marq2 <- marqLevAlg(b=c(3,-1,0,1),m=4,maxiter=100,epsa=0.001,epsb=0.001,epsd=0.001,fn=f2)
test.marq2


## third test
f<-function(x){(x[1]^42*(x[2]-3)^22)} 
test.marq3 <- marqLevAlg(b=c(10,10),maxiter=10000,epsa=0.1^5,epsb=0.1^5,epsd=0.1^5,fn=f)
test.marq3


## test 4
f<-function(x){x^120} 
test.marq4 <- marqLevAlg(b=10,maxiter=10000,epsa=0.1^5,epsb=0.1^5,epsd=0.1^5,fn=f)
test.marq4

## tests were extracted from optimx help
## test 5

fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
    x1 <- x[1]
    x2 <- x[2]
    c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
       200 *      (x2 - x1 * x1))
}

test.marq5 <- marqLevAlg(b=c(-1.2,1),maxiter=10000,epsa=0.1^5,epsb=0.1^5,epsd=0.1^5,fn=fr)

## test 6

test.marq6 <- marqLevAlg(b=c(-1.2,1),maxiter=10000,epsa=0.1^5,epsb=0.1^5,epsd=0.1^5,fn=fr,gr=grr)


# genrose function code

genrose.f<- function(x, gs=NULL){ # objective function
## One generalization of the Rosenbrock banana valley function (n parameters)
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
        return(fval)
}

genrose.g <- function(x, gs=NULL){
# vectorized gradient for genrose.f
# Ravi Varadhan 2009-04-03
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	gg <- as.vector(rep(0, n))
	tn <- 2:n
	tn1 <- tn - 1
	z1 <- x[tn] - x[tn1]^2
	z2 <- 1 - x[tn]
	gg[tn] <- 2 * (gs * z1 - z2)
	gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1
	return(gg)
}

genrose.h <- function(x, gs=NULL) { ## compute Hessian
   if(is.null(gs)) { gs=100.0 }
	n <- length(x)
	hh<-matrix(rep(0, n*n),n,n)
	for (i in 2:n) {
		z1<-x[i]-x[i-1]*x[i-1]
		z2<-1.0-x[i]
                hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
                hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
                hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
                hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
	}
        return(hh)
}
## initial parameters
startx<-4*seq(1:10)/3.

test.marq7 <- marqLevAlg(b=startx,maxiter=10000,epsa=0.1^5,epsb=0.1^5,epsd=0.1^5,fn=genrose.f,gr=genrose.g, hess=genrose.h)



}

\keyword{marqLevAlg}
\keyword{algorithm}
\keyword{optimization}
\keyword{minimization}
\keyword{maximisation}
\keyword{package}