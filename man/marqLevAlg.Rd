\name{marqLevAlg}
\alias{marqLevAlg}
\title{
An algorithm for least-squares curve fitting.
}
\description{
This algorithm provides a numerical solution to the problem of minimizing a function. This is more efficient than the Gauss-Newton-like algorithm when starting from points vey far from the final minimum. A new convergence test is implemented (RDM) in addition to the usual stopping criterion : stopping rule is when the gradients are small enough in the parameters metric (GH-1G). 
}
\usage{
marqLevAlg(b, m = FALSE, fn, gr=NULL, hess = NULL, maxiter = 500,
epsa = 0.001, epsb = 0.001, epsd = 0.01, digits = 8,
print.info = FALSE, blinding = TRUE, multipleTry = 25)
}

\arguments{
  \item{b}{
an optional vector containing the initial values for the parameters. Default is 0.1 for every parameter.
}
  \item{m}{
an optional parameter if the vector of parameter is not missing compulsory if b is not given.  
}
  \item{fn}{
The function to be minimized (or maximized), with argument the vector of parameters over which minimization is
to take place.  It should return a scalar result.
}
  \item{gr}{a function to return the gradient value for a specific point. If missing, finite-difference approximation will be used.
}
  \item{hess}{a function to return the hessian matrix for a specific point. If missing, finite-difference approximation will be used.
}
  \item{maxiter}{
optional maximum number of iterations for the marqLevAlg iterative algorithm. Default is 500.
}
  \item{epsa}{
optional threshold for the convergence criterion based on the parameter stability. Default is 0.001.
}
  \item{epsb}{
optional threshold for the convergence criterion based on the log-likelihood stability. Default is 0.001.
}
  \item{epsd}{
optional threshold for the relative distance to minimum. This criterion has  the nice interpretation of estimating the ratio of the approximation error over the statistical error, thus it can be used for stopping the iterative process whathever the problem. Default is 0.01.
}

  \item{digits}{Number of digits to print in outputs. Default value is 8.
}
  \item{print.info}{Logical.Equals to TRUE if report (parameters at iteration, function value, convergence criterion ...) at each iteration is requested. Default value is FALSE.
}
  \item{blinding}{Logical. Equals to TRUE if the algorithm is allowed to go on in case of an infinite or not definite value of function. Default value is FALSE.}
  \item{multipleTry}{Integer, different from 1 if the algorithm is allowed to go for the first iteration in case of an infinite or not definite value of gradients or hessian. This account for a starting point to far from the definition set. As many tries as requested in multipleTry will be done by changing the starting point of the algorithm. Default value is 25.}
}

\details{
Convergence criteria are very strict as they are based on derivatives of the log-likelihood in addition to the parameter and log-likelihood stability.
In some cases, the program may not converge and reach the maximum number of iterations fixed at 500. 
In this case, the user should check that parameter estimates at the last iteration are not on the boundaries of the parameter space. 
If the parameters are on the boundaries of the parameter space, the identifiability of the model should be assessed. 
If not, the program should be run again with other initial values, with a higher maximum number of iterations or less strict convergence tolerances.
}
\value{
  \item{cl}{
  summary of the call to the function marqLevAlg.
  }
   \item{ni}{
    number of marqLevAlg iterations before reaching stopping criterion.
  }
  \item{istop}{
 status of convergence: =1 if the convergence criteria were satisfied, =2 if the maximum number of iterations was reached, =4 if the algorithm encountered a problem in the function computation. 
  }
  \item{v}{
vector containing the upper triangle matrix of variance-covariance estimates at the stopping point.  
  }
  \item{fn.value}{
  function evaluation at the stopping point.
  }
  \item{b}{
stopping point value.
  }
  \item{ca}{
    convergence criteria for parameters stabilisation.
  }
  \item{cb}{
    convergence criteria for function stabilisation.  
  } 
  \item{rdm}{
    convergence criteria on the relative distance to minimum.  
  }   
  \item{time}{
    a running time.  
  }    
}
\references{
 \emph{marqLevAlg Algorithm}
 
Donald W. marquardt An algorithm for Least-Squares Estimation of Nonlinear Parameters. Journal of the Society for Industrial and Applied Mathematics, Vol. 11, No. 2. (Jun, 1963), pp. 431-441.

\emph{Convergence criteria : Relative distance to Minimim}

Commenges D. Jacqmin-Gadda H. Proust C. Guedj J. A Newton-like algorithm for likelihood maximization the robust-variance scoring algorithm arxiv:math/0610402v2 (2006)

}
\author{
D. Commenges - M. Prague - A. Diakite
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