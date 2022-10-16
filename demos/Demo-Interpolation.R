# Load libraries
library(pracma)
library(knitr)
library(tictoc)

# What is interpolation?
# Suppose "true" function is cosine
x <- seq(from=0,to=2*pi,length=100)
y <- cos(x)
xdata <- seq(from=0,to=2*pi,length=4)
ydata <- cos(xdata)
yinterp <- lagrangeInterp(xdata, ydata, x)
plot(x,y,type="l",col="black")
lines(x,yinterp,col="red")

# Function vs. Taylor vs. Interpolating
f <- function(x){1/x}
t <- function(x){3-3*x+x^2}
p <- function(x){11/6 - x + 1/6*x^2}
xdata <- seq(from=1,to=3,length=3)
x <- seq(from=1,to=3,length=200)
plot(x,f(x),type="l",lwd=2)
lines(x,t(x),col="red",lwd=2)
points(xdata,f(xdata),cex=2)
py <- lagrangeInterp(xdata,f(xdata),x)
lines(x,py,col="green",lwd=2)

# Conditioning of the Vandermonde matrix
nvals <- 2^(1:8)
kappavals <- NULL
for (n in nvals){
  x <- seq(from=1,to=3,length=n)
  kappavals <- c(kappavals,kappa(vander(x)))
}
kable(cbind(nvals,kappavals), format = "simple")

# Solution with Vandermonde
x <- c(1,2,3,4)
y <- c(10,6,4,10)
C <- solve(vander(x),y)

# Vandermonde vs. Lagrange for 20 pts
set.seed(123)
n <- 20
x <- 1:n
y <- runif(n)
x0 <- 1.5
C <- solve(vander(x),y)
lagrangeInterp(x,y,x0)

# Vandermonde vs. Lagrange speed test for 10 pts
set.seed(123)
numTrials <- 100000
n <- 10
x <- 1:n
y <- runif(n)
tic()
for (i in 1:numTrials){
  C <- solve(vander(x),y)
  horner(C,x0)
}
T1 <- toc()
t1 <- T1$toc - T1$tic
tic()
for (i in 1:numTrials){
  lagrangeInterp(x,y,x0)
}
T2 <- toc()
t2 <- T2$toc - T2$tic
as.numeric(t2/t1)

# Set up function to do interpolation and calculate error
interperror <- function(f,xsample,x,n,plotflag=FALSE){
  ysample <- f(xsample)
  y <- f(x)
  yinterp <- lagrangeInterp(xsample,ysample,x)
  if (plotflag==TRUE){
    plot(x,y,type="l")
    points(xsample,ysample)
    lines(x,yinterp,col="blue")
  }
  error  <- max(abs(y-yinterp))
  return(error)
}

# Test for sin(x)
f <- function(x) {sin(x)}
n <- 5
xsample <- seq(from=0,to=2*pi,length=n)
x <- seq(from=0,to=2*pi,length=1000)
interperror(f,xsample,x,n,plotflag=TRUE)

# Interpolation error for sin(x) as a function of n
f <- function(x) {sin(x)}
x <- seq(from=0,to=2*pi,length=1000)
nvec <- 3:20
errorvec <- NULL
for (n in nvec){
  xsample <- seq(from=0,to=2*pi,length=n)
  errorvec <- c(errorvec,interperror(f,xsample,x,n,FALSE))
}
orderofmag <- round(log10(errorvec))
plot(nvec,orderofmag)

# Interpolation error for 1/(1+x^2) as a function of n
f <- function(x) {1/(1+x^2)}
x <- seq(from=-5,to=5,length=1000)
nvec <- seq(from = 5, to = 80, by = 5)
errorvec <- NULL
for (n in nvec){
  xsample <- seq(from=-5,to=5,length=n)
  errorvec <- c(errorvec,interperror(f,xsample,x,n,FALSE))
}
orderofmag <- round(log10(errorvec))
plot(nvec,orderofmag)

# Comparison of polynomial portion of error for equal spacing and Chebyshev
nvec <- 3:30
x <- seq(from=-1,to=1,length=5000)
equallyspaced <- NULL
for (n in nvec){
  nodes <- seq(from=-1,to=1,length=n)
  prod <- 1
  for (i in 1:n){
    prod <- prod*(x-nodes[i])
  }
  equallyspaced <- c(equallyspaced,unique(max(abs(prod))))
}
chebyshev <- 1/2^(nvec-1)
plot(nvec,log10(chebyshev),col="green",pch=16,ylim=c(-10,0),xlab="points",ylab="bound on portion of error")
points(nvec,log10(equallyspaced),col="red",pch=16)

# Maximum of derivative of function on interval
n <- 0:10
maxderiv <- c(0.398942,0.241971,0.178032, 0.550588,1.19683,2.30711,4.24061,14.178,41.8889,115.091,302.425) # Computed in Mathematica
plot(n,log(maxderiv))

# Interpolation with equal spacing
a <- -10
b <- 10
xexact <- seq(from=a,to=b,length=10000)
f <- function(x){exp(-x^2/2)/sqrt(2*pi)}
yexact <- f(xexact)
plot(xexact,yexact,type="l",lwd=3,xlim=c(a,b),ylim=c(-0.2,0.5))
n <- 30
xequal <- seq(from=a,to=b,length=n)
yequal <- interp(xequal,f(xequal),xexact)
lines(xexact,yequal,col="red",lwd=2)

# Comparison: interpolation with linear lookup table for given tolerance
n <- 1
error <- Inf
while (error > 0.5e-4){
  n <- n+1
  xlookup <- seq(from=a,to=b,length=n)
  ylookup <- f(xlookup)
  ytable <- approx(xlookup,ylookup,xexact)$y
  error <- max(abs(yexact-ytable))
}
nlookup <- n
print(nlookup)

# Comparison: interpolation with Chebyshev for given tolerance
n <- 1
error <- Inf
while (error > 0.5e-4){
  n <- n+1
  odds <- seq(from=1,to=2*n-1,by=2)
  xcheb <- (b+a)/2 + (b-a)/2*cos(odds*pi/2/n)
  ycheb <- interp(xcheb,f(xcheb),xexact)
  error <- max(abs(yexact-ycheb))
}
ncheb <- n
print(ncheb)

# Compression factor
nlookup/ncheb
