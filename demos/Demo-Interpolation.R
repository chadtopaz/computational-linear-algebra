# Load libraries
library(pracma)
library(knitr)
library(matlib)

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

# Vandermonde (regularized) vs. Lagrange speed test for 100 pts
set.seed(123)
numTrials <- 1000
n <- 10
x <- 1:n
y <- runif(n)
x0 <- 1.5
tic()
for (i in 1:numTrials){
  c <- echelon(vander(x), y)[, n+1]
  horner(C,x0)
}
t1 <- toc(echo = FALSE)
tic()
for (i in 1:numTrials){
  lagrangeInterp(x,y,x0)
}
t2 <- toc(echo = FALSE)
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