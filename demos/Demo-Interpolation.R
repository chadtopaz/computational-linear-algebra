# Load libraries
library(pracma)
library(knitr)
library(tictoc)

# What is interpolation?
# Suppose "true" function is cosine
x <- seq(from=0,to=2*pi,length=100)
y <- cos(xgrid)
xdata <- seq(from=0,to=2*pi,length=4)
ydata <- cos(xdata)
yinterp <- lagrangeInterp(xdata, ydata, x)
plot(x,ytrue,type="l",col="black")
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

# Interpolation error for sin(x) with 5 pts
interperror <- function(n,plotflag=FALSE){
  xsample <- seq(from=0,to=2*pi,length=n)
  ysample <- sin(xsample)
  x <- seq(from=0,to=2*pi,length=1000)
  y <- sin(x)
  yinterp <- lagrangeInterp(xsample,ysample,x)
  if (plotflag==TRUE){
    plot(x,y,type="l")
    points(xsample,ysample)
    lines(x,yinterp,col="blue")
  }
  error  <- max(abs(y-yinterp))
  return(error)
}
interperror(5,plotflag=TRUE)

# Interpolation error for sin(x) as a function of n
nvec <- 3:20
errorvec <- NULL
for (n in nvec){
  errorvec <- c(errorvec,interperror(n))
}
orderofmag <- round(log10(errorvec))
plot(nvec,orderofmag)