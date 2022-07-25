# Load libraries
library(pracma)

# Setting up least squares problem
# Advertising example from lecture notes
a <- 3:6
s <- c(105,117,141,152)
plot(a, s, pch = 19)
vandermonde <- function(x) {
  A <- outer(x, (length(x)):1 - 1, "^")
  return(A)
}
A <- vandermonde(a)[,c(3,4)]

# Least squares with one basis vector
# Find best solution x (scalar) to ax = b
a <- c(2,1)
b <- c(6,8)
pseudoinv <- solve(t(a) %*% a) %*% t(a)
x <- pseudoinv %*% b
bhat <- a %*% x
r <- b - bhat
print(x)
print(bhat)
dot <- function(v1,v2){sum(v1*v2)}
dot(a,r)

# Least squares with two basis vectors
# Advertising example from lecture notes
a <- 3:6
A <- vandermonde(a)[,c(3,4)]
b <- c(105,117,141,152)
pseudoinv <- solve(t(A) %*% A) %*% t(A)
x <- pseudoinv %*% b
bhat <- A %*% x
r <- b - bhat
print(x)
print(bhat)
dot(r,A[,1])
dot(r,A[,2])
plot(a1,b,xlab="advertising",ylab="sales")
xx <- seq(from=0,to=6,length=200)
lines(xx,horner(as.numeric(x),xx)$y)

# Best-fit parabola
x <- c(0,1,2,3)
A <- vandermonde(x)[,2:4]
b <- c(6,5,2,2)
pseudoinv <- solve(t(A) %*% A) %*% t(A)
c <- pseudoinv %*% b
bhat <- A %*% c
r <- b - bhat
print(c)
plot(x,b)
xx <- seq(from=0,to=3,length=200)
lines(xx,horner(c,xx)$y)
print(r)

# Data compression
x <- seq(from=0,to=1,length=100)
y <- 3*x + 0.4*(2*runif(100))
plot(x,y)
A <- vandermonde(x)[,(length(x)-1):length(x)]
c <- qr.solve(A,y)
lines(x, horner(c,x)$y, col = "red")