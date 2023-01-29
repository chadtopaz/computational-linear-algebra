# Load libraries
library(tictoc)

# Set up x values
xvec <- runif(10^6)

# Polynomial evaluation - naive
tic()
for (x in xvec){
  2*x^4 + 3*x^3 - 3*x^2 + 5*x - 1
}
T1 <- toc()
t1 <- T1$toc - T1$tic

# Polynomial evaluation - Horner
tic()
for (x in xvec){
  -1 + x*(5 + x*(-3 + x*(3 + 2*x)))
}
T2 <- toc()
t2 <- T2$toc - T2$tic

# Compare times
t1/t2

# Flavors of vector multiplication
x <- c(1,2,3)
y <- c(4,5,6)
z <- c(7,8,9,10)
x*y         # element-wise
t(x) %*% y  # dot product
sum(x*y)    # dot product
x %*% t(z)  # outer product
x %o% z     # outer product
outer(x,z)  # outer product

# Taylor Polynomial
P1 <- function(x) {x}
P3 <- function(x) {x-x^3/6}
x <- seq(from = -pi/2, to = pi/2, length = 200)
plot(x, sin(x), col = "black", type = "l")
lines(x, P1(x), col = "green")
lines(x, P3(x), col = "red")