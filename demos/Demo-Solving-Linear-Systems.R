# Load libraries
library(pracma)

# LU decomposition
A <- matrix(c(1,3,2,4,6,5,9,8,7), nrow = 3, byrow = TRUE)
LU <- lu(A, scheme = "ijk")
L <- LU$L
U <- LU$U
L
U
(L %*% U) - A

# LU solve with forward solve / back solve
b <- c(1,3,5)
y <- forwardsolve(L,b)
x1 <- backsolve(U,y)
x2 <- solve(A,b)
Norm(x1 - x2, p = Inf)

# Fixed point iteration
x <- -2
for (i in 1:20){
  x <-  3/(x-2)
}
x <- 3.0001
for (i in 1:20){
  x <-  3/(x-2)
}

# Jacobi iteration
A <- matrix(c(3,1,-1,2,6,-3,3,5,-10),nrow = 3, byrow = TRUE)
b <- c(1,2,3)
x1 <- solve(A,b)
D <- diag(A)
R <- A - diag(D)
Dinv <- 1/D
x2 <- c(0,0,0)
for (i in 1:30) { 
  x2 <- (b - R%*%x2)*Dinv
}
Norm(x1 - x2, p = Inf)