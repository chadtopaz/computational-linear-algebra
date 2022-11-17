# Load libraries
library(pracma)
library(expm)

# What is power iteration?
A <- matrix(c(1,-12,-12,-6), nrow = 2)
x <- c(1,1)
for (i in 1:100){
  x <- A%*%x
  x <- x/Norm(x)
}
print(x)

# How do we understand the answer?
sol <- eigen(A)
Lambda <- diag(sol$values)
S <- sol$vectors
S%*%Lambda%*%solve(S)
v1 <- S[,1]
v2 <- S[,2]
x <- c(1,1)
solve(S)%*%x
-0.2*v1 + 1.4*v2

# Page Rank
A <- matrix(c(0,1,0,0,1/2,0,1/2,0,0,0,0,1,0,1,0,0), nrow = 4)
(A%^%100) %*% rep(0.25,4)
v <- eigen(A)$vectors[,1]
v <- v/sum(v)

# Absorbing Markov chain example

# Define transition matrix
P <- matrix(c(0,0,0.2,0,0,0.2,1,0.5,0,0,0.4,0,0,0,0.6,0.3,0,0.2,1,0.4,0.1,0,0.1,0,0), nrow = 5, byrow = TRUE)
rownames(P) <- 1:5
colnames(P) <- 1:5

# Plot Markov chain diagram
A <- t(P)
G <- graph_from_adjacency_matrix(A, weighted = TRUE)
plot(G, edge.width = A[as_edgelist(G)])

# Reorder the states
absorbing <- c(2,4)
transient <- c(1,3,5)
P[c(absorbing, transient), c(absorbing, transient)]

# Extract/create submatrices for computation
S <- P[absorbing, transient]
R <- P[transient, transient]
Itilde <- diag(nrow(R))

# Compute absorption probabilities
S %*%solve(Itilde - R)