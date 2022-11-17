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