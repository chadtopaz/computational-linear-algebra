# Load libraries
library(tidyverse)

# Define and plot data
u <- c(2.5, 0.5, 2.2, 1.9, 3.1, 2.3, 2, 1, 1.5, 1.1)
v <- c(2.4, 0.7, 2.9, 2.2, 3.0, 2.7, 1.6, 1.1, 1.6, 0.9)
plot(u, v, pch = 19)

# Check means and adjust so that data has mean zero
mean(u)
mean(v)
u <- u - mean(u)
v <- v - mean(v)
summary(u)
summary(v)
plot(u, v, pch = 19)

# Create covariance matrix
m <- matrix(c(cov(u, u), cov(u, v), cov(v, u),cov(v, v)),nrow=2,ncol=2,byrow=TRUE,dimnames=list(c("u","v"),c("u","v")))
# Note: m is symmetric and positive definite
m
# This means all eigenvalues are positive and
# Eigenbasis is orthogonal

# Calculate eigenvectors and eigenvalues
e <- eigen(m)
v1 <- e$vectors[,1]
v2 <- e$vectors[,2]

# Show that they are orthogonal
sum(v1*v2)

# Plot eigendirections
abline(0,v1[2]/v1[1],col="red")
abline(0,v2[2]/v2[1],col="blue")

# See which direction is most important
lambda1 <- e$values[1]
lambda2 <- e$values[2]
lambda1/(lambda1+lambda2)
lambda2/(lambda1+lambda2)

# Express data in the basis of v1, v2
V <- cbind(v1,v2)
solve(V) %*% rbind(u,v)
# However, note that solve(V) = t(v), so
t(V) %*% rbind(u,v)
# We usually store this in columns, so I will transpose the whole thing
x <- t(t(V) %*% rbind(u,v))

# Let's keep just first principal component
# That means use the left column
# and interpret it as amount of first eigenvector to take
lowdim <- x[,1] %*% t(V[,1])
points(lowdim,col="red",pch=19)

# Now do it all using prcomp
# prcomp expects observations to be ROWS
U <- cbind(u,v)
pca <- prcomp(U)

# Look at amount of variance explained
summary(pca)
lambda1/(lambda1+lambda2)
lambda2/(lambda1+lambda2)
summary(pca)
importance <- summary(pca)$importance
barplot(importance[3,])

# Look at principal components
pca$x
x

# Look at eigenvectors
pca$rotation
V

# Now an actual cool example
worldbank <- read.csv("worldbank.csv")
X <- worldbank[,2:7]
pca <- prcomp(X)
summary(pca)
plot(pca$x[,1:2],pch=19)
data <- data.frame(country=worldbank$country,x=pca$x[,1],y=pca$x[,2])
keeps <- c("Bolivia","Venezuela, RB","Mexico","United States","Canada","England","France","Germany","Russian Federation","Turkey","Afghanistan","Saudi Arabia","India","Myanmar","China","South Korea","Japan","Russia","New Zealand")
Europe <- c("Albania","Armenia","Austria","Azerbaijan","Belarus","Belgium","Bulgaria","Croatia","Cyprus","Denmark","Estonia","Finland","France")
ggplot(subset(data,country %in% union(keeps,Europe)),aes(x=x,y=y,label=country)) + geom_point() + geom_label()

set.seed(2)
mysample <- sample(1:nrow(data),20)
ggplot(data[mysample,],aes(x=x,y=y,label=country)) + geom_point() + geom_label()
