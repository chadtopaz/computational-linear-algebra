# Gram-Schmidt example
# Note vectors are not linearly independent
# since 1/3(v1+v2) = v3
# So we should only generate two vectors
v1 <- c(1,2,2)
v2 <- c(2,1,-2)
v3 <- c(1,1,0)
y1 <- v1
r11 <- Norm(y1)
q1 <- y1/r11
y2 <- v2 - (q1%*%t(q1))%*%v2
r22 <- Norm(y2)
q2 <- y2/r22
q1
q2
Norm(q1)
Norm(q2)
dot(q1,q2)

# QR Example
# Define matrix
v1 <- c(2,3,6)
v2 <- c(10,8,9)
v3 <- c(19,-3,-13)
A <- cbind(v1,v2,v3)
A
# Step 1
y1 <- v1
r11 <- Norm(y1)
q1 <- y1/r11
# Step2
y2 <- v2 - (q1%*%t(q1))%*%v2
r22 <- Norm(y2)
q2 <- y2/r22
r12 <- dot(q1,v2)
# Step 3
y3 <- v3 - (q1%*%t(q1))%*%v3 - (q2%*%t(q2))%*%v3
r33 <- Norm(y3)
q3 <- y3/r33
r13 <- dot(q1,v3)
r23 <- dot(q2,v3)
# Assemble into Qbar and Rbar
Qbar <- cbind(q1,q2,q3)
Rbar <- rbind(c(r11,r12,r13),c(0,r22,r23),c(0,0,r33))
Qbar
Rbar
Qbar%*%Rbar - A
QRbarcheck <- qr(A)
Qbarcheck <- qr.Q(QRbarcheck, complete=TRUE)
Rbarcheck <- qr.R(QRbarcheck, complete=TRUE)
Qbar
Qbarcheck
S <- diag(c(-1,-1,-1))
Qbarcheck <- Qbarcheck%*%S
Rbarcheck <- S%*%Rbarcheck
Qbarcheck - Qbar
Rbarcheck - Rbar

# Next QR example
# Define matrix
v1 <- c(1,2,2)
v2 <- c(2,1,-2)
v3 <- c(1,1,0)
A <- cbind(v1,v2,v3)
# Step 1
y1 <- v1
r11 <- Norm(y1)
q1 <- y1/r11
# Step2
y2 <- v2 - (q1%*%t(q1))%*%v2
r22 <- Norm(y2)
q2 <- y2/r22
r12 <- dot(q1,v2)
# Step 3
y3 <- v3 - (q1%*%t(q1))%*%v3 - (q2%*%t(q2))%*%v3
y3
V3 <- c(1,1,1)
y3 <- V3 - (q1%*%t(q1))%*%V3 - (q2%*%t(q2))%*%V3
q3 <- y3/Norm(y3) 
r13 <- dot(q1,v3)
r23 <- dot(q2,v3)
r33 <- 0
# Assemble into Qbar and Rbar
Qbar <- cbind(q1,q2,q3)
Rbar <- rbind(c(r11,r12,r13),c(0,r22,r23),c(0,0,r33))
Qbar
Rbar
Qbar%*%Rbar - A
Q <- Qbar[1:3,1:2]
R <- Rbar[1:2,1:3]
Q
R
Q%*%R - A

# Last QR example
v1 <- c(1,1,1,1)
v2 <- c(0,1,1,1)
v3 <- c(0,0,1,1)
A <- cbind(v1,v2,v3)
A
# Step 1
y1 <- v1
r11 <- Norm(y1)
q1 <- y1/r11
# Step 2
y2 <- v2 - (q1%*%t(q1))%*%v2
r22 <- Norm(y2)
q2 <- y2/r22
r12 <- dot(q1,v2)
# Step 3
y3 <- v3 - (q1%*%t(q1))%*%v3 - (q2%*%t(q2))%*%v3
r33 <- Norm(y3) 
q3 <- y3/(r33)
r13 <- dot(q1,v3)
r23 <- dot(q2,v3)
# Choose a vector not in the span of q1, q2, q3
v4 <- c(1,2,3,4)
y4 <- v4 - (q1%*%t(q1))%*%v4 - (q2%*%t(q2))%*%v4  - (q3%*%t(q3))%*%v4
r44 <- Norm(y4) 
q4 <- y4/(r44)
r14 <- dot(q1,v4)
r24 <- dot(q2,v4)
r34 <- dot(q3,v4)
# Assemble into Qbar and Rbar, check answer
Qbar <- cbind(q1,q2,q3,q4)
Rbar <- rbind(c(r11,r12,r13),c(0,r22,r23),c(0,0,r33),c(0,0,0))
Qbar
Rbar
Qbar%*%Rbar - A
# Make reduced QR and check
Q <- Qbar[1:4,1:3]
R <- Rbar[1:3,1:3]
Q%*%R - A