**Remember to show your work / give your reasoning, as applicable.**

### Problem 1

Which of the following sets of vectors are bases for ℝ<sup>2</sup>?

-   {(0,1), (1,1)}
-   {(1,0), (0,1), (1,1)}
-   {(1,0), (−1,0)}
-   {(1,1), (1,−1)}
-   {(1,1), (2,2)}
-   {(1,2)}

### Problem 1 solution

Your solution goes here.

### Problem 2

What is the dimension of the intersection of the following two planes in
ℝ<sup>3</sup>?

*x* + 2*y* − *z* = 0,   3*x* − 3*y* + *z* = 0

### Problem 2 Solution

Your solution goes here.

### Problem 3

Solve (by hand) the system below.

$$
\begin{pmatrix}1 & 2 & 0 \\\\ 3 & 2 & 4 \\\\ -2 & 1 & -2 \end{pmatrix} \begin{pmatrix} x \\\\ y \\\\ z \end{pmatrix} = \begin{pmatrix} 1 \\\\ 7\\\\ -1 \end{pmatrix}
$$

### Problem 3 Solution

Your solution goes here.

### Problem 4

Put the augmented coefficient matrix for the system of equations

$$
\begin{aligned}
x+y+z&=2\\\\
x+3y+3z&=0\\\\
x+3y+6z&=3
\end{aligned}
$$

into row echelon form.

### Problem 4 Solution

Your solution goes here.

### Problem 5

I’ve written a function called `eliminate` that performs Gaussian
elimination in order to transform a matrix to row echelon form.

    eliminate <- function(A, tol = 10^-8) {
      n <- nrow(A)
      for ( j in 1:(n-1) ) {
        pivot <- A[j,j]
        if (abs(pivot) < tol) stop('zero pivot encountered')
        for ( i in (j+1):n ) {
          A[i,] <- A[i,] - A[i,j]/pivot * A[j,]
        }
      }
      return(A)
    }

Demonstrate 𝒪(*n*<sup>3</sup>) complexity by doing the following. Run
the `eliminate` command on an *n* × *n* random matrix (generated using
the `runif` command) for each of *n* = 270, 540, 810, 1080, 1350, saving
the run time for each value of *n*. To save the run time, you can use
something like

    library(tictoc)
    tic()
    2+2 # Just as an example
    T1 <- toc()
    t1 <- T1$toc - T1$tic

Then, use the `lm` command to fit a line to the points (log*n*,log*t*)
and find the slope.

### Problem 5 solution

Your solution goes here.

### Problem 6

Compute by hand the eigenvalues of the matrices

$$
\mathbf{A}=\begin{pmatrix} 1 & 1000 \\\\ 0 & 1\end{pmatrix}, \qquad \widetilde{\mathbf{A}}=\begin{pmatrix} 1 & 1000 \\\\ 0.001 & 1\end{pmatrix}.
$$

Would you say that the problem of computing eigenvalues is
well-conditioned or ill-conditioned? Check the condition number of **A**
using the `kappa` command with the option `exact = TRUE`. For
concreteness, use the 2-norm, though there is nothing special about that
choice.

### Problem 6 Solution

Your solution goes here.
