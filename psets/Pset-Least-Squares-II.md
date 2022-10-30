Load necessary packages:

    library(pracma)

### Problem 1

Define

$$
\mathbf{A} = \begin{pmatrix} 1 & 0 & 5\\\\ 2 & 1 & -5 \\\\ 3 & 3 & 1 \\\\ 4 & 0 & 0 \end{pmatrix}.
$$

Find the matrix **R** that reflects a given vector **b** across the
column space of **A**. Express it both symbolically in terms of **A**
and **R**, and n numerically (by printing out the matrix). You don’t
need to include a drawing as part of your answer, but I urge you to make
yourself a drawing and use vector/matrix arithmetic to help you solve
this problem.

### Problem 1 Solution

Your solution goes here.

### Problem 2

Define the *d*-degree polynomial

*f*(*x*) = 1 + *x*<sup>2</sup> + … + *x*<sup>*d*</sup>,
that is, the *d*-degree polynomial with all coefficients equal to one.
Let
(*x*<sub>1</sub>,*y*<sub>1</sub>), …, (*x*<sub>9</sub>,*y*<sub>9</sub>)
be 9 points sampled from this function with *x* values equally spaced on
\[2,4\]. Assume *d* &lt; 8. In that case, we (theoretically) have more
than enough points to recover the coefficients in the polynomial by
plugging in the values to a generic *d*-degree polynomial and solving
for the unknown coefficients.This is a least squares problem. For each
of *d* = 4, *d* = 5, and *d* = 6, try to solve the least squares problem
two ways: (1) by setting up the normal equations and solving them with
the `solve` command, and (2) by using the `qr.solve` command. Compare
the results of the two methods to each other, and also, across the
different values of *d*.

### Problem 2 Solution

Your solution goes here.
