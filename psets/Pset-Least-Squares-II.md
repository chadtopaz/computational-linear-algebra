Pset - Least Squares II
================

------------------------------------------------------------------------

When you begin work during class, work with your assigned partner.
Please have only one electronic device open and work on it jointly. When
writing up this assignment, please remember that showing all of your
work and giving your reasoning are critical parts of achieving mastery.
If the course staff cannot tell how you solved a problem or finds leaps
in explanation or logic, the problem is not mastered. Finally, as a
matter of academic integrity, please make sure that you are positioned
to honestly answer yes to these questions:

- Have I disclosed everyone with whom I collaborated on this work? (Even
  if it is only my assigned partner.)

- Have I made a substantive intellectual contribution to the solution of
  every problem?

- Am I making sure not to pass off as my own work any work that belongs
  to someone else?

Whether intentional or unintentional, any potential violations of
academic integrity will be referred to the Honor Committee.

------------------------------------------------------------------------

Load necessary packages:

``` r
library(pracma)
```

### Problem 1

Define

$$
\mathbf{A} = \begin{pmatrix}
1 & 0 & 5\\
2 & 1 & -5 \\
3 & 3 & 1 \\
4 & 0 & 0
\end{pmatrix}.
$$

Find the matrix $\mathbf{R}$ that reflects a given vector $\mathbf{b}$
across the column space of $\mathbf{A}$. Express it both symbolically in
terms of $\mathbf{A}$ and $\mathbf{R}$, and n numerically (by printing
out the matrix). You don’t need to include a drawing as part of your
answer, but I urge you to make yourself a drawing and use vector/matrix
arithmetic to help you solve this problem.

### Problem 1 Solution

Your solution goes here.

### Problem 2

Define the $d$-degree polynomial

$$
f(x) = 1 + x^2 + \ldots + x^d,
$$

that is, the $d$-degree polynomial with all coefficients equal to one.
Let $(x_1, y_1), \ldots, (x_9, y_9)$ be 9 points sampled from this
function with $x$ values equally spaced on $[2,4]$. Assume $d < 8$. In
that case, we (theoretically) have more than enough points to recover
the coefficients in the polynomial by plugging in the values to a
generic $d$-degree polynomial and solving for the unknown
coefficients.This is a least squares problem. For each of $d = 4$,
$d = 5$, and $d = 6$, try to solve the least squares problem two ways:
(1) by setting up the normal equations and solving them with the `solve`
command, and (2) by using the `qr.solve` command. Compare the results of
the two methods to each other, and also, across the different values of
$d$.

### Problem 2 Solution

Your solution goes here.
