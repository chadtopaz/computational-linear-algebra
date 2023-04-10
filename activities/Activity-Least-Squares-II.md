Activity - Least Squares II
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

a\. Use the classical Gram-Schmidt orthogonalization algorithm to find
the reduced $\mathbf{QR}$ factorization and full $\mathbf{QR}$
factorization of the matrix

$$
\mathbf{A} = \begin{pmatrix}
2 & 3 \\
-2 & -6 \\
1 & 0
\end{pmatrix}.
$$

Don’t use R’s built-in `qr` function. Do the calculation by hand.
However, you can certainly use R to carry out the algebra if you want.

b\. First check your answers by verifying that
$\mathbf{A}=\mathbf{Q}\mathbf{R}$ and
$\mathbf{A}=\overline{\mathbf{Q}}\overline{\mathbf{R}}$. Then check
whether you computed the same factorization as the `qr` function, which
you can now use.

Note that the `qr` algorithm does not ensure that all of the diagonal
entries of $\mathbf{R}$ are nonnegative (in which case the factorization
is not unique). If you want to force the $\mathbf{R}$ matrix to have
positive diagonals, you can form a diagonal matrix $\mathbf{S}$ whose
$i^{th}$ diagonal is equal to the sign of the $i^{th}$ diagonal of
$\mathbf{R}$. Then let $\widetilde{\mathbf{Q}}=\mathbf{Q}\mathbf{S}$ and
$\widetilde{\mathbf{R}}=\mathbf{S}\mathbf{R}$, so that
$\widetilde{\mathbf{Q}}\widetilde{\mathbf{R}}=\mathbf{Q}\mathbf{S}^2 \mathbf{R}=\mathbf{Q}\mathbf{R}=\mathbf{A}$
(since $\mathbf{S}^2=\mathbf{I}$).

c\. Use the reduced $\mathbf{QR}$ factorization of $\mathbf{A}$ from
part (a) to find the least squares solution to

$$
\begin{pmatrix}
2 & 3 \\
-2 & -6 \\
1 & 0
\end{pmatrix}
\begin{pmatrix}
x_1 \\
x_2
\end{pmatrix} = \begin{pmatrix}
3 \\
-3\\
6
\end{pmatrix}.
$$

d\. Now solve the same system using R’s built-in `qr.solve` function.
You’ll get the same answer as in part (c) above. The point of this
question is just to show you a new, useful R command that shortcuts some
steps for you.

e\. If $\mathbf{A}$ is an $m\times n$ matrix, the null space of
$\mathbf{A}^{\top}$ is the orthogonal complement of the column space of
$\mathbf{A}$. Use the full $\mathbf{QR}$ factorization of $\mathbf{A}$
above to find a basis for the null space of $\mathbf{A}^{\top}$.

### Problem 1 Solution

a\. Your solution goes here.

b\. Your solution goes here.

c\. Your solution goes here.

d\. Your solution goes here.

e\. Your solution goes here.
