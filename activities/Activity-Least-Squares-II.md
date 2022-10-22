Remember to show your work / give your reasoning, as applicable.

    library(pracma)

### Problem 1

a\. Use the classical Gram-Schmidt orthogonalization algorithm to find
the reduced **QR** factorization and full **QR** factorization of the
matrix

$$
\mathbf{A} = \begin{pmatrix}
2 & 3 \\\\
-2 & -6 \\\\
1 & 0
\end{pmatrix}.
$$

Don’t use R’s built-in `qr` function. Do the calculation by hand.
However, you can certainly use R to carry out the algebra if you want.

b\. First check your answers by verifying that **A** = **Q****R** and
$\mathbf{A}=\overline{\mathbf{Q}}\overline{\mathbf{R}}$. Then check
whether you computed the same factorization as the `qr` function, which
you can now use.

Note that the `qr` algorithm does not ensure that all of the diagonal
entries of **R** are nonnegative (in which case the factorization is not
unique). If you want to force the **R** matrix to have positive
diagonals, you can form a diagonal matrix **S** whose
*i*<sup>*t**h*</sup> diagonal is equal to the sign of the
*i*<sup>*t**h*</sup> diagonal of **R**. Then let
$\widetilde{\mathbf{Q}}=\mathbf{Q}\mathbf{S}$ and
$\widetilde{\mathbf{R}}=\mathbf{S}\mathbf{R}$, so that
$\widetilde{\mathbf{Q}}\widetilde{\mathbf{R}}=\mathbf{Q}\mathbf{S}^2 \mathbf{R}=\mathbf{Q}\mathbf{R}=\mathbf{A}$
(since **S**<sup>2</sup> = **I**).

c\. Use the reduced **QR** factorization of **A** from part (a) to find
the least squares solution to

$$
\begin{pmatrix}
2 & 3 \\\\
-2 & -6 \\\\
1 & 0
\end{pmatrix}
\begin{pmatrix}
x\_1 \\\\
x\_2
\end{pmatrix}
=
\begin{pmatrix}
3 \\\\
-3\\\\
6
\end{pmatrix}.
$$

d\. Now solve the same system using R’s built-in `qr.solve` function.
You’ll get the same answer as in part (c) above. The point of this
question is just to show you a new, useful R command that shortcuts some
steps for you.

e\. If **A** is an *m* × *n* matrix, the null space of **A**<sup>⊤</sup>
is the orthogonal complement of the column space of **A**. Use the full
**QR** factorization of **A** above to find a basis for the null space
of **A**<sup>⊤</sup>.

### Problem 1 Solution

a\. Your solution goes here.

b\. Your solution goes here.

c\. Your solution goes here.

d\. Your solution goes here.

e\. Your solution goes here.
