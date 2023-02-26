Activity - Solving Linear Systems
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

a\. By hand (without using a computer) find the **LU** decomposition of
the matrix

$$
\mathbf{A} = \begin{pmatrix}
4 & 5 & 6\\
8 & 14 & 17\\
12 & 23 & 32
\end{pmatrix}
$$

b\. Check your answer by entering your **L** and **U** matrices into R
and carrying out matrix multiplication (with `%*%`) to recover the
original matrix.

c\. Use the built in command `lu` with option `scheme = "ijk"` to obtain
the same result as you did by hand.

### Problem 1 Solution

a\. Your solution goes here.

b\. Your solution goes here.

c\. Your solution goes here.

### Problem 2

a\. Generate a 1000 x 1000 random matrix **A** with entries drawn
uniformly from the interval \[0,1\]. Use R’s build in `solve` command to
solve **Ax**=**b** for 100 randomly generated **b** vectors that are
1000 x 1 vectors with random entries drawn in the same manner as for A.
Time the calculation.

b\. Using the same matrix **A**, compute (a single time) its **LU**
decomposition and then use the `backsolve` and `forwardsolve` commands
to solve **LUx**=**b** for 100 randomly generated **b** vectors. Time
the calculation. When you time it, for fairness, you should include the
time to calculate the single **LU** decomposition as well as the 100
solutions for different **b**.

c\. Compare the times.

### Problem 2 Solution

a\. Your solution goes here.

b\. Your solution goes here.

c\. Your solution goes here.

### Problem 3

Consider **Ax**=**b** for

$$
\mathbf{A} = \begin{pmatrix}
4 & -1 & -1 \\
-2 & 6 & 1\\
-1 & 1 & 7
\end{pmatrix}, \quad
\mathbf{b} = \begin{pmatrix}
3 \\ 9 \\ -6
\end{pmatrix}
$$

a\. Show that **A** is strictly diagonally dominant.

b\. Choose any initial guess for **x** that you want (really! or you can
try different ones!) and use Jacobi iteration to solve until the
backwards error is less than $10^{-6}$ in the infinity norm. Report your
answer.

### Problem 3 Solution

a\. Your solution goes here.

b\. Your solution goes here.
