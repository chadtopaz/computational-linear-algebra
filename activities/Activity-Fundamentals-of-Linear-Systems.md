Activity - Fundamentals of Linear Systems
================

Remember to show your work / give your reasoning, as applicable.

### Problem 1

Which of the following sets of vectors are bases for
![\\mathbb{R}^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BR%7D%5E2 "\mathbb{R}^2")?

1.  ![\\{(0, 1), (1, 1)\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7B%280%2C%201%29%2C%20%281%2C%201%29%5C%7D "\{(0, 1), (1, 1)\}")
2.  ![\\{(1, 0), (0, 1), (1, 1)\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7B%281%2C%200%29%2C%20%280%2C%201%29%2C%20%281%2C%201%29%5C%7D "\{(1, 0), (0, 1), (1, 1)\}")
3.  ![\\{(1, 0), (−1, 0\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7B%281%2C%200%29%2C%20%28%E2%88%921%2C%200%5C%7D "\{(1, 0), (−1, 0\}")
4.  ![\\{(1, 1), (1, −1)\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7B%281%2C%201%29%2C%20%281%2C%20%E2%88%921%29%5C%7D "\{(1, 1), (1, −1)\}")
5.  ![\\{((1, 1), (2, 2)\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7B%28%281%2C%201%29%2C%20%282%2C%202%29%5C%7D "\{((1, 1), (2, 2)\}")
6.  ![\\{(1, 2)\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7B%281%2C%202%29%5C%7D "\{(1, 2)\}")

### Problem 1 solution

Your solution goes here.

### Problem 2

What is the dimension of the intersection of the following two planes in
![\\mathbb{R}^3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BR%7D%5E3 "\mathbb{R}^3")?

![
x + 2y − z = 0, \\quad 3x − 3y + z = 0
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Ax%20%2B%202y%20%E2%88%92%20z%20%3D%200%2C%20%5Cquad%203x%20%E2%88%92%203y%20%2B%20z%20%3D%200%0A "
x + 2y − z = 0, \quad 3x − 3y + z = 0
")

### Problem 2 Solution

Your solution goes here.

### Problem 3

Solve (by hand) the system below.

![
\\begin{pmatrix}1 & 2 & 0 \\\\ 3 & 2 & 4 \\\\ -2 & 1 & -2 \\end{pmatrix} \\begin{pmatrix} x \\\\ y \\\\ z \\end{pmatrix} = \\begin{pmatrix} 1 \\\\ 7\\\\ -1 \\end{pmatrix}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Bpmatrix%7D1%20%26%202%20%26%200%20%5C%5C%203%20%26%202%20%26%204%20%5C%5C%20-2%20%26%201%20%26%20-2%20%5Cend%7Bpmatrix%7D%20%5Cbegin%7Bpmatrix%7D%20x%20%5C%5C%20y%20%5C%5C%20z%20%5Cend%7Bpmatrix%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%201%20%5C%5C%207%5C%5C%20-1%20%5Cend%7Bpmatrix%7D%0A "
\begin{pmatrix}1 & 2 & 0 \\ 3 & 2 & 4 \\ -2 & 1 & -2 \end{pmatrix} \begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} 1 \\ 7\\ -1 \end{pmatrix}
")

### Problem 3 Solution

Your solution goes here.

### Problem 4

Put the augmented coefficient matrix for the system of equations

![
\\begin{aligned}
x+y+z&=2\\\\
x+3y+3z&=0\\\\
x+3y+6z&=3
\\end{aligned}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Baligned%7D%0Ax%2By%2Bz%26%3D2%5C%5C%0Ax%2B3y%2B3z%26%3D0%5C%5C%0Ax%2B3y%2B6z%26%3D3%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
x+y+z&=2\\
x+3y+3z&=0\\
x+3y+6z&=3
\end{aligned}
")

into row echelon form.

### Problem 4 Solution

Your solution goes here.

### Problem 5

I’ve written a function called `eliminate` that performs Gaussian
elimination in order to transform a matrix to row echelon form.

``` r
eliminate <- function(A, tol=10^-8) {
  n = nrow(A)
  for ( j in 1:(n-1) ) {
    pivot = A[j,j]
    if (abs(pivot) < tol) stop('zero pivot encountered')
    for ( i in (j+1):n ) {
      A[i,] = A[i,] - A[i,j]/pivot * A[j,]
    }
  }
  return(A)
}
```

Demonstrate
![\\mathcal{O}(n^3)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathcal%7BO%7D%28n%5E3%29 "\mathcal{O}(n^3)")
complexity by doing the following. Run the `eliminate` command on an
![n \\times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%20n "n \times n")
random matrix for each of
![n=270,\\ 540,\\ 810,\\ 1080,\\ 1350](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%3D270%2C%5C%20540%2C%5C%20810%2C%5C%201080%2C%5C%201350 "n=270,\ 540,\ 810,\ 1080,\ 1350"),
saving the run time for each value of
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n").
Then, use the `lm` command to fit a line to the points
![(\\log n,\\log t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28%5Clog%20n%2C%5Clog%20t%29 "(\log n,\log t)")
and find the slope.

### Problem 5 solution

Your solution goes here.

### Problem 6

Compute by hand the eigenvalues of the matrices

![
\\mathbf{A}=\\begin{pmatrix} 1 & 1000 \\\\ 0 & 1\\end{pmatrix}, \\qquad \\widetilde{\\mathbf{A}}=\\begin{pmatrix} 1 & 1000 \\\\ 0.001 & 1\\end{pmatrix}.
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cmathbf%7BA%7D%3D%5Cbegin%7Bpmatrix%7D%201%20%26%201000%20%5C%5C%200%20%26%201%5Cend%7Bpmatrix%7D%2C%20%5Cqquad%20%5Cwidetilde%7B%5Cmathbf%7BA%7D%7D%3D%5Cbegin%7Bpmatrix%7D%201%20%26%201000%20%5C%5C%200.001%20%26%201%5Cend%7Bpmatrix%7D.%0A "
\mathbf{A}=\begin{pmatrix} 1 & 1000 \\ 0 & 1\end{pmatrix}, \qquad \widetilde{\mathbf{A}}=\begin{pmatrix} 1 & 1000 \\ 0.001 & 1\end{pmatrix}.
")

Would you say that the problem of computing eigenvalues is
well-conditioned or ill-conditioned? Check the condition number of
![\\mathbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D "\mathbf{A}")
using the `kappa` command. For concreteness, use the 2-norm, though
there is nothing special about that choice.

### Problem 6 Solution

Your solution goes here.
