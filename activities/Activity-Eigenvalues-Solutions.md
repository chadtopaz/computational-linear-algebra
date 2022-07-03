Activity - Eigenvalues
================
Solutions

### Procedure

-   Please collaborate on this activity with one person as scribe
    (typist). There should only be one laptop open. We’re aiming for
    actual collaboration here, not two or three people working in
    parallel.

-   It doesn’t matter if you finish the whole activity during class.

-   It **does** mater if you focus on helping each other learn. Helping
    each other learn is the point of this activity. Completing the
    problems is merely a mechanism to guide you towards learning.

-   It **does** mater if everyone participates. The expectation is that
    every group member participates, where the definition of
    “participation” can include asking questions, answering questions,
    brainstorming, or any other form of meaningful engagement.

-   In the last three minutes of class, before you leave, please knit
    your .Rmd file to html, print the html to .pdf, and have the scribe
    submit it as a group assignment through GLOW.

-   Don’t forget that the solutions are available on GLOW.

### Problem 0

Recalling that academic integrity policies apply even in an ungraded
course, please disclose any irregular circumstances related to
participation in your group (for instance, a group member was missing
and did not contribute to this document).

### Problem 0 Solution

Your solution goes here.

### Problem 1

Let

![
\\mathbf{A}=\\frac{1}{9}\\begin{pmatrix} \~\~\~83 & \~\~\~296 & -128 \\\\ \~\~\~296 & \~\~\~473 & -152 \\\\ -128 & -152 & \~\~\~335\\end{pmatrix}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cmathbf%7BA%7D%3D%5Cfrac%7B1%7D%7B9%7D%5Cbegin%7Bpmatrix%7D%20~~~83%20%26%20~~~296%20%26%20-128%20%5C%5C%20~~~296%20%26%20~~~473%20%26%20-152%20%5C%5C%20-128%20%26%20-152%20%26%20~~~335%5Cend%7Bpmatrix%7D%0A "
\mathbf{A}=\frac{1}{9}\begin{pmatrix} ~~~83 & ~~~296 & -128 \\ ~~~296 & ~~~473 & -152 \\ -128 & -152 & ~~~335\end{pmatrix}
")

Using any vector
![\\mathbf{v}\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_0 "\mathbf{v}_0")
in
![\\mathbb{R}^3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BR%7D%5E3 "\mathbb{R}^3")
as a starting guess, perform power iteration on **A** to find the
doiminant eigenvector
![\\mathbf{v}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D "\mathbf{v}").
Stop computing when
![\|\|\\mathbf{v}\_{i} - \\mathbf{v}\_{i-1}\|\|\_2 \< 0.5 \\times 10^{-6}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%7C%5Cmathbf%7Bv%7D_%7Bi%7D%20-%20%5Cmathbf%7Bv%7D_%7Bi-1%7D%7C%7C_2%20%3C%200.5%20%5Ctimes%2010%5E%7B-6%7D "||\mathbf{v}_{i} - \mathbf{v}_{i-1}||_2 < 0.5 \times 10^{-6}").
Then estimate the largest (in magnitude) eigenvalue
![\\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda")
as follows. Take the definition of an eigenpair, left multiply each side
by
![\\mathbf{v}^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D%5ET "\mathbf{v}^T"),
and solve for
![\\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda"),
that is:

![
\\begin{eqnarray\*}
\\mathbf{A} \\mathbf{v} & = & \\lambda \\mathbf{v}\\\\
\\mathbf{v}^T \\mathbf{A} \\mathbf{v} & = & \\lambda \\mathbf{v}^T \\mathbf{v} \\\\
\\lambda & = & \\frac{\\mathbf{v}^T \\mathbf{A} \\mathbf{v}}{\\mathbf{v}^T \\mathbf{v}}.
\\end{eqnarray\*}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Beqnarray%2A%7D%0A%5Cmathbf%7BA%7D%20%5Cmathbf%7Bv%7D%20%26%20%3D%20%26%20%5Clambda%20%5Cmathbf%7Bv%7D%5C%5C%0A%5Cmathbf%7Bv%7D%5ET%20%5Cmathbf%7BA%7D%20%5Cmathbf%7Bv%7D%20%26%20%3D%20%26%20%5Clambda%20%5Cmathbf%7Bv%7D%5ET%20%5Cmathbf%7Bv%7D%20%5C%5C%0A%5Clambda%20%26%20%3D%20%26%20%5Cfrac%7B%5Cmathbf%7Bv%7D%5ET%20%5Cmathbf%7BA%7D%20%5Cmathbf%7Bv%7D%7D%7B%5Cmathbf%7Bv%7D%5ET%20%5Cmathbf%7Bv%7D%7D.%0A%5Cend%7Beqnarray%2A%7D%0A "
\begin{eqnarray*}
\mathbf{A} \mathbf{v} & = & \lambda \mathbf{v}\\
\mathbf{v}^T \mathbf{A} \mathbf{v} & = & \lambda \mathbf{v}^T \mathbf{v} \\
\lambda & = & \frac{\mathbf{v}^T \mathbf{A} \mathbf{v}}{\mathbf{v}^T \mathbf{v}}.
\end{eqnarray*}
")

### Problem 1 Solution

``` r
A <- 1/9*cbind(c(83,296,-128),c(296,473,-152),c(-128,-152,335))
PI <- function(A,tol=.000005) {
  err <- Inf
  v <- runif(ncol(A))
  v <- v/Norm(v)
  while (err > tol) {
    vold <- v
    v <- A %*% v
    v <- v/Norm(v)
    err <- Norm(v-vold)
  }
  return(v)
}
print(v <- PI(A))
```

    ##            [,1]
    ## [1,]  0.4444446
    ## [2,]  0.7777784
    ## [3,] -0.4444432

``` r
print(lambda <- (t(v)%*%A%*%v)/(t(v)%*%v))
```

    ##      [,1]
    ## [1,]   81

### Problem 2

Every position in the United States House of Representatives is up for
election every two years. A linear algebra student built the following
(admittedly oversimplified, but still fun) model for predicting the
political party of the winning candidate in a local congressional
district:

-   If our current representative is a Democrat, there is a 75% chance
    that a Democrat will win the next election, and a 25% chance that a
    Republican will win the next election.

-   If our current representative is a Republican, there is a 50% chance
    that a Democrat will win the next election, and a 50% chance that a
    Republican will win the next election.

Note that this model assumes that the winner will always be a Democrat
or Republican.

1.  Let’s model the political party of the current representative as a
    Markov chain with two states. Write down the transition matrix
    ![\\mathbf{P}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BP%7D "\mathbf{P}").
    Let the first state be Democrat.

2.  We can represent the Markov chain as the dynamical system
    ![\\mathbf{x}\_{k+1}=Ax_k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bx%7D_%7Bk%2B1%7D%3DAx_k "\mathbf{x}_{k+1}=Ax_k"),
    where  
    ![x_k \\in \\mathbb{R}^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_k%20%5Cin%20%5Cmathbb%7BR%7D%5E2 "x_k \in \mathbb{R}^2")
    is a vector whose first component represents the probability that a
    Democrat wins the
    ![k^{th}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%5E%7Bth%7D "k^{th}")
    election, and the second component represents the probability that a
    Republican wins the
    ![k^{th}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%5E%7Bth%7D "k^{th}")
    election. What is the matrix
    ![\\mathbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D "\mathbf{A}")
    so that
    ![x\_{k+1}=\\mathbf{A}x_k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_%7Bk%2B1%7D%3D%5Cmathbf%7BA%7Dx_k "x_{k+1}=\mathbf{A}x_k")?
    Hint: Double check your
    ![\\mathbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D "\mathbf{A}")
    on a simple example.

3.  Let
    ![k=0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%3D0 "k=0")
    be the year 2020 election. A Democrat won the election last fall, so
    we’ll take
    ![\\mathbf{x}\_0=(1,0)^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bx%7D_0%3D%281%2C0%29%5ET "\mathbf{x}_0=(1,0)^T").
    According to the model, what is the probability that a Democrat will
    win the election in 2026? Hint: to raise a matrix
    ![A](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A "A")
    to, for example, the power 2, make sure the `expm` package is
    installed, and type `A%^%2`. This function can’t handle sparse
    matrices, however, so if
    ![A](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A "A")
    were stored that way, you would need to write `as.matrix(A)%^%2`.

4.  Let’s introduce some helpful terminology. A matrix
    ![\\mathbf{M}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BM%7D "\mathbf{M}")
    is **primitive** if there exists a positive integer
    ![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k")
    such that
    ![\\mathbf{M}^k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BM%7D%5Ek "\mathbf{M}^k")
    has all positive entries.
    ![\\mathbf{M}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BM%7D "\mathbf{M}")
    is **stochastic** if its columns sum to one. There is a theorem
    called the **Perron-Frobenius Theorem** that says: if
    ![\\mathbf{M}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BM%7D "\mathbf{M}")
    is a nonnegative square primitive matrix, then there is a dominant
    eigenvalue
    ![\\lambda_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_1 "\lambda_1")
    with eigenvector
    ![\\mathbf{v}\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_1 "\mathbf{v}_1")
    that has all positive entries. If
    ![\\mathbf{M}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BM%7D "\mathbf{M}")
    is stochastic, then additionally,
    ![\\lambda_1 = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_1%20%3D%201 "\lambda_1 = 1"),
    so that
    ![\\mathbf{M}\\mathbf{v}\_1 = \\mathbf{v}\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BM%7D%5Cmathbf%7Bv%7D_1%20%3D%20%5Cmathbf%7Bv%7D_1 "\mathbf{M}\mathbf{v}_1 = \mathbf{v}_1").
    In the context of Markov chains, you can imagine
    ![\\mathbf{v}\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_1 "\mathbf{v}_1")
    as being the (possibly scaled) stationary distribution of
    probabilities, that is, the probabilities of being in each state
    after many, many iterations of the Markov chain. Now apply this
    concept. In the limit of very long time, what is the probability
    that a Democrat will win the election? Would your answer be
    different if a Republican had won the last election?

### Problem 2 Solution

1.  

``` r
P <- rbind(c(.75,.25),c(.5,.5))
```

2.  

``` r
A <- t(P)
```

After the current election, there is probability 1 that a democrat is in
office, so we can check the matrix
![\\mathbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D "\mathbf{A}")
by letting
![\\mathbf{x}\_0 = (1,0)^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bx%7D_0%20%3D%20%281%2C0%29%5ET "\mathbf{x}_0 = (1,0)^T")
and computing
![\\mathbf{A}\\mathbf{x}\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D%5Cmathbf%7Bx%7D_0 "\mathbf{A}\mathbf{x}_0"):

``` r
x0 <- c(1,0)
A%*%x0
```

    ##      [,1]
    ## [1,] 0.75
    ## [2,] 0.25

Indeed it shows there is a 75% chance of a democrat winning the next
election.

3.  The election in 2026 is the third election, so we want to examine
    ![\\mathbf{A}^3 \\mathbf{x}\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D%5E3%20%5Cmathbf%7Bx%7D_0 "\mathbf{A}^3 \mathbf{x}_0"):

``` r
A%^%3%*%x0
```

    ##          [,1]
    ## [1,] 0.671875
    ## [2,] 0.328125

So there is a 67.2% chance that a democrat wins the 2026 election
according to this model.

4.  Note that
    ![\\mathbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D "\mathbf{A}")
    is nonnegative. It is primitive because
    ![\\mathbf{A}^1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D%5E1 "\mathbf{A}^1")
    is a positive matrix. It is stochastic because its columns sum to 1:

``` r
colSums(A)
```

    ## [1] 1 1

To find the long-term probabilities, we can look at the dominant
eigenvector of
![\\mathbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D "\mathbf{A}").
Note that `R` normalizes this vector to be length 1, but we want to
normalize it to sum to one to get a probability distribution. Since our
matrix is small, feel free to use R’s built-in `eigen` command to find
the eigenvectors.

``` r
out <- eigen(A)
long.term.prob <- out$vectors[,1]/sum(out$vectors[,1])
print(long.term.prob)
```

    ## [1] 0.6666667 0.3333333

So there is an approximately 2/3 chance a democrat wins the election far
in the future. By the Perron-Frobenius theorem, our answer wouldn’t
depend on the initial state of the Markov chain.

### Problem 3

1.  Show that if
    ![\\mathbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D "\mathbf{A}")
    is an
    ![n \\times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%20n "n \times n")
    invertible matrix and
    ![\\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda")
    is an eigenvalue of
    ![\\mathbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D "\mathbf{A}")
    with eigenvector
    ![\\mathbf{v}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D "\mathbf{v}"),
    then
    ![1/\\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%2F%5Clambda "1/\lambda")
    is an eigenvalue of
    ![\\mathbf{A}^{-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D%5E%7B-1%7D "\mathbf{A}^{-1}")
    with the same eigenvector
    ![\\mathbf{v}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D "\mathbf{v}").

2.  What happens if you apply the power iteration to
    ![\\mathbf{A}^{-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D%5E%7B-1%7D "\mathbf{A}^{-1}")?

3.  Let
    ![\\mathbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D "\mathbf{A}")
    be an
    ![n \\times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%20n "n \times n")
    matrix, and let
    ![\\mathbf{C}=\\mathbf{A}-s\\mathbf{I}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BC%7D%3D%5Cmathbf%7BA%7D-s%5Cmathbf%7BI%7D "\mathbf{C}=\mathbf{A}-s\mathbf{I}"),
    where
    ![\\mathbf{I}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BI%7D "\mathbf{I}")
    is the
    ![n \\times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%20n "n \times n")
    identity matrix and
    ![s](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s "s")
    is a scalar. Show that if
    ![\\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda")
    is an eigenvalue of
    ![\\mathbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D "\mathbf{A}")
    with eigenvector
    ![\\mathbf{v}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D "\mathbf{v}"),
    then
    ![\\lambda-s](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda-s "\lambda-s")
    is an eigenvalue of
    ![\\mathbf{C}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BC%7D "\mathbf{C}")
    with the same eigenvector
    ![\\mathbf{v}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D "\mathbf{v}").
    Note:
    ![s](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s "s")
    is often called a **shift**.

4.  Let’s say you had a guess
    ![\\bar{\\lambda}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbar%7B%5Clambda%7D "\bar{\lambda}")
    for an eigenvalue of
    ![A](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A "A")
    and wanted to find the associated eigenvector. Use the previous two
    results to come up with a strategy.

### Problem 3 Solution

1.  

![ \\mathbf{A}\\mathbf{v}=\\lambda \\mathbf{v} \\Rightarrow \\mathbf{A}^{-1}\\mathbf{A}\\mathbf{v} = \\lambda \\mathbf{A}^{-1}\\mathbf{v} \\Rightarrow \\frac{1}{\\lambda}\\mathbf{v}=\\mathbf{A}^{-1}\\mathbf{v}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%20%5Cmathbf%7BA%7D%5Cmathbf%7Bv%7D%3D%5Clambda%20%5Cmathbf%7Bv%7D%20%5CRightarrow%20%5Cmathbf%7BA%7D%5E%7B-1%7D%5Cmathbf%7BA%7D%5Cmathbf%7Bv%7D%20%3D%20%5Clambda%20%5Cmathbf%7BA%7D%5E%7B-1%7D%5Cmathbf%7Bv%7D%20%5CRightarrow%20%5Cfrac%7B1%7D%7B%5Clambda%7D%5Cmathbf%7Bv%7D%3D%5Cmathbf%7BA%7D%5E%7B-1%7D%5Cmathbf%7Bv%7D. " \mathbf{A}\mathbf{v}=\lambda \mathbf{v} \Rightarrow \mathbf{A}^{-1}\mathbf{A}\mathbf{v} = \lambda \mathbf{A}^{-1}\mathbf{v} \Rightarrow \frac{1}{\lambda}\mathbf{v}=\mathbf{A}^{-1}\mathbf{v}.")

2.  You’ll converge to the eigenvector associated with the eigenvalue of
    ![\\mathbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D "\mathbf{A}")
    that is smallest in magnitude.

3.  

![\\mathbf{C}\\mathbf{v}=(A-s\\mathbf{I})\\mathbf{v}=\\mathbf{A}\\mathbf{v}-s\\mathbf{v}=\\lambda \\mathbf{v} - s\\mathbf{v} = (\\lambda-s)\\mathbf{v}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BC%7D%5Cmathbf%7Bv%7D%3D%28A-s%5Cmathbf%7BI%7D%29%5Cmathbf%7Bv%7D%3D%5Cmathbf%7BA%7D%5Cmathbf%7Bv%7D-s%5Cmathbf%7Bv%7D%3D%5Clambda%20%5Cmathbf%7Bv%7D%20-%20s%5Cmathbf%7Bv%7D%20%3D%20%28%5Clambda-s%29%5Cmathbf%7Bv%7D. "\mathbf{C}\mathbf{v}=(A-s\mathbf{I})\mathbf{v}=\mathbf{A}\mathbf{v}-s\mathbf{v}=\lambda \mathbf{v} - s\mathbf{v} = (\lambda-s)\mathbf{v}.")

4.  Let
    ![\\bar{\\lambda}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbar%7B%5Clambda%7D "\bar{\lambda}")
    be your shift
    ![s](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s "s")
    from part (c) above, run the inverse power iteration on
    ![\\mathbf{C}=\\mathbf{A}-\\bar{\\lambda}\\mathbf{I}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BC%7D%3D%5Cmathbf%7BA%7D-%5Cbar%7B%5Clambda%7D%5Cmathbf%7BI%7D "\mathbf{C}=\mathbf{A}-\bar{\lambda}\mathbf{I}"),
    and you’ll find the associated eigenvector. To find the actual
    eigenvalue, you’ll have to add
    ![\\bar{\\lambda}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbar%7B%5Clambda%7D "\bar{\lambda}")
    back to whatever eigenvalue you found in the inverse power
    iteration.

### Problem 4

Assume that
![\\mathbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D "\mathbf{A}")
is a
![5 \\times 5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;5%20%5Ctimes%205 "5 \times 5")
matrix with eigenvalues
![-5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-5 "-5"),
![-2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-2 "-2"),
![1/2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%2F2 "1/2"),
![3/2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;3%2F2 "3/2"),
![4](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;4 "4").

1.  What eigenvalue is expected if you apply power iteration?

2.  What eigenvalue is expected if you apply inverse power iteration?

3.  What eigenvalue is expected if you apply inverse power iteration
    with shift
    ![s=2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s%3D2 "s=2")?

### Problem 4 Solution

1.  Power iternation will converge to
    ![-5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-5 "-5")
    since it is the eigenvalue largest in magntidue.

2.  We consider the reciprocals of the eigenvalues, namely
    ![-1/5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-1%2F5 "-1/5"),
    ![-1/2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-1%2F2 "-1/2"),
    ![2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2 "2"),
    ![2/3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2%2F3 "2/3"),
    ![1/2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%2F2 "1/2").
    The largest is
    ![2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2 "2").
    Therefore, power iteration on
    ![\\mathbf{A}^{-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D%5E%7B-1%7D "\mathbf{A}^{-1}")
    will converge to
    ![2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2 "2"),
    which corresponds to an eigenvalue of
    ![\\mathbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D "\mathbf{A}")
    equal to
    ![1/2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%2F2 "1/2").

3.  First, the shift of
    ![2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2 "2")
    shifts the eigenvalues of
    ![\\mathbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D "\mathbf{A}")
    to
    ![-7](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-7 "-7"),
    ![-4](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-4 "-4"),
    ![-3/2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-3%2F2 "-3/2"),
    ![-1/2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-1%2F2 "-1/2"),
    ![2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2 "2").
    The largest of the reciprocals of these is
    ![-2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-2 "-2"),
    so inverse power iteration on the shifted matrix will converge to
    ![-1/2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-1%2F2 "-1/2").
    After shifting back by 2, we find an eigenvalue of
    ![3/2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;3%2F2 "3/2").
