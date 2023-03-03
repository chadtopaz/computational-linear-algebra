Pset - Solving Linear Systems
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
library(Matrix)
library(knitr)
library(tidyverse)
library(pracma)
library(data.table)
library(R.utils)
library(igraph)
```

### Problem 1

In solving linear systems, we will see algorithms that (hopefully)
produce a sequence of approximations to a true solution. We want to know
if that sequence converges, and if so, how fast.

Let’s review some key concepts about convergence. We say a sequence
${a_n}$ *converges* if it has a finite limit as $n \to \infty$. For
instance, if $a_n = 4n^2/(2n^2+1)$ then the sequence converges to $2$.
On the other hand, if $a_n = (-1)^n + 1$, the limit does not exist so
the sequence diverges.

In addition to asking whether or not a sequence converges, for a
convergent sequence, we can ask about how fast (in some sense) it
converges. As a quick example, think of the sequences ${(1/2)^n}$ and
${(1/2)^{2^n}}$. Both converge to 0 but look at how differently they do
so.

``` r
n <- 0:6
an <- (0.5)^n
bn <- (0.5)^(2^n-1)
kable(cbind(n,an,bn))
```

|   n |       an |        bn |
|----:|---------:|----------:|
|   0 | 1.000000 | 1.0000000 |
|   1 | 0.500000 | 0.5000000 |
|   2 | 0.250000 | 0.1250000 |
|   3 | 0.125000 | 0.0078125 |
|   4 | 0.062500 | 0.0000305 |
|   5 | 0.031250 | 0.0000000 |
|   6 | 0.015625 | 0.0000000 |

The *order of convergence* of a sequence ${a_n}$ is defined as the
number $q \geq 1$ satisfying

$$
\lim_{n \to \infty} \frac{|a_{n+1}-a|}{|a_n - a|^q} = C
$$

such that $C$ is a finite, nonzero number. Here, $a$ is the limit of the
sequence, that is, the value of $a_n$ as $n \to \infty$. Sometimes you
might see this expression rearranged and written more compactly as as

$$
e_{n+1} = C e_n^q
$$

as $n \to \infty$, where $e_n \equiv |a_n - a|$.

To bring it back to computational linear algebra… when solving a problem
with an iterative algorithm, the relevant sequence is the sequence of
errors produced. If we don’t know the true solution, we can’t calculate
the forward error, but we can calculate the backward error. Let’s call
the backward error on the nth iteration $e_n$. We can estimate the
convergence of the algorithm by making a log-log plot of error data.
Taking the log of each side of the equation above, we have

$$
\log e_{n+1} = \log C + q \log e_n.
$$

If we plot $(\log\,e_n,\log\,e_{n+1})$ and observe linear behavior for
large enough $n$, then the slope is the convergence rate $q$.

a\. Calculate (analytically) the order of convergence of the sequence
defined by

$$
a_n = 3 + (1/2)^{2^n - 1}, \quad n \geq 1.
$$

b\. Numerically measure the order of convergence of the sequence from
(a) by making a base-10 log-log plot with nine points, fitting a line to
it, and extracting the slope.

c\. Now let’s investigate numerical convergence of Jacobi’s method for
solving $\mathbf{A}\mathbf{x} = \mathbf{b}$ be a 1000 x 1 column vector
with every entry equal to one. Let $\mathbf{A}$ be a 1000 x 1000 matrix
that I have created for you below.

``` r
n <- 1000
ThreeBanded <- function(n,offset){
  spMatrix(n,n,i=c(1:n,1:(n-1),2:n,(offset+1):n,1:(n-offset)),j=c(1:n,2:n,1:(n-1),1:(n-offset),(offset+1):n),x=c(.5+sqrt(1:n),rep(1,(2*(2*n-1-offset)))))
}
A <- as.matrix(ThreeBanded(n,100))
```

Use Jacobi’s method to solve the system. Run 35 iterations of Jacobi’s
method taking $\mathbf{x}$ as the zero vector for your initial guess and
find the order of convergence of the backwards error in the 2-norm using
the same idea you used in part (b).

### Problem 1 Solution

a\. Your solution goes here.

b\. Your solution goes here.

c\. Your solution goes here.

### Problem 2

Read [this student
project](https://home.csulb.edu/~jchang9/m247/m247_fa11_David_Diego_Alissa_Daniel.pdf)
about linear algebra and ciphers, from [Prof. Jen-Mei
Chang](https://home.csulb.edu/~jchang9/index.html)’s Linear Algebra
course at California State University, Long Beach.

For convenience, I have written a function that converts a string of
text (letters a through z, ignoring case, as well as space, comma, and
period) into a matrix with entries 1 - 29, where the letters are 1 - 26,
the space is 27, and the comma is 28, and the period is 29. We’ll work
with text in groups of N characters, so the matrix produced will have N
rows. There’s no actual ciphering going on yet. We’re just implementing
a way of representing letters (and space and period) straightforwardly
as numbers.

``` r
texttomatrix <- function(text,N){
  key <- c(LETTERS," ",",",".")
  text <- toupper(text)
  newlength <- ceiling(nchar(text)/N)*N
  text <- str_pad(text, newlength, "right")
  tmp <- unlist(strsplit(text, split = ""))
  result <- matrix(match(tmp, key), nrow = N)
  return(result)
}
```

For example, using blocks of 20 letters,

``` r
M <- texttomatrix("Once upon a time there was a very cute puppy.",20)
M
```

    ##       [,1] [,2] [,3]
    ##  [1,]   15   18   21
    ##  [2,]   14    5   16
    ##  [3,]    3   27   16
    ##  [4,]    5   23   25
    ##  [5,]   27    1   29
    ##  [6,]   21   19   27
    ##  [7,]   16   27   27
    ##  [8,]   15    1   27
    ##  [9,]   14   27   27
    ## [10,]   27   22   27
    ## [11,]    1    5   27
    ## [12,]   27   18   27
    ## [13,]   20   25   27
    ## [14,]    9   27   27
    ## [15,]   13    3   27
    ## [16,]    5   21   27
    ## [17,]   27   20   27
    ## [18,]   20    5   27
    ## [19,]    8   27   27
    ## [20,]    5   16   27

I have also written a function that takes this type of matrix and turns
it back into text.

``` r
matrixtotext <- function(M){
  key <- c(LETTERS," ",",",".")
  M <- as.numeric(M)
  result <- trimws(paste(key[M], collapse = ""))
  return(result)
}
```

For example,

``` r
matrixtotext(M)
```

    ## [1] "ONCE UPON A TIME THERE WAS A VERY CUTE PUPPY."

a\. Experiment with each line of code in the `texttomatrix` and
`matrixtotext` functions I wrote for you. Explain what the built-in R
functions `LETTER`, `toupper`, `ceiling`, `nchar`, `str_pad`,
`strsplit`, `unlist`, `match`, `as.numeric` (applied to something that
is already a numerical matrix), `paste` (with `collapse = ""`), and
`trimws` do in their contexts above.

b\. In the example usage of `texttomatrix` above, why are there repeated
values of `27` at the end?

c\. Suppose you are receiving continuous transmissions of ciphertext (as
numbers, stored in a vector) one sentence at a time, and you want to
decipher them real time as they come in. I’ve created a 20 x 20 cipher
matrix for you below.

``` r
set.seed(123)
N <- 20
c <- sample(-300:300, N*N, replace = TRUE)
C <- matrix(c, nrow = N)
```

What is a computationally efficient way to decode the incoming text?
(Hint: in this context, “computationally efficient” means that you
should do as few Gaussian eliminations / matrix inversions as possible.)
Implement your method on the four ciphered sentences stored in the
variables `s1`, `s2`, `s3`, `s4`. These sentences are defined in a code
block below but I have suppressed printing of that code block since
there are very long lists of numbers. Note: before solving, do the
command `options(digits = 10)` because it will help you see when
something unexpected happens that you will need to correct.

### Problem 2 Solution

a\. Your solution goes here.

b\. Your solution goes here.

c\. Your solution goes here.

### Problem 3

Read the introduction to [this Wikipedia article about the infamous
Enron scandal](https://en.wikipedia.org/wiki/Enron_scandal). The network
of email exchanges is available from the [Stanford Large Network Dataset
Collection](http://snap.stanford.edu/data/email-Enron.html). Make sure
you read the brief descrioption of the data set.

The code below creates the network’s adjacency matrix for you.

``` r
url <- "http://snap.stanford.edu/data/email-Enron.txt.gz"
data <- fread(url)
data <- as.matrix(data) + 1
G <- graph_from_edgelist(data)
A <- as_adjacency_matrix(G)
```

One measure of the importance of a node is the Katz centrality. The
vector of Katz centralities **x** satisfies

$$
\mathbf{x} = \bigl((\mathbf{I} - \alpha \mathbf{A}^T)^{-1} - \mathbf{I} \bigr) \mathbf{1}.
$$

where $\mathbf{1}$ is a vector of ones.

You are going to compute the Katz centrality without inverting a matrix.
Some hints:

- Write the system in the form $\mathbf{M}\mathbf{x}=\mathbf{b}$.
- At some point, you will need to create the identity matrix that
  appears in the definition of Katz centrality. Do this using the
  command `Diagonal` (as opposed to the usal `diag`, which is not
  optimal for large matrices).

a\. State what method you should use (remember, no matrix inversion or
Gaussian elimination allowed) and justify mathematically/numerically why
it should work.

b\. For $\alpha = 0.0005$, find the biggest Katz centrality. Solve such
that the backwards error is less than $10^{-10}$ in the infinity norm.
Report the largest Katz centrality, along with the index of the node
that has it. By the way, though the indices of the nodes are not
identified with a persons’s name in the data set, you can at least go
[browse the actual Enron emails](http://www.enron-mail.com/email/) if
you are interested.

### Problem 3 Solution

a\. Your solution goes here.

b\. Your solution goes here.
