Computational Linear Algebra Course Book
================
Chad M. Topaz

- <a href="#r-bootcamp" id="toc-r-bootcamp">R Bootcamp</a>
- <a href="#how-computers-store-numbers"
  id="toc-how-computers-store-numbers">How Computers Store Numbers</a>
- <a href="#fundamentals-of-linear-systems"
  id="toc-fundamentals-of-linear-systems">Fundamentals of Linear
  Systems</a>
- <a href="#solving-linear-systems"
  id="toc-solving-linear-systems">Solving Linear Systems</a>
- <a href="#interpolation" id="toc-interpolation">Interpolation</a>
- <a href="#least-squares-i" id="toc-least-squares-i">Least Squares I</a>
- <a href="#least-squares-ii" id="toc-least-squares-ii">Least Squares
  II</a>
- <a href="#eigenvalues-i" id="toc-eigenvalues-i">Eigenvalues I</a>

# R Bootcamp

## Big Picture

It’s time to review some math and start getting experience with R
programming.

## Goals

- Implement Horner’s method
- Explain advantage of Horner’s method
- Define different ways to multiply vectors
- Compute Taylor polynomials and bound the error

## Horner’s Method

Evaluating a polynomial might seem boring, and admittedly, it sort of
is. But even a simple task like evaluating a polynomial can involve
complexities. Suppose that evaluating a polynomial on your comptuer
takes 0.005 seconds. That seems fast, but what if you are solving a
heavy-duty industrial problem that requires millions of evaluations?

Consider some methods for evaluating the polynomial
$2x^4+3x^3-3x^2+5x-1$:

- Naively,
  $2 \cdot x \cdot x \cdot x \cdot x + 3 \cdot x \cdot x \cdot x \ldots$.

  - 10 multiplications, 4 additions (includes subtractions)

- Store $(1/2)(1/2)$, bootstrap to $(1/2)(1/2)(1/2)$, and so forth.

  - 7 multiplications, 4 additions

- Horner’s method, factoring out successive powers of $x$, as in
  $-1+x(5+x(-3+x(3+2x)))$.

  - 4 multiplications, 4 additions

We can test, using commands from the `tictoc` library.

``` r
xvec <- runif(10^6)
tic()
  for (x in xvec){
    2*x^4 + 3*x^3 - 3*x^2 + 5*x - 1
  }
T1 <- toc()
```

    ## 0.112 sec elapsed

``` r
tic()
for (x in xvec){
    -1 + x*(5 + x*(-3 + x*(3 + 2*x)))
}
T2 <- toc()
```

    ## 0.045 sec elapsed

``` r
t1 <- T1$toc - T1$tic
t2 <- T2$toc - T2$tic
t1/t2
```

    ##  elapsed 
    ## 2.488889

## Inner and Outer Products

There are several different ways to “multiply” vectors $\mathbf{x}$ and
$\mathbf{y}$:

1.  Element-wise. Just multiply each element of $\mathbf{x}$ with the
    element in the corresponding position in $\mathbf{y}$. Note that
    $\mathbf{x}$ and $\mathbf{y}$ must have the same number of elements,
    and the result is a vector with the same number of elements.

2.  Dot product, also known as the inner product. To compute, calculate
    $\mathbf{x} \cdot \mathbf{y} \equiv \mathbf{x}^T \mathbf{y}$ where
    the superscript $T$ indicates transpose. Dot product is related to
    the angle $\theta$ between the two vectors as
    $\mathbf{x} \cdot \mathbf{y} = |\mathbf{x}||\mathbf{y}| \cos \theta$.
    Rearranging this as
    $\mathbf{x} \cdot \mathbf{y}/|\mathbf{y}| = |\mathbf{x}| \cos \theta$
    suggests the intuition of the dot product. It calculates a
    *projection* of one vector onto the other, or restated, it tells us
    how much of one vector is pointing in the direction of the other
    vector. Note that $\mathbf{x}$ and $\mathbf{y}$ must have the same
    number of elements, and the result is a scalar.

3.  Outer product. To compute, calculate
    $\mathbf{x} \otimes \mathbf{y} \equiv \mathbf{x}\,\mathbf{y}^T$.
    Note that $\mathbf{x}$ and $\mathbf{y}$ do not need to have the same
    number of elements, and the result is a matrix. If $\mathbf{x}$ is
    $m \times 1$ and $\mathbf{y}$ is $n \times 1$ then
    $\mathbf{x} \otimes \mathbf{y}$ is $m \times n$.

4.  Cross product. We are not going to worry about this one for now.

Let’s calculate some examples.

``` r
x <- c(1,2,3)
y <- c(4,5,6)
z <- c(7,8,9,10)
x*y         # element-wise
```

    ## [1]  4 10 18

``` r
t(x) %*% y  # dot product
```

    ##      [,1]
    ## [1,]   32

``` r
sum(x*y)    # dot product
```

    ## [1] 32

``` r
x %*% t(z)  # outer product
```

    ##      [,1] [,2] [,3] [,4]
    ## [1,]    7    8    9   10
    ## [2,]   14   16   18   20
    ## [3,]   21   24   27   30

``` r
x %o% z     # outer product
```

    ##      [,1] [,2] [,3] [,4]
    ## [1,]    7    8    9   10
    ## [2,]   14   16   18   20
    ## [3,]   21   24   27   30

``` r
outer(x,z)  # outer product
```

    ##      [,1] [,2] [,3] [,4]
    ## [1,]    7    8    9   10
    ## [2,]   14   16   18   20
    ## [3,]   21   24   27   30

## Taylor’s Theorem

This theorem will be especially useful for error analysis of some
algorithms we use. The basic idea of Taylor’s theorem is that for many
functions, we can approximate the function around a point $x_0$ as a
polynomial of finite degree, called the Taylor polynomial, plus an error
term that accounts for all the higher degree terms we ignored.

More formally, suppose $f$ is $n+1$ times continuously differentiable on
the interval between $x_0$ and $x$. Then

$$
f(x) = \left( \displaystyle \sum_{k=0}^n \frac{f^{(k)}(x_0)}{k!}(x-x_0)^k \right) + \frac{f^{(n+1)}(c)}{(n+1)!}(x-x_0)^{n+1}
$$

where $c$ is an (unknown) number between $x_0$ and $x$.

For example, suppose $f(x)=\sin(x)$ and $x_0=0$. We can pre-compute some
derivatives, $f(0)=\sin(0)=0$, $f^{(1)}(0)=\cos(0)=1$,
$f^{(2)}(0)=-\sin(0)=0$, $f^{(3)}(0)=-\cos(0)=-1$,
$f^{(4)}(0)=\sin(0)=0$, and so on. Then the Taylor polynomial to 4th
order and the error term are:

$$
0 + x + 0x^2 - \frac{x^3}{6} + 0x^4 + \frac{x^5}{120}\cos c.
$$

Here is a plot of the function and the subsequent Taylor approximations.

``` r
P1 <- function(x) {x}
P3 <- function(x) {x-x^3/6}
x <- seq(from = -pi/2, to = pi/2, length = 200)
plot(x, sin(x), col = "black", type = "l")
lines(x, P1(x), col = "green")
lines(x, P3(x), col = "red")
```

![](coursebook_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

We can also calculate a bound on the error if we use the fourth degree
polynomial to approximate $\sin(0.1)$. Note $\sin$ is at most one in
magnitude, so the error term is at most
$(0.1)^5/120\approx 8.333 \times 10^{-8}$. The actual error achieved is

``` r
abs(sin(0.1) - P3(0.1))
```

    ## [1] 8.331349e-08

# How Computers Store Numbers

## Big Picture

Many mathematics and applied mathematics problems can only be solved
numerically, on computers. There is, potentially, error in the way that
computers store numbers and perform arithmetic, and it is crucial to
understand this important source of error.

## Goals

- Convert between decimal and binary
- Explain the form of double precision floating point numbers
- Convert between real numbers and floating point numbers
- Implement nearest neighbor rounding
- Recognize loss of numerical significance
- Implement strategies to reduce loss of significance
- Calculate/estimate the error arising in numerical computations

## Binary Numbers

We usually work in the base 10 system, where each place in the number
represents a power of 10 with a coefficient of 0 through 9, *e.g.*,

$$
13.25 = 1 \times 10^1 + 3\times 10^0 + 2 \times 10^{-1} + 5 \times 10^{-2}.
$$

In the binary system, we instead use powers of 2 with coefficients of 0
through 1, *e.g.*,

$$
1101.01 = 1 \times 2^3 + 1\times 2^2 + 0 \times 2^1 + 1 \times 2^0 + 0 \times 2^{-1} + 1 \times 2^{-2}.
$$

## Machine Numbers

Arithmetic performed on a machine can be different from what we do in
our heads. One fundamental difference is that computers can’t perfectly
store most numbers.

Most computers use a standard representation of real numbers called a
*floating point* number. Sometimes we just say *float* for short. We
represent the floating point version of a number $x$ as
$\mathop{\mathrm{fl}}(x)$, and to reiterate the point above, it is not
necessarily true that $\mathop{\mathrm{fl}}(x)=x$.

We focus on double precision floats. A double precision float consists
of 64 bits (1’s and 0’s) and, expressed in base 10, is equal to
$$(-1)^s 2^{c-1023}(1+f)$$ The parts of this number are:

| Part of number | Name     | Bits              | Notes                                             |
|----------------|----------|-------------------|---------------------------------------------------|
| s              | sign     | 1                 | Determines positive or negative                   |
| c              | exponent | 2 - 12 (11 bits)  | Subtract 1023 to correct exponential bias         |
| f              | mantissa | 13 - 64 (52 bits) | Interpreted as a binary fraction (after radix pt) |

For example, suppose our 64 bits are
$$1 | 10000000111 | 1101000000000000000000000000000000000000000000000000$$
Then expressed in base 10, $s=1$, $c=1024+4+2+1=1031$, and
$f=1/2+1/4+1/16=0.8125$. So the floating point number is
$(-1)^1 \cdot 2^{1031-1023} \cdot 1.8125=-208$.

Now let’s go the other way, namely converting from a decimal number to a
floating point number. In doing so, we’ll see that there can be error in
the computer storage of numbers. Let’s take the number $59\ 1/3$.

1.  Convert to binary. The number is $111011.01010101\ldots$.
2.  Write in normalized form, similar to scientific notation for
    base 10. In normalized form, there is a 1 before the radix point
    (nothing more, nothing less). This yelds
    $1.11011\overline{01}\times 2^5$.
3.  Work backwords to figure out $s$, $c$, and $f$. Since the number is
    positive, the sign is $s=0$. Since we have $2^5$, the exponent must
    equal $5+1023 = 1028$ in base 10, which is $10000000100$. The
    mantissa $f$ is $11011\overline{01}$.

The problem is that $f$ is repeating, but on a computer only holds 52
bits. The computer truncates to 52 bits by nearest neigbor rounding,
which works like this:

- If the 53rd bit is a 0, drop the 53rd bit and everything after it.
- If the 53rd bit is a 1 and there is any nonzero bit after it, add 1 to
  the 52nd bit (and carry as necessary).
- If the 53rd bit is a 1 and there are only zeros after it,
- If the 52nd bit is 1, add 1 to it and carry as necessary
- If the 52nd bit is a 0, do nothing (truncate after 52nd bit)

For our example above, in the mantissa $f$ we have $11011$ followed by
23 copies of the pattern $01$ followed by a $0$. This means the 53rd bit
is a $1$, and there are nonzero bits after it, so we have to add $1$ to
the 52nd. In base 10, we have
$$f = 2^{-1} + 2^{-2} + 2^{-4} + 2^{-5} + 2^{-7} + 2^{-9} + + \ldots + 2^{-51} + 2^{-52} = 0.85416666666666674068.$$
Therefore, on the machine, our number is stored as
$2^5(1.85416666666666674068)=59.333333333333335702 \neq 59\ 1/3$.

The relative error in the storage of this number is
$$\frac{|59\ 1/3 -  59.333333333333335702|}{59\ 1/3} \approx -4 \times 10^{-17}.$$
An important number is *machine epsilon*, denoted $\epsilon_{mach}$,
which is the distance from the number 1 to the next smallest number that
can be represented exactly in floating point form. This distance is
$\epsilon_{mach}=2^{-52}\approx 2.2 \times 10^{-16}$.

## Machine Addition

Machine addition is defined as

$$
\mathop{\mathrm{fl}}(x+y) = \mathop{\mathrm{fl}}(\mathop{\mathrm{fl}}(x)+\mathop{\mathrm{fl}}(y)).
$$

That is, take $x$ and $y$, convert each to a machine number, add them
(exactly, since more registers are available for this operation), and
convert the result to a machine number. Subtraction is just addition
with a negative sign.

For example, let’s add $1$ and $2^{-53}$. Well,
$\mathop{\mathrm{fl}}(1)=1$ and $\mathop{\mathrm{fl}}(2^{-53})=2^{-53}$.
The (exact) sum is $1+2^{-53}$ but due to the rounding rules on
machines, $\mathop{\mathrm{fl}}(1+2^{-53})=1$. Therefore, on a machine,
the sum is $1$. Let’s try:

``` r
1 + 2^(-53)
```

    ## [1] 1

## Loss of Significance

We’ve seen that computer storage of numbers can have error, and
therefore arithmetic can have error. This error is sometimes called
*loss of significance* and the most dangerous operation is subtraction
of nearly equal numbers. Consider the expression $(1-\cos x)/\sin^2 x$,
which can also be written as $1/(1+\cos x)$. We compute this both ways
for $x$ decreasing from $1$.

``` r
options(digits=12)
x <- 10^(-(0:12))
E1 <- (1 - cos(x))/sin(x)^2
E2 <- 1/(1 + cos(x))
kable(cbind(x, E1, E2))
```

<table>
<thead>
<tr>
<th style="text-align:right;">
x
</th>
<th style="text-align:right;">
E1
</th>
<th style="text-align:right;">
E2
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
1e+00
</td>
<td style="text-align:right;">
0.649223205205
</td>
<td style="text-align:right;">
0.649223205205
</td>
</tr>
<tr>
<td style="text-align:right;">
1e-01
</td>
<td style="text-align:right;">
0.501252086289
</td>
<td style="text-align:right;">
0.501252086289
</td>
</tr>
<tr>
<td style="text-align:right;">
1e-02
</td>
<td style="text-align:right;">
0.500012500208
</td>
<td style="text-align:right;">
0.500012500208
</td>
</tr>
<tr>
<td style="text-align:right;">
1e-03
</td>
<td style="text-align:right;">
0.500000124992
</td>
<td style="text-align:right;">
0.500000125000
</td>
</tr>
<tr>
<td style="text-align:right;">
1e-04
</td>
<td style="text-align:right;">
0.499999998628
</td>
<td style="text-align:right;">
0.500000001250
</td>
</tr>
<tr>
<td style="text-align:right;">
1e-05
</td>
<td style="text-align:right;">
0.500000041387
</td>
<td style="text-align:right;">
0.500000000012
</td>
</tr>
<tr>
<td style="text-align:right;">
1e-06
</td>
<td style="text-align:right;">
0.500044450291
</td>
<td style="text-align:right;">
0.500000000000
</td>
</tr>
<tr>
<td style="text-align:right;">
1e-07
</td>
<td style="text-align:right;">
0.499600361081
</td>
<td style="text-align:right;">
0.500000000000
</td>
</tr>
<tr>
<td style="text-align:right;">
1e-08
</td>
<td style="text-align:right;">
0.000000000000
</td>
<td style="text-align:right;">
0.500000000000
</td>
</tr>
<tr>
<td style="text-align:right;">
1e-09
</td>
<td style="text-align:right;">
0.000000000000
</td>
<td style="text-align:right;">
0.500000000000
</td>
</tr>
<tr>
<td style="text-align:right;">
1e-10
</td>
<td style="text-align:right;">
0.000000000000
</td>
<td style="text-align:right;">
0.500000000000
</td>
</tr>
<tr>
<td style="text-align:right;">
1e-11
</td>
<td style="text-align:right;">
0.000000000000
</td>
<td style="text-align:right;">
0.500000000000
</td>
</tr>
<tr>
<td style="text-align:right;">
1e-12
</td>
<td style="text-align:right;">
0.000000000000
</td>
<td style="text-align:right;">
0.500000000000
</td>
</tr>
</tbody>
</table>

The result becomes drastically wrong for small $x$ because of
subtraction of nearly equal numbers. In such problems, be aware and seek
alternative ways to represent the necessary computation, as we did with
the second option above.

# Fundamentals of Linear Systems

## Big Picture

We begin considering the solution of systems of linear equations,
$\mathbf{A} \mathbf{x}=\mathbf{b}$. Linear systems arise in analysis of
different equations (modeling macromolecules, electromagnetics, heat
flow, wave motion, structural engineering, a million other examples), in
curve fitting, in optimization, and many other applications. To
understand solution of linear systems, it’s helpful to recall some
fundamental ideas of linear algebra. Additionally, as we consider the
solution of problems on a computer, we have to think about the effect of
small errors on the solution of the problem. This is called
conditioning. In the linear algebra setting, conditioning is intimately
related to matrix and vector norms.

## Goals

- Interpret linear systems geometrically
- Define linear algebra terms including inverse, determinant,
  eigenvalues, nullspace, linear independence, span, and image
- Give intuitive explanations for the “big theorem” of linear algebra
- Perform the steps of Gaussian elimination
- Establish the computational complexity of Gaussian elimination
- Define forward and backward error
- Define vector norms
- Define matrix norms
- Define condition number for solution of
  $\mathbf{A} \mathbf{x}=\mathbf{b}$

## Linear Systems

For concreteness let’s work in three dimensions and let’s consider the
system

$$ \begin{pmatrix}
2 & 4 & -2\\
1 & 4 & -3\\
-2 & -6 & 7
\end{pmatrix} \begin{pmatrix}
x_1 \\
x_2 \\
x_3
\end{pmatrix} = \begin{pmatrix}
8 \\
8 \\
-3
\end{pmatrix}.$$

Let’s interpret this system two ways.

**Intersecting hyperplanes.** Carrying out the matrix multiplication, we
write the equations

$$
\begin{align*} 2x_1 + 4x_2 -2x_3 & = 8\\
x_1 + 4x_2 -3x_3 & = 8\\
-2x_1 -6x_2 + 7x_3 & = -3 \end{align*}
$$

This form suggests thinking of the set of points that are at the
intersection of these three planes, which could be the empty set, or a
point, or a line, or a plane.

**Vector spans.** Write the equation in terms of its columns, as

$$ x_1 \begin{pmatrix}
2 \\
1 \\
-2
\end{pmatrix} + x_2 \begin{pmatrix}
4 \\
4 \\
-6
\end{pmatrix} + x_3 \begin{pmatrix}
-2 \\
-3 \\
7 \end{pmatrix} = \begin{pmatrix}
8 \\
8 \\
-3
\end{pmatrix} $$

This form suggests thinking of the linear combination of three basis
vectors necessary to reach a particular target vector. There could be 0,
1, or inifinity depending on the arrangement of those vectors.

## Linear Algebra Review and the Big Theorem

Let’s consider an $n \times n$ matrix $\mathbf{A}$.

1\. **Invertible** means that the inverse $\mathbf{A}^{-1}$ exists. This
matrix satisfies $\mathbf{A}\,\mathbf{A}^{-1} = \mathbf{I}$, where
$\mathbf{I}$ is the $n \times n$ identity matrix.

2\. The notation $\det$ means **determinant**. Think of it as a scaling
factor for the transformation defined by a matrix. That is,
multiplication by $\mathbf{A}$ can cause a region to contract
($|\det \mathbf{A}| < 1$) or expand ($|\det \mathbf{A}| > 1$) and/or
reflect ($\det \mathbf{A} < 0$). As an example, let

$$
\mathbf{A} =
\begin{pmatrix}
a & b\\c & d
\end{pmatrix}.
$$

For this example,

$$
\mathbf{A}^{-1} = \frac{1}{ad-bc}
\begin{pmatrix}
d & -b\\
-c & a
\end{pmatrix}.
$$

So $\mathbf{A}$ invertible $\iff \det \mathbf{A} \neq 0$.

3\. The **eigenvalues** $\lambda_i$ of $\mathbf{A}$ satisfy
$\mathbf{A}\,\mathbf{v}_i=\lambda_i \mathbf{v}_i$, where $\mathbf{v}_i$
are the **eigenvectors**. There is a nice property of eigenvalues,
namely $\prod_i \lambda_i = \det \mathbf{A}$, from which it follows that
no $\lambda_i = 0 \iff \det \mathbf{A} \neq 0$.

4\. $\mathbf{A}\,\mathbf{z} \neq 0$ for all
$\mathbf{z} \in \mathbb{R}^n$ except $\mathbf{z}=0 \iff \mathbf{A}$ is
invertible. Why? If $\mathbf{A}\,\mathbf{z} = 0$ for
$\mathbf{z} \neq 0$, then $\mathbf{A}\,\mathbf{z} = 0\,\mathbf{z}$, so 0
is an eigenvalue. But for $\mathbf{A}$ to be invertible, we know $0$
can’t be an eigenvalue. (If you are proof oriented you might notice that
the implication needs to be shown both ways, but I am trying purposely
not to prove here – just to give you some intuition.)

5\. **Null space** just means the parts of $\mathbb{R}^n$ that gets
mapped to $\mathbf{0}$ by $\mathbf{A}$. Another name for nullspace is
**kernel**. Mathematically, the nullspace is all the vectors
$\mathbf{v}$ for which $\mathbf{A}\,\mathbf{v}=0$. So for an invertible
matrix $\mathbf{A}$, the nullspace is $\mathbf{0}$. This is basically
assigning a definition to the previous point (above).

6\. If vectors are **linearly independent**, it means that none of the
vectors can be written as a linear combination of the others. If
$\mathbf{A}$ is invertible, then its columns are linearly independent.
Why? If the columns were linearly dependent, you could take a linear
combination of them to reach $\mathbf{0}$ nontrivially, for instance,

$$
\begin{pmatrix}1 & -2\\2 & -4\end{pmatrix}\begin{pmatrix}x_1 \\ x_2 \end{pmatrix}=\begin{pmatrix}0 \\ 0\end{pmatrix}.
$$

This violates our previous condition (above) about the nullspace only
being $\mathbf{0}$.

7\. **Span** means the set of points reachable by taking linear
combinations of a set of vectors. If you have $n$ linearly independent
vectors in $\mathbb{R}^n$, they span $\mathbb{R}^n$. The **rank** of
$\mathbf{A}$ is just the dimension of the space spanned by the column
vectors.

8\. The **image** of $\mathbf{A}$ means all the points that $\mathbf{A}$
can map to. This is synonymous with the point above: it’s the span of
the columns of $\mathbf{A}$. If the rank of $\mathbf{A}$ is $n$, then
the column vectors are linearly independent so they span $\mathbb{R}^n$,
so the image of $\mathbf{A}$ is $\mathbb{R}^n$.

9\. $\mathbf{A}\,\mathbf{x}=\mathbf{b}$ has a unique solution for all
$\mathbf{b} \iff \mathbf{A}$ is invertible, since you can left multiply
by $\mathbf{A}^{-1}$.

## Gaussian Elimination

First, we define **row echelon form**. A matrix is in row echelon for if

- all nonzero rows (rows with at least one nonzero element) are above
  any rows of all zeroes
- the leading coefficient (the first nonzero number from the left, also
  called the **pivot**) of a nonzero row is always strictly to the right
  of the leading coefficient of the row above it

Gaussian elimination refers to transforming a matrix to row echelon form
by applying the following operations, which do not change the solution
set:

- swap two rows
- multiply a row by a nonzero scalar
- Add one row to a scalar multiple of another

To solve a linear system $\mathbf{A}\,\mathbf{x}=\mathbf{b}$, write it
as an augmented matrix, reduce it to row echelon form, and then use back
substitution to solve.

For example, take

$$
\mathbf{A} = 
\begin{pmatrix}
1 & 3 & 1\\
1 & 1 & -1 \\
3 & 11 & 5
\end{pmatrix}, \quad
\mathbf{b} = \begin{pmatrix}
9 \\ 1 \\ 35
\end{pmatrix}.
$$

Write the augmented matrix

$$
\mathbf{A_a} = \begin{pmatrix}
1 & 3 & 1 & 9 \\
1 & 1 & -1 & 1\\
3 & 11 & 5 & 35
\end{pmatrix}.
$$

Apply $\mathrm{II \leftarrow II - I}$ and
$\mathrm{III \leftarrow III - 3 I}$ to obtain:

$$
\mathbf{A_a} = \begin{pmatrix}
1 & 3 & 1 & 9 \\
0 & -2 & -2 & -8\\
0 & 2 & 2 & 8
\end{pmatrix}.
$$

Apply $\mathrm{III \leftarrow III + II}$ to obtain:

$$
\mathbf{A_a} = \begin{pmatrix}
1 & 3 & 1 & 9 \\
0 & -2 & -2 & -8\\
0 & 0 & 0 & 0
\end{pmatrix}.
$$

The bottom row tells us nothing. The second row tells us there is a free
variable, which we take to be $x_3$. So we solve this equation for
$x_2$, finding $x_2 = 4 -x_3$. This is called back substitution. Then we
do back substitution on the top row, from which we find
$x_1 = 9 - x_3 - 3 x_2 = 9 - x_3 - 3(4 - x_3) = -3 + 2x_3$. Therefore,
the solution is

$$
\mathbf{x} = \begin{pmatrix}  -3 + 2x_3 \\ 4 - x_3 \\ x_3 \end{pmatrix} = \begin{pmatrix}  -3 \\ 4\\ 0 \end{pmatrix} + x_3 \begin{pmatrix}  2 \\ -1 \\ 1 \end{pmatrix}.
$$

By the way, we can go ahead and use a routine I’ve written to perform
the elimination.

``` r
A <- matrix(c(1,3,1,9,1,1,-1,1,3,11,5,35), nrow = 3, byrow = TRUE)
eliminate(A)
```

    ##      [,1] [,2] [,3] [,4]
    ## [1,]    1    3    1    9
    ## [2,]    0   -2   -2   -8
    ## [3,]    0    0    0    0

## Operation Counts and Complexity

Because we are solving problems on computers, we should care about how
long solution takes, which in turn depends on the number of
computational operations carried out. This is called the *complexity* of
the method. For solving $\mathbf{A}\,\mathbf{x}=\mathbf{b}$, with
$\mathbf{A}$ an $n \times n$ matrix, we write the complexity in terms of
$n$. Then, we are usually concerned with the behavior of the operation
count for $n$ large, so we might retain just the leading term in $n$ as
an approximation, or even ignore the coefficient in front that leading
term.

For Gaussian elimination, we have to compute the complexity of the two
stages.

1.  Reduce to echelon form. The number of operations is
    $\frac{2}{3}n^3 + \frac{1}{2}n^2-\frac{7}{6}n = \mathcal{O}(n^3)$.

2.  Back substitute. This takes $n^2 = \mathcal{O}(n^2)$ operations.
    Back substitution is comptutationally cheap compared to row
    reduction. For large enough $n$, the back substitution step is
    negligible since $n^3 \gg n^2$. We can use these operation counts to
    make estimates of how long calculations should take.

For example, suppose row reduction on a $500 \times 500$ matrix takes 1
second. How long does back substitution take? Well, we just use leading
terms in the complexity. Since $n=500$, we have $(2/3)*500^3 t =1$ where
$t$ is time per operation. Solving, $t = 1.2 \times 10^{-8}$. Back
substitution takes $n^2=500^2$ operations, so the total time is
$n^2t=500^2\times 1.2 \times 10^{-8}=0.003$ seconds.

Let’s test scaling of the reduction step for a random matrix.

``` r
set.seed(123)
n1 <- 200
A1 <- matrix(runif(n1^2), ncol = n1)
t1 <- system.time(eliminate(A1))[3]
n2 <- 2*n1
A2 <- matrix(runif(n2^2), ncol = n2)
t2 <- system.time(eliminate(A2))[3]
t2/t1
```

    ##       elapsed 
    ## 5.86956521739

## Forward and Backward Error

In the solution of a computational problem, **forward error** is the
difference between the exact and computed solution, and **backwards
error** is the difference between the original problem and the so-called
modified problem that the approximate solution satisfies. This probably
sounds abstract, so let’s make it concrete in the cases of a
root-finding problem and a linear algebra problem.

Suppose we want to solve $\mathbf{A}\,\mathbf{x}=\mathbf{b}$. The true
solution is $\mathbf{x}$ but our computational method finds an
approximate solution $\mathbf{x}_a$. The forward error is the distance
between the two solutions, that is, $||\mathbf{x}-\mathbf{x_a}||$. The
backward error is the distance between what the matrix outputs when
applied to those solutions, that is,
$||\mathbf{A}\,\mathbf{x}-\mathbf{A}\,\mathbf{x}_a||=||\mathbf{b}-\mathbf{A}\,\mathbf{x}_a||$.
Distance here is the length of the difference between two quantities.

Notice that we haven’t specified what distance means! This is why we
need to define vector and matrix norms.

## Vector Norms

A vector norm is a rule that assigns a real number to every vector.
Intuitively, it measures length. There are a bunch of requirements that
this rule must satisfy in order to be a norm, but rather than stating
those requirements, I’m going to just tell you some practicalities.

The vector norm we’ll work with is called the $p$-norm. The $p$-norm for
$1 \le p \le \infty$ is defined as $$
|| \mathbf{x} ||_p = \left({|x_1|^p + |x_2|^p + \cdots + |x_n|^p} \right)^{1/p}.
$$

The three most common $p$-norms are $p = 1, 2, \infty$ since they are
the easiest to compute with and in some sense are the most natural:

- $p=1$ (the Manhattan or taxicab norm)
  $$|| \mathbf{x} ||_1 =|x_1| + |x_2| + \cdots + |x_n|$$

- $p =2$ (the Euclidean norm)
  $$|| \mathbf{x} ||_2 = \sqrt{x_1^2 + x_2^2 + \cdots + x_n^2\ } = \sqrt{\mathbf{x} \cdot \mathbf{x}}$$

- $p=\infty$
  $$|| \mathbf{x} ||_\infty = \max{\left(| x_1|, | x_2|,  \ldots, |x_n|\right)}$$

For $p = \infty$ it takes a little analysis to show why the
computational definition is what it is, but a numerical study is usually
convincing. We can use the `Norm` command

``` r
v <- c(3,-2,2,3,1,4,1,2,3)
pvals <- c(1,1.5,2,3,4,5,6,7,20)
res <- NULL
for (p in pvals){
  res <- c(res, Norm(v,p))
}
res <- c(res, Norm(v, Inf))
kable(cbind(c(pvals, "Infinity"), res), col.names=c("p", "norm"))
```

<table>
<thead>
<tr>
<th style="text-align:left;">
p
</th>
<th style="text-align:left;">
norm
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
21
</td>
</tr>
<tr>
<td style="text-align:left;">
1.5
</td>
<td style="text-align:left;">
10.5102535215316
</td>
</tr>
<tr>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
7.54983443527075
</td>
</tr>
<tr>
<td style="text-align:left;">
3
</td>
<td style="text-align:left;">
5.55049910291155
</td>
</tr>
<tr>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
4.84053189512475
</td>
</tr>
<tr>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
4.50278575773901
</td>
</tr>
<tr>
<td style="text-align:left;">
6
</td>
<td style="text-align:left;">
4.31746656321528
</td>
</tr>
<tr>
<td style="text-align:left;">
7
</td>
<td style="text-align:left;">
4.20717405025799
</td>
</tr>
<tr>
<td style="text-align:left;">
20
</td>
<td style="text-align:left;">
4.00189474866413
</td>
</tr>
<tr>
<td style="text-align:left;">
Infinity
</td>
<td style="text-align:left;">
4
</td>
</tr>
</tbody>
</table>

To further build intuition, we can plot the unit sphere in
$\mathbb{R}^2$ for various values of $p$, that is, the set of points
that are unit distance from the origin. ![](UnitCircleGrid.png)

## Matrix Norms

The matrix $p$-norm is closely related to the vector $p$-norm, and is
given by

$$
||\mathbf{A}||_p = \max_{\mathbf{x} \not = \mathbf{ 0}} \frac{ || \mathbf{A} \mathbf{x}||_p} { ||\mathbf{x}||_p} = \max_{||\mathbf{x}||_p  = 1}  || \mathbf{A} \mathbf{x}||_p
$$

The matrix $p$-norm says: apply $\mathbf{A}$ to the unit sphere, and
$||\mathbf{A}||_p$ is the length of the vector that is farthest from the
origin. This is not trivial to calculate! You have an infinite number of
vectors to consider.

Fortunately, just like for the vector case, the matrix $p$-norm has a
few special values of $p$ for which it is easy to compute. We have:

- $p=1$

$$
||\mathbf{A} ||_1 = \displaystyle{\max_{1 \le j \le n} \sum_{i=1}^n |a_{ij}|} = \text{maximum absolute column sum}
$$

- $p=2$

$$
|| \mathbf{A} ||_2 = \sqrt{\max\{ \text{eigenvalue}(A^TA) \} }
$$

- $p=\infty$

$$
|| \mathbf{A} ||_\infty = \displaystyle{\max_{1 \le i \le n} \sum_{j=1}^n |a_{ij}|} = \text{maximum absolute row sum}
$$

To see why these definitions are true requires some analysis. If you are
interested, I am happy to point you to proofs.

You can calculate the $1$, $2$, and $\infty$ matrix norms using the R
command `norm`.

``` r
A <- matrix(c(2,-1,1,1,0,1,3,-1,4), byrow = TRUE, nrow = 3)
norm(A, "1")
```

    ## [1] 6

``` r
norm(A, "2")
```

    ## [1] 5.72292695333

``` r
norm(A, "I")
```

    ## [1] 8

There is one really useful identity you should know about matrix norms:

$$||\mathbf{A}\,\mathbf{x}||_p \leq ||\mathbf{A}||_p ||\mathbf{x}||_p.$$

To see this, we start with the right hand side and note

$$
||\mathbf{A}||_p ||\mathbf{x}||_p = \left( \max_{\mathbf{y} \not = \mathbf{ 0}} \frac{ || \mathbf{A} \mathbf{y}||_p} { ||\mathbf{y}||_p} \right) ||\mathbf{x}||_p \geq \frac{ || \mathbf{A} \mathbf{x}||_p} { ||\mathbf{x}||_p} ||\mathbf{x}||_p = || \mathbf{A} \mathbf{x}||_p.
$$

## Condition Number for $\mathbf{A}\,\mathbf{x} = \mathbf{b}$

Let us consider solving $\mathbf{A}\,\mathbf{x}=\mathbf{b}$. Suppose we
find approximate solution $\mathbf{x}_a$. The **relative forward error**
is

$$
\frac{||\mathbf{x}-\mathbf{x}_a||}{||\mathbf{x}||}
$$

and the **relative backward error** is

$$
\frac{||\mathbf{A}\mathbf{x}-\mathbf{A}\mathbf{x}_a||}{||\mathbf{A}\mathbf{x}||}=\frac{||\mathbf{b}-\mathbf{A}\mathbf{x}_a||}{||\mathbf{b}||}.
$$

We define **error magnification** as the ratio

$$
\frac{\text{relative forward error}}{\text{relative backward error}}=\frac{\frac{||\mathbf{x}-\mathbf{x}_a||}{||\mathbf{x}||}}{\frac{||\mathbf{b}-\mathbf{A}\mathbf{x}_a||}{||\mathbf{b}||}}.
$$

The **condition number** $\kappa(\mathbf{A})$ is the largest possible
error magnification (over all possible $\mathbf{x}$). Or restated, it’s
the worst possible ratio of relative forward error to relative backward
error.

Why do we care about this? It tells us whether we expect a small
residual to imply a small error in the solution or not. Let’s make this
concrete with an example.

Consider:

- Let
  $\mathbf{A} = \begin{pmatrix}0.913 & 0.659 \\ 0.457 & 0.330 \end{pmatrix}$
- Then $\kappa_2(A) = 1.25\times10^4$
- Let $\mathbf{b} = \begin{pmatrix} 0.254 \\ 0.127 \end{pmatrix}$
- Then $\mathbf{x} =(1,-1)$.
- Consider two approximate solutions $\mathbf{x}_{1,2}$

$$
\begin{array}{lcl}
\mathbf{x}_1 = (-0.0827,0.5) && \mathbf{x}_2 = (0.999,-1.001) \\
\triangle \mathbf{x}_1 = (1.0827, -1.5) && \triangle \mathbf{x}_2 = (0.001,0.001) \\
||\triangle \mathbf{x}_1 || = 1.85&& ||\triangle \mathbf{x}_2|| = .0014\\
||\triangle \mathbf{x}_1 ||/||\mathbf{x}|| = 1.308&& ||\triangle \mathbf{x}_2||/||\mathbf{x}|| = .001\\
\mathbf{b}_1 = (0.2539949, 0.1272061) &&  \mathbf{b}_2 = (0.252428, 0.126213) \\
\triangle\mathbf{b}_1 = (0.0000051,- 0.0002061) &&  \triangle\mathbf{b}_2 = (0.001572, 0.000787) \\
||\triangle \mathbf{b}_1 || = 0.000206&& ||\triangle \mathbf{b}_2|| = .00176\\
||\triangle \mathbf{b}_1 ||/||\mathbf{b}|| = 0.000726&& ||\triangle \mathbf{b}_2||/||\mathbf{b}|| = .0062\\
mag = 1.8 \times 10^3 && mag = 1.6 \times 10^1
\end{array}
$$

We can go ahead and calculate the actual condition number of the matrix.
R has a command called `kappa` that computes the condition number
approximately by default, or exactly if specified.

``` r
A <- matrix(c(0.913,0.659,0.457,0.330), nrow = 2, byrow = TRUE)
kappa(A)
```

    ## [1] 14132.0316376

``` r
kappa(A, exact = TRUE)
```

    ## [1] 12485.031416

## Calculating the Condition Number

Remember that the condition number isn’t merely an error magnification –
it’s the maximum possible error magnification. Computing $\kappa$
exactly using this definition is impossible because there are an
infinite number vectors one must consider $\mathbf{A}$ acting on.

Fortunately, there’s another way to calculate condition number:

$$
\kappa_p(\mathbf{A})=||\mathbf{A}||_p||\mathbf{A}^{-1}||_p
$$

The derivation of this identity is a few lines of linear algebra that I
am happy to show you if you are interested. We can check it numerically
for now.

``` r
set.seed(123)
N <- 10
A <- matrix(runif(N^2), nrow = N)
norm(A, "2") * norm(solve(A), "2")
```

    ## [1] 47.2939732454

``` r
kappa(A, norm = "2", exact = TRUE)
```

    ## [1] 47.2939732454

# Solving Linear Systems

## Big Picture

We will learn two methods of solving systems of linear equations. First,
depending on the context, solution by Gaussian elimination can be
computationally costly. Sometimes it is better to decompose (factor) the
matrix. There are a number of useful ways to do this but the one we will
focus on is called **LU** decomposition. Second, every method you have
seen thus far for solving linear systems provides exact solutions
(excluding numerical error). These are called direct methods. While we
like exact solutions, a disadvantage is that these methods involve at
least one step that is $\mathcal{O}(n^3)$. If this computational cost is
prohibitive, we can consider a potentially faster iterative method at
the expense of giving up having an exact solution. The one we will learn
is called Jacobi’s method.

## Goals

- Define and implement **LU** decomposition and state the potential
  advantages
- Explain how fixed point iteration relates to solution of an equation
- Derive and implement Jacobi’s method, and explain advantages and
  limitations
- Use convergence criteria for Jacobi’s method
- State other iterative solution methods

## LU Decomposition

Recall that when solving an $n \times n$ system
$\mathbf{A}\,\mathbf{x} = \mathbf{b}$ with Gaussian elimination, the
elimination step is $\mathcal{O}(n^3)$ and back substitution is
$\mathcal{O}(n^2)$. In some applications, it is necessary to solve

$$
\mathbf{A}\,\mathbf{x} = \mathbf{b}_1, \quad \mathbf{A}\,\mathbf{x} = \mathbf{b}_2, \quad \mathbf{A}\,\mathbf{x} = \mathbf{b}_3, \quad \ldots, \quad \mathbf{A}\,\mathbf{x} = \mathbf{b}_M
$$

where $\mathbf{A}$ is the same each time and $M$ is large. $\mathbf{A}$
itself needs the same row reductions each time. Only the augmented part
$\mathbf{b}$ changes. It would be a waste of computation to run Gaussian
elimination $M$ times.

**LU** decomposition is a way of storing the Gaussian elimination steps
in matrix form so that they can be applied to many $\mathbf{b}$. We take
$\mathbf{A}$ and decompose (or factorize) it as the product

$$
A  = 
\underbrace{\begin{array}{|cccccc|}
\hline 1 &&&&&\\
\ast & 1 &&&&\\
\ast & \ast  & 1 &&&\\
\ast & \ast & \ast & 1 &&\\
\ast & \ast & \ast & \ast & 1  &\\
\ast & \ast & \ast & \ast & \ast & 1 \\ \hline
\end{array}}_{\mathbf{L}=\text{Lower unit triangular}} \
\underbrace{\begin{array}{|cccccc|}
\hline \ast & \ast & \ast & \ast & \ast & \ast \\ 
&\ast & \ast & \ast & \ast & \ast  \\
&&\ast & \ast & \ast & \ast \\
&&&\ast & \ast  & \ast \\
&&&&\ast & \ast \\
&&&&&\ast \\ \hline
\end{array}}_{\mathbf{U}=\text{Upper triangular}}
$$

The $\mathbf{L}$ matrix encodes the multipliers used to eliminate
elements during Gaussian elimination and the $\mathbf{U}$ matrix is the
result of the elimination process. Therefore, putting $\mathbf{A}$ into
its **LU** factorization takes one application of Gaussian elimination,
or approximately $\frac{2}{3} n^3$ operations. Solving **LUx** = **b**
requires two back substitutions, namely one to solve
$\mathbf{L}\,\mathbf{y} = \mathbf{b}$ for $\mathbf{y}$ and one to solve
$\mathbf{U}\,\mathbf{x}=\mathbf{y}$ for $\mathbf{x}$. This takes $2n^2$
operations. So, to solve
$\mathbf{A}\,\mathbf{x} = \mathbf{b}_1, \ldots, \mathbf{A}\,\mathbf{x} = \mathbf{b}_M$
takes approximately $\frac{2}{3}n^3 + 2 M n^2$ operations, in contrast
to $\frac{2}{3}Mn^3 + Mn^2$ for Gaussian elimination.

The **LU** decomposition exists if and only if the upper-left sub-blocks
$\mathbf{A}_{1:k,1:k}$ are non-singular for all $1\leq k \leq n$ (not
proven here). If the decomposition exists, it is unique.

Let’s do an example of how **LU** decomposition works. Take

$$
\mathbf{A} = 
\begin{pmatrix}
1 & 3 & 1\\
1 & 1 & -1 \\
3 & 11 & 5
\end{pmatrix}.
$$

Start by defining $\mathbf{U} = \mathbf{A}$ (it’s not upper triangluar
yet, but I am still calling it $\mathbf{U}$) and

$$
\mathbf{L} = \mathbf{I}_3 = 
\begin{pmatrix}
1 & 0 & 0\\
0 & 1 & 0 \\
0 & 0 & 1
\end{pmatrix}.
$$

Apply $\mathrm{II} \leftarrow \mathrm{II} - 1\cdot \mathrm{I}$ and
$\mathrm{III} \leftarrow \mathrm{III} - 3\cdot \mathrm{I}$ to
$\mathbf{U}$, so

$$
\mathbf{U} = \begin{pmatrix}
1 & 3 & 1  \\
0 & -2 & -2\\
0 & 2 & 2
\end{pmatrix}.
$$

Also, store the multupliers $1$ and $3$ in their corresponding rows in
the first column of $\mathbf{L}$, so that

$$
\mathbf{L} =
\begin{pmatrix}
1 & 0 & 0\\
1 & 1 & 0 \\
3 & 0 & 1
\end{pmatrix}.
$$

Now apply $\mathrm{III} \leftarrow \mathrm{III} + \mathrm{II}$, which
yields

$$
\mathbf{U} = \begin{pmatrix}
1 & 3 & 1  \\
0 & -2 & -2\\
0 & 0 & 0
\end{pmatrix}
$$

and

$$
\mathbf{L} =
\begin{pmatrix}
1 & 0 & 0\\
1 & 1 & 0 \\
3 & -1 & 1
\end{pmatrix}.
$$

Since $\mathbf{U}$ is in echelon form, we are done! We can check that
our decomposition worked.

``` r
A <- matrix(c(1,3,1,1,1,-1,3,11,5),nrow=3,byrow=TRUE)
L <- matrix(c(1,0,0,1,1,0,3,-1,1),nrow=3,byrow=TRUE)
U <- matrix(c(1,3,1,0,-2,-2,0,0,0),nrow=3,byrow=TRUE)
A - L %*% U
```

    ##      [,1] [,2] [,3]
    ## [1,]    0    0    0
    ## [2,]    0    0    0
    ## [3,]    0    0    0

We can also use the `lu` command.

``` r
sol <-  lu(A, scheme = "ijk")
L <- sol$L
U <- sol$U
A - L%*%U
```

    ##      [,1] [,2] [,3]
    ## [1,]    0    0    0
    ## [2,]    0    0    0
    ## [3,]    0    0    0

Let’s use Gaussian elimination with the `echelon` command and **LU**
decomposition with the to compare the time for solving a
$100 \times 100$ system for $100$ different right hand sides.

``` r
n <- 10
A <- matrix(runif(n^2),nrow=n)
set.seed(123)
t1 <- system.time(
  for (i in 1:n){
    b <- runif(n)
    x <- echelon(A,b)
  }
)[3]
t2 <- system.time(
  LU <- lu(A, scheme = "ijk")
)[3]
L <- LU$L
U <- LU$U
set.seed(123)
t3 <- system.time(
  for (i in 1:n){
    b <- runif(n)
    y <- forwardsolve(L,b)
    x <- backsolve(U,y)
  }
)[3]
as.numeric(t1/(t2 + t3))
```

    ## [1] 6.66666666667

## Fixed Point Iteration

Now we turn to iterative methods. Sometimes you can solve a problem by a
method called **fixed point iteration** whereby you just keep plugging
into an expression until the output equals the input. For example,
suppose you want to solve $(x-3)(x+1)=x^2-2x-3=0$. Pretend you don’t
know where the roots are but you think there is one near x = -2, so you
start out with that guess. You also notice you can write
$x^2 - 2x - 3 =0$ as $x = 3/(x-2)$. So you define an iteration
$x_{i+1} = 3/(x_i-2)$.

``` r
x <- -2
for (i in 1:10){
  x <-  3/(x-2)
}
```

The iteration seems to be converging to a root. What if we try a guess
near the other root?

``` r
x <- 3.00001
for (i in 1:20){
  x <-  3/(x-2)
}
```

The iteration does not have to converge to the solution we want, and in
fact, it doesn’t have to converge at all.

But if it converges, it’s pretty nifty. It’s computationally cheap – all
we have to do is evaluate the right hand side of our iteration
repeatedly.

## Jacobi Iteration

Let’s take this idea and apply it to solving
$\mathbf{A}\,\mathbf{x}=\mathbf{b}$. Let
$\mathbf{A} = \mathbf{D} + \mathbf{R}$ where $\mathbf{D}$ contains the
diagonal elements of $\mathbf{A}$ and $\mathbf{R}$ contains everything
else. Then we can write

$$
\begin{align}
\mathbf{A}\,\mathbf{x} &= \mathbf{b} \\
(\mathbf{D}+\mathbf{R})\mathbf{x} &= \mathbf{b} \\
\mathbf{D}\mathbf{x} + \mathbf{R} \mathbf{x} &= \mathbf{b}\\
\mathbf{D} \mathbf{x} &= \mathbf{b} - \mathbf{R} \mathbf{x}\\
\mathbf{x} &= \mathbf{D}^{-1} (\mathbf{b}-\mathbf{R}\mathbf{x})
\end{align}
$$

We can consider this a fixed point iteration,

$$
\mathbf{x}_{i+1} = \mathbf{D}^{-1} (\mathbf{b}-\mathbf{R}\,\mathbf{x}_i).
$$

Here are some computational advantages of this method.

1.  The matrix $\mathbf{D}$ is very cheap to invert, because it is a
    diagonal matrix. The inverse is simply the diagonal matrix with the
    reciprocals of the original elements.
2.  Each iteration is only $\mathcal{O}(n^2)$, and if the matrix is
    sparse, it is even cheaper.

The usual way we stop iterating is that we choose in advance a threshold
for the relative backwards error and stop when we fall below it.

Here’s an example.

``` r
set.seed(123)
n <- 1000
A <- matrix(runif(n^2), nrow = n)
diag(A) <- n + 1 + diag(A)
b <- runif(n)
t1 <- system.time(
  xexact <- solve(A,b)
)[3]
xapprox <- rep(0,n)
d <- diag(A)
R <- A - diag(d)
err <-  Inf
tol <- 1e-10
t2 <- system.time(
  while (err > tol){
    xapprox <- (b - R%*%xapprox)/d
    err <- Norm(b-A%*%xapprox,Inf)/Norm(b,Inf)
  }
)[3]
Norm(xapprox-xexact,Inf)
```

    ## [1] 3.64232569886e-14

``` r
t1/t2
```

    ##       elapsed 
    ## 2.47297297297

## Convergence of Jacobi’s Method

There’s no reason at all to expect that Jacobi’s method converges.
There’s a really useful theorem that says it it will converge if and
only if the eigenvalues of $-\mathbf{D}^{-1}\mathbf{R}$ are all less
than one in magnitude. If this criterion is met, then the closer to one
the eigenvalues are in magnitude, the slower convergence will be. This
takes about half a page to prove and is a worthwhile exercise to
understand, so feel free to ask me for the proof.

However, calculating the eigenvalues of that matrix could be really hard
and costly! There’s a sufficient condition for convergence that is much
easier to check computationally, namely that $\mathbf{A}$ is **strictly
diagonally dominant**. This means

$$
|a_{i,i}| > \sum_{j \not=i} |a_{i,j}| \quad\text{in each row $i$.}
$$

## Other Iterative Methods

There are other iterative solution methods for linear systems that all
are inspired by Jacobi’s method. Some of these include Gauss-Seidel
iteration and Successive Over-Relaxation. DYou are welcome to ask me
about these.

# Interpolation

## Big Picture

Now we enter into the part of this course that is about data. As
scientists, often we will have access only to noisy or partial data such
as a sound signal, a visual image, geograhpic data, network data, and so
forth. We may wish to fill in the missing data, make predictions, smooth
the data, compress it, differentiate or integrate it, visualize it, and
more. These are data processing tasks. The first approach we present is
polynomial interpolation and approximation. Interpolation refers to
constructing new data points within the range of a discrete set of known
points. Approximation refers to finding an approximation to given
function by choosing a function from a predetermined class, which in
this case, is polynomials.

## Goals

- Explain the advantage of using polynomials to describe data
- Implement Vandermonde interpolation and explain its pros/cons
- Implement Lagrange inteprolation and explain its pros/cons
- Explain and implement data compression via interpolation
- Describe and recognize Runge’s phenomenon
- State the error term for polynomial interpolation and bound it
- Explain the advantages of Chebyshev integration and implement the
  technique
- Compare approaches to interpolation
- Explain advantages of interpolating data with cubic splines
- Implement cubic spline interpolation

## Why Polynomials?

Suppose we have incomplete data and we’d like to estimate a piece of
information that we don’t have.

``` r
# I created some mystery data and hid it from you
plot(x,y)
```

![](coursebook_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

We might try to use polynomials to describe the data, and then glean
information from the polynomial. Polynomials are convenient for several
reasons.

First, for any function that is defined and continuous on a closed,
bounded interval, there exists a polynomial that is as “close” to the
given function as desired. Stated more precisely, this is the
**Weierstrass Approximation Theorem**. Suppose $f$ is defined and
continuous on $[a,b]$. For each $\epsilon > 0$, there exists a
polynomial $P(x)$ with the property that $|f(x)-P(x)| < \epsilon$ for
all $x \in [a,b]$.

Second, polynomials are easy to handle computationally. It is fast to
evaluate them using Horner’s method and it is straightforward to
integrate and differentiate them.

Though polynomials are convenient, we won’t use Taylor polynomials.
Taylor polynomials can do a great job accurately describing a function
near a single point, but we are now instead interested in getting a good
approximation over an interval.

In contrast, a polynomial that passes through every data point is called
the interpolating polynomial.

Even though I haven’t yet told you how to find an interpolating
polynomial, let’s compare the Taylor and interpolating approaches.
Consider $f(x) = 1/x$ on the interval $[1,3]$. The second degree Taylor
polynomal through $x=1$ (as an example) is $T(x)=3-3x+x^2$. The
interpolating polynomial using three points sampled equally across the
interval, namely $(1,1),(2,1/2),(3,1/3)$, is $11/6 - x + (1/6)x^2$.

``` r
f <- function(x){1/x}
tee <- function(x){3 - 3*x + x^2}
p <- function(x){11/6 - x + 1/6*x^2}
xdata <- seq(from = 1, to = 3, length = 3)
x <- seq(from = 1, to = 3,length = 200)
plot(x, f(x), type = "l", lwd = 2)
lines(x, tee(x), col = "red", lwd = 2)
points(xdata, f(xdata), cex = 2)
lines(x, p(x), col = "green", lwd = 2)
```

![](coursebook_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

A few important points to know about the interpolating polynomial
through $n$ points with distinct $x$ coordinates include:

- It exists
- It is unique
- It has degree at most $n-1$

Let’s examine some different ways to compute the interpolating
polynomial.

## Vandermonde Matrix

I also call this method of interpolation “brute force.” Let’s start with
an example. Suppose we have three data points $(x_1,y_1)$, $(x_2,y_2)$,
$(x_3,y_3)$. The interpolating polynomial is
$P(x) = c_0 + c_1 x + c_2 x^2$, where we need to determine the
coefficients. Let’s determine them simply by plugging in. We have the
equations

$$
\begin{align*}
c_0 + c_1 x_1 + c_2 x_1^2 &= y_1 \\
c_0 + c_1 x_2 + c_2 x_2^2 &= y_2 \\
c_0 + c_1 x_3 + c_2 x_3^2 &= y_3
\end{align*}.
$$

We can write this in matrix form as

$$
\begin{pmatrix}
1 & x_1 & x_1^2 \\
1 & x_2 & x_2^2 \\
1 & x_3 & x_3^2
\end{pmatrix} \begin{pmatrix}
c_0 \\ c_1 \\ c_2
\end{pmatrix} =
\begin{pmatrix}
y_1 \\ y_2 \\ y_3
\end{pmatrix}.
$$

and solve the linear system to find the coefficients $c_i$. More
generally, for $n$ data points, the problem is

$$
\begin{pmatrix}
1 & x_1 & x_1^2 &  & x_1^{n-1} \\ 
1 & x_2 & x_2^2 & \cdots & x_2^{n-1} \\ 
1 & x_3 & x_3^2 &  & x_3^{n-1} \\ 
&  \vdots &  &  & \vdots \\ 
1 & x_n & x_n^2 & \cdots & x_n^{n-1} \\ 
\end{pmatrix} \begin{pmatrix} c_0 \\ c_1 \\ c_2 \\ \vdots \\ c_{n-1} \end{pmatrix} = 
\begin{pmatrix} y_1 \\ y_2 \\ y_3 \\ \vdots \\ y_n \end{pmatrix}.
$$

It can be proven that if the $x_i$ are distinct, the matrix has nonzero
determinant, and hence the system is solvable with a unique solution.

In solving problems using the Vandermonde matrix, we are using a basis
for the interpolating polynomial that is $\{1,x,x^2,\ldots\}$. This
seems very natural since this is how we usually think of polynomials!
The problems is that to find the coefficients in this basis, we have to
solve a linear problem whose matrix is very ill-conditioned. Let’s see
what happens if we sample more and more points from our function and
construct the interpolating polynomial. We’ll look at $\kappa$ for the
Vandermonde matrix.

``` r
nvals <- 2^(1:8)
kappavals <- NULL
for (n in nvals){
  x <- seq(from = 1, to = 3, length = n)
  kappavals <- c(kappavals, kappa(vander(x)))
}
kable(cbind(nvals, kappavals))
```

<table>
<thead>
<tr>
<th style="text-align:right;">
nvals
</th>
<th style="text-align:right;">
kappavals
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
6.25000000000e+00
</td>
</tr>
<tr>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1.67854715706e+03
</td>
</tr>
<tr>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
3.21899099387e+08
</td>
</tr>
<tr>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
2.08214238608e+19
</td>
</tr>
<tr>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
2.00880822269e+29
</td>
</tr>
<tr>
<td style="text-align:right;">
64
</td>
<td style="text-align:right;">
5.35079511091e+45
</td>
</tr>
<tr>
<td style="text-align:right;">
128
</td>
<td style="text-align:right;">
2.74099311576e+76
</td>
</tr>
<tr>
<td style="text-align:right;">
256
</td>
<td style="text-align:right;">
2.55777379858e+137
</td>
</tr>
</tbody>
</table>

## Lagrange Interpolating Polynomial

Let’s seek an alternative method that gets around these issues. There’s
actually a way to simply write down the interpolating polynomial without
doing any calculation of coefficients at all. Let’s do an example.

For points $(1,10), (2, 6), (3,4), (4,10)$, consider the polynomial

$$\begin{aligned} p(x) &= 10 \frac{(x-2)(x-3)(x-4)}{(1-2)(1-3)(1-4)} + 6 \frac{(x-1)(x-3)(x-4)}{(2-1)(2-3)(2-4)} \\ &+ 4 \frac{(x-1)(x-2)(x-4)}{(3-1)(3-2)(3-4)} + 10 \frac{(x-1)(x-2)(x-3)}{(4-1)(4-2)(4-3)} \end{aligned}$$

First, notice the pattern of how each term is constructed: there’s a
coefficient times a polynomial. The coefficient is the $y$-value of an
interpolation point. The polynomial is one that is equal to $1$ at that
the $x$-coordinate of that point, and $0$ at the $x$-coordinate of every
other point. So now we ask…

- Is the expression a polynomial?
- Is it degree at most $n-1$?
- Does it pass through each data point?

Since the answer to all these questions is yes, it is the interpolating
polynomial.

If we expand out our Lagrange polynomial we find
$p(x) = x^3 - 5 x^2 + 4 x+10$. We can use `R` to verify this result by
solving the Vandermonde problem.

``` r
x <- c(1,2,3,4)
y <- c(10,6,4,10)
c <- solve(vander(x), y)
c
```

    ## [1]  1 -5  4 10

Following the pattern we established above, the Lagrange polynomial for
points $(x_1,y_1),\ldots,(x_n,y_n)$ is

$$
p(x) = \sum_{i = 1}^n y_i \prod_{j \not = i} \frac{(x - x_j)}{(x_i - x_j)}.
$$

The advantage of this method is that it doesn’t require any numerical
solution… just evaluation. Let’s try a comparison: Vandermonde
vs. Lagrange. Here, I’ll use R’s `lagrangeInterp` function.

``` r
set.seed(123)
n <- 10
x <- 1:n
y <- runif(n)
x0 <- 1.5
c <- solve(vander(x),y)
horner(c, x0)$y
```

    ## [1] -0.157060132098

``` r
lagrangeInterp(x,y,x0)
```

    ## [1] -0.157060132133

We can also do a speed comparison test.

``` r
set.seed(123)
numTrials <- 1000
n <- 10
x <- 1:n
tic()
for (i in 1:numTrials){
  y <- runif(n)
  c <- echelon(vander(x), y)[, n+1]
  horner(c, x0)
}
T1 <- toc()
```

    ## 1.984 sec elapsed

``` r
t1 <- T1$toc - T1$tic
tic()
for (i in 1:numTrials){
  y <- runif(n)
  lagrangeInterp(x, y, x0)
}
T2 <- toc()
```

    ## 0.029 sec elapsed

``` r
t2 <- T2$toc - T2$tic
as.numeric(t2/t1)
```

    ## [1] 0.0146169354839

## Data Compression

One of the powerful things interpolation can do is compress data. Let’s
do an example. Suppose that we need to know values for the function
$\sin x$. A computer doesn’t magically know this function, so it has to
have some way to compute/evaluate it. One option would be to store a
giant look-up table. There are an infinite number of numbers to store,
though, even for the interval $[0,2\pi)$. Another option is to fit a
polynomial based on a finite number of points, store the coefficients,
and evaluate the polynomial as needed. Let’s do this using 5 points to
begin with. We’ll write a function that takes a specified number of
points, samples them from the function, constructs the interpolating
polynomial, plots the function and the polynomial, and calculates the
maximum error. The inputs are your $x$ data, your $y$ data, and the $x$
values at which you’d like interpolated values.

``` r
interperror <- function(n, plotflag = FALSE){
  xdata <- seq(from = 0, to = 2*pi, length = n)
  ydata<- sin(xdata)
  x <- seq(from = 0, to = 2*pi, length = 1000)
  y <- sin(x)
  yinterp <- lagrangeInterp(xdata, ydata, x)
  if (plotflag == TRUE){
    plot(x, y, type = "l")
    points(xdata, ydata)
    lines(x, yinterp, col = "blue")
  }
  error <- max(abs(y - yinterp))
  return(error)
}
interperror(5, plotflag = TRUE)
```

![](coursebook_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

    ## [1] 0.180757796555

Not bad for just 5 points. Let’s examine how the error changes as a
function of $n$.

``` r
nvec <- 2:20
errorvec <- NULL
for (n in nvec){
  errorvec <- c(errorvec, interperror(n))
}
```

    ## Warning in max(abs(y - yinterp)): no non-missing arguments to max; returning
    ## -Inf

``` r
orderofmag <- round(log10(errorvec))
```

    ## Warning: NaNs produced

``` r
plot(nvec, orderofmag)
```

![](coursebook_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

This means that we can represent the sine function with $10^{-13}$ error
using only 20 pieces of information, instead of storing a huge lookup
table.

## Runge’s Phenomenon

We have been thinking about data and talking about polynomial
interpolation as a way of estimating it, compressing it, and
representing it conveniently to do other mathematical operations on it.
So far, we have only examined error numerically. Now it is time to look
at the error in more detail, including finding out when it is
potentially large and thinking about how we can reduce it.

Before thinking about error in polynomial interpolation, let’s start out
with an example that is about Taylor polynomials. We’ll construct Taylor
polynomials of increasing degree to estimate the function
$f(x) = \cos(x)$ around the point $x_0=0$ on $[0,2\pi]$. The $n$th
degree taylor polynomial is

$$
\sum_{i=1}^n (-1)^{i/2} \frac{x^i}{i!}.
$$

``` r
mytaylor1 <- function(x,n){
  ans <- 0
  for (i in seq(from = 0, to = n, by = 2)){
    ans <- ans + (-1)^(i/2)*x^i/factorial(i)
  }
  return(ans)
}
x <- seq(from = 0, to = 2*pi, length = 1000)
plot(x, cos(x), type = "l", col = "red", lwd = 5, xlim = c(0,2*pi), ylim = c(-1.5,1.5), xlab = "x", ylab = "y")
for (n in seq(from = 0,to = 14,by = 2)){
  lines(x, mytaylor1(x, n))
}
```

![](coursebook_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

So more terms are better, right? Let’s try again with the function
$f(x) = 1/x$ around the point $x_0=1$. The $n$th degree Taylor
polynomial is

$$
\sum_{i=1}^n (-1)^i (x-1)^i.
$$

``` r
mytaylor2 <- function(x,n){
  ans <- 0
  for (i in seq(from = 0, to = n, by = 1)){
    ans <- ans + (-1)^(i)*(x-1)^(i)
  }
  return(ans)
}
x <- seq(from = 1, to = 2.5, length = 1000)
plot(x,1/x, type = "l", col = "red", lwd = 5, xlim = c(1,2.5), ylim = c(0,1.5), xlab = "x", ylab = "y")
for (n in seq(from = 0, to = 40, by = 4)){
  lines(x, mytaylor2(x, n))
}
```

![](coursebook_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

Oh! I guess that more isn’t always better.

Now let’s think now about interpolating polynomials. Let’s consider
$\cos(x)$ with $n$ equally sampled points across $[0,2\pi]$ for
different values of $n$.

``` r
x <- seq(from = 0, to = 2*pi, length = 1000)
y <- cos(x)
plot(x, y, type = "l", col = "red", lwd = 5, xlim = c(0,2*pi), ylim = c(-1.1,1.1))
nvec <- 3:20
error <- NULL
for (n in nvec){
  xdata <- seq(from = 0, to = 2*pi, length = n)
  ydata <- cos(xdata)
  yinterp <- lagrangeInterp(xdata, ydata, x)
  lines(x, yinterp)
  error <- c(error, max(abs(y - yinterp)))
}
```

![](coursebook_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

``` r
plot(nvec, log10(error), xlab = "n", ylab = "log10 of error")
```

![](coursebook_files/figure-gfm/unnamed-chunk-32-2.png)<!-- -->

Looks good! Let’s try again with a different function,
$f(x) = (1+x^2)^{-1}$ on $[-1,1]$.

``` r
x <- seq(from = -5, to = 5, length = 1000)
y <- 1/(1+x^2)
plot(x, y, type = "l", col = "red", lwd = 5, xlim = c(-5,5), ylim = c(-3,3))
nvec <- seq(from = 3, to = 33, by = 6)
error <- NULL
for (n in nvec){
  xdata <- seq(from = -5, to = 5, length = n)
  ydata <- 1/(1 + xdata^2)
  yinterp <- lagrangeInterp(xdata, ydata, x)
  lines(x, yinterp)
  error <- c(error, max(abs(y - yinterp)))
}
```

![](coursebook_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

``` r
plot(nvec, log10(error), xlab = "n", ylab = "log10 of error")
```

![](coursebook_files/figure-gfm/unnamed-chunk-33-2.png)<!-- -->

Not good! The error goes up as we take more and more points.
Equally-spaced nodes are very natural in many applications (scientific
measurements, audio/video signals, etc.). But sadly, it turns out that
in some circumstances, inteprolating with polynomials through
equally-spaced nodes leads to very undesirable oscillations like those
above, called **Runge’s phenomenon**.

## Interpolation Error

You might remember that for an $n-1$st degree Taylor series (that is,
$n$ coefficients) of $f(x)$ centered around $x=x_0$, the magnitude of
the error term is

$$
\left|\frac{f^{(n)}(c)}{n!}(x-x_0)^{n}\right|
$$

where $c$ is a number between $x$ and $x_0$. For an interpolating
polynomial constructed from $n$ points on an interval $[x_1,x_n]$, the
interpolation error is

$$
\left|\frac{f^{(n)}(c)}{n!}(x-x_1)(x-x_2)\cdots (x-x_n)\right|
$$

for some $c \in [x_1,x_n]$. If you are interested in the proof, I am
happy to go over it with you. Honestly, though, I find the proof to be
rather nonintuitive and unenlightening.

Let’s use the error term to do an example, where we again ask the
question HOW DO WE EVEN KNOW ANYTHING? For instance, suppose you want to
calculate the value of an exponential function. [This turns out to be
hard](https://math.stackexchange.com/questions/1239352/how-do-pocket-calculators-calculate-exponents)!
Let’s restrict the problem a bit. Suppose we want to calculate
$\mathrm{e}^{x}$ on $[0,1]$ with 5 digits of accuracy. One way to do
this would be to choose equally-spaced points on $[0,1]$ and make a
table of $\mathrm{e}^{x}$ for those values using some very accurate
method. Then, we let the user choose a value of $x$. We find the two
values that surround it in our table, construct the linear interpolating
polynomial between them, and use it to estimate our answer. For
instance, if we chose $x$ values of $0.1, 0.2, \ldots, 1$ and the user
input $x=0.3268$, we would construct the inteprolating polynomial
through $(0.3,\mathrm{e}^{0.3})$ and $(0.4,\mathrm{e}^{0.4})$ and then
plug in our $x$. Let’s ask the question: to achieve 5-digit accuracy
with this method, how many points must we take on $[0,1]$? (Note:
crucially, this is not at all the same as creating a single
interpolating polynomial that goes through all the points.) We have

$$\begin{align*}
|f(x) - P_1(x)| &= \left|\frac{f^{(n)}(c)}{n!}(x-x_1)(x-x_2)\cdots (x-x_n)\right|\ \\
&= \left|\frac{\left(\mathrm{e}^{x}\right)^{''}|_c}{2!}(x-x_1)(x-x_2)\right| \\ \ 
& = \left|\frac{\mathrm{e}^{c}}{2} (x-x_1)(x-x_1-h)\right|
\end{align*}$$

Here, I’ve called the spacing between the two points $h$, that is
$x_2 = x_1 + h$. Now, we can ask what is the worst (biggest) that the
term $\left|(x-x_1)(x-x_1-h)\right|$ can be, and via calculus, we can
show it is $h^2/4$. Also, on our interval of interest, we know that
$\mathrm{e}^{c} \leq \mathrm{e}^{1}$. Therefore, we can write

$$
|f(x) - P_1(x)| \leq \frac{\mathrm{e}h^2}{8}.
$$

For our desired level of accuracy, we need

$$
\frac{\mathrm{e}h^2}{8} \leq 0.5 \times 10^{-5}.
$$

Solving for $h$, we find $h \lessapprox 3.8 \times 10^{-3}$. That is the
space between points, so on $[0,1]$ this corresponds to just over 260
points.

## Chebyshev Nodes

Look at the error expression again. There are three parts of it.
Factorial is the same for every function. If the derivative gets smaller
with more $n$ then we’re good. If it gets bigger, the error term could
potentially grow overall with $n$. Since there’s nothing we can do about
the derivative (we can’t choose it!) we should try to do what we can to
minimize what we can control, namely the remaining part of the error
term. Chebyshev nodes have special x-values that do this! Much calculus
leads us to the following result. For $n$ points on $[0,1]$, we should
choose nodes

$$
x_i = \cos \frac{(2i-1)\pi}{2n}, \quad i = 1,\ldots,n
$$

and the maximum magnitude of the relevant portion of the error term is

$$
\max_{-1\leq x\leq 1}\left|\prod_{i=1}^n (x-x_i)\right| = \frac{1}{2^{n-1}}.
$$

If we generalize these results to the interval $[a,b]$, we find

$$
x_i = \frac{b+a}{2} + \frac{b-a}{2}\cos \frac{(2i-1)\pi}{2n}, \quad i = 1,\ldots,n
$$

with

$$
\max_{a\leq x\leq b}\left|\prod_{i=1}^n (x-x_i)\right| = \frac{\left(\frac{b-a}{2}\right)^n}{2^{n-1}}.
$$

The proof is not direct, which is why I have eliminated it here.

## Comparing Interpolation Methods

Let’s do an example comparing interpolation approaches. Consider
$f(x)=(1/\sqrt{2\pi})\mathrm{e}^{-x^2/2}$. This is the standard normal
distribution which plays a key role in probability and statistics.

One thing to know is that the derivatives of this function grow with
$n$.

``` r
n <- 0:10
maxderiv <- c(0.398942,0.241971,0.178032, 0.550588,1.19683,2.30711,4.24061,14.178,41.8889,115.091,302.425) # Computed in Mathematica
plot(n, log(maxderiv))
```

![](coursebook_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

Given this, using equally-spaced nodes seems reckless, but we can try it
anyway, say, on $[-10,10]$ with 30 data points to start with.

``` r
a <- -10
b <- 10
x <- seq(from = a, to = b, length = 10000)
f <- function(x){exp(-x^2/2)/sqrt(2*pi)}
y <- f(x)
plot(x, y, type = "l", lwd = 3, xlim = c(a,b), ylim = c(-0.2,0.5))
n <- 30
xdataequal <- seq(from = a, to = b, length = n)
ydataequal <- lagrangeInterp(xdataequal, f(xdataequal), x)
lines(x, ydataequal, col = "red", lwd = 2)
```

![](coursebook_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

Ok, that approach is not going to work! Let’s try the approach of a
lookup table with linear interpolation. Let’s suppose we wish to achieve
4 digit accuracy.

``` r
n <- 1
error <- Inf
while (error > 0.5e-4){
  n <- n+1
  xlookup <- seq(from = a, to = b, length = n)
  ylookup <- f(xlookup)
  ytable <- approx(xlookup, ylookup, x)$y
  error <- max(abs(y - ytable))
}
nlookup <- n
print(nlookup)
```

    ## [1] 632

Now we can go back to polynomial interpolation with Chebyshev nodes.

``` r
n <- 30
error <- Inf
while (error > 0.5e-4){
  n <- n+1
  odds <- seq(from = 1, to = 2*n-1, by = 2)
  xdatacheb <- sort((b + a)/2 + (b - a)/2 * cos(odds*pi/2/n))
  ydatacheb <- lagrangeInterp(xdatacheb, f(xdatacheb), x)
  error <- max(abs(y - ydatacheb))
}
ncheb <- n
print(ncheb)
```

    ## [1] 41

This is an improvement in compression by a factor of 632/41 =
15.4146341.

## Cubic Splines

I’ve tried to convince you that it can be problematic to construct
interpolating polynomials of high degree. When dealing with a lot of
data, an alternative approach can be to construct low degree
interpolating polynomials through successive sets of points. Typically
we use cubics, and these are called cubic splines.

Before progressing to real data later on, let’s do an illustrative
example with a small amount of fake data. We make some data points and
connect them with linear splines. We can do this using the built-in
`approxfun` command which returns a function representing the linear
spline.

``` r
x <- c(-2,-1.5,-1,0.25,1,2,3.75,4,5)
y <- c(4,4.2,3,5,0,-2,2,1,1)
xplot <- seq(from = -2, to = 5, length = 200)
linearspline <- approxfun(x, y)
plot(x, y, ylim = c(-2.5,5.5))
lines(xplot, linearspline(xplot), col = "red", lwd = 2)
```

![](coursebook_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

Your eyeball might be telling you that this is a very jagged graph. Most
things in nature and society are not this jagged, so it might feel
desirable to represent the data with smoother functions. Let me show you
what this looks like.

``` r
cubicspline <- splinefun(x, y)
plot(x, y, ylim = c(-2.5,5.5))
lines(xplot, linearspline(xplot), col = "red", lwd = 2)
lines(xplot, cubicspline(xplot), col = "green", lwd = 2)
```

![](coursebook_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

This is smoother. What smoothness means here is continuity of
derivatives from one spline to the next. Using cubic rather than linear
splines gives us more coefficients, and using these coefficients we can
make the derivatives of successive splines match up. How exactly does
this work though?

## Mathematical Conditions for Cubic Splines

Suppose we have data points $(x_1,y_1),\ldots,(x_n,y_n)$. We will
connect successive points with cubic curves. With $n$ points, there are
$n-1$ curves. Each curve is cubic, having the form

$$
S_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3.
$$

With $n-1$ curves each having four coefficients, there are $4(n-1)$
coefficients.

First, we force each spline to pass through its two endpoints (which
also makes the overall spline curve continuous). By inspection, this
forces $a_i=y_i$, and it also enforces a condition on the relationship
between $b_i, c_i, d_i$. So doing the bookkeeping, that is $2(n-1)$
conditions to enforce, which leaves $2n-2$ coefficients undetermined.

Next, we force continuity of the first derivatives of successive
splines. Note that at $x_1$ there’s no condition to enforce, since
there’s no spline to the left of it. And at $x_n$ there’s no condition
to enforce, since there’s no spline to the right of it. Therefore, we
enforce conditions at the $n-2$ points $x_2,\ldots,x_{n-1}$.
Specifically, the condition is $S_i^'(x_i)=S_{i-1}^'(x_i)$. Subtracting
these $n-2$ conditions from our previous count of $2n-2$, there are $n$
conditions left.

Finally, we enforce continuity of the second derivatives of successive
splines. This is very similar to enforcing the previous condition. The
condition is $S_i^{''}(x_i)=S_{i-1}^{''}(x_i)$ at the points
$x_2,\ldots,x_{n-1}$. Substracting these $n-2$ conditions from our
previous total of $n$ conditions, there are two conditions left.

In short, thus far, the procedure for finding spline coefficients is
underetermined because there are $4n-4$ unknowns but only $4n-2$
equations. To make the system solvable, we have to make a choice about
what conditions to enforce at the leftmost and rightmost points. Some of
the most common choices are:

- **Natural spline**. The concavity at the left and right endpoints is
  zero, that is $S_1^{''}(x_1)=S_{n-1}^{''}(x_n)=0$.
- **Clamped sline**. The concavity at the left and right endpoints is
  set to a user specified value, that is, $S_1^{''}(x_1)=m_1$ and
  $S_{n-1}^{''}(x_n)=m_n$.
- **FMM (not-a-knot)**. $S_1=S_2$ is a single cubic equation that is run
  through the first 3 points, and $S_{n-2}=S_{n-1}$ is, similarly, a
  single cubic run through the last 3 points.

Crucially, since we find the coefficients by solving an square linear
system of dimension $4n-4$, these boundary conditions don’t merely
affect the first and last splines, but influence ALL of the splines.

You don’t need to memorize the details of these different types of
splines. My main goals for you are to understand what they mean and to
be able to implement them in R.

## Splines and Linear Algebra

As we have been discussing, to find spline coefficients, we have to
solve a linear system. I won’t write down the whole system here because
we are going to use built-in tools to solve it. However, it’s good
mathematical literacy to know that by writing down the system of
equations, you can see that it is tridiagonal and strictly diagonally
dominant, which are nice numerical properties.

## Implementing Cubic Spline Interpolation

Just to emphasize how splines avoid the problem of high-degree
polynomial interpolation, let’s do a cooked example.

``` r
set.seed(123)
n <- 30
x <- sort(runif(n))
y <- cumsum(abs(rnorm(n)))
plot(x, y, pch = 19, ylim = c(-2,35), cex = 1.5)
xx <- seq(from = min(x), to = max(x), length = 1000)
yy <- lagrangeInterp(x, y, xx)
lines(xx, yy, col = "red", lwd = 3)
cubicspline <- splinefun(x, y, method = "natural")
lines(xx, cubicspline(xx), col = "blue", lwd = 3)
```

![](coursebook_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

Now let’s work with some real data, let’s say, Tesla stock price for the
last 100 days (we won’t actually get 100 days because of days when the
markets were closed, including weekends).

``` r
# Set dates and stock symbol
first.date <- Sys.Date() - 100
last.date <- Sys.Date()
tickers <- c('TSLA')
# Acquire data
l.out <- BatchGetSymbols(tickers = tickers, first.date = first.date, last.date = last.date)
```

    ## Warning: `BatchGetSymbols()` was deprecated in BatchGetSymbols 2.6.4.
    ## ℹ Please use `yfR::yf_get()` instead.
    ## ℹ 2022-05-01: Package BatchGetSymbols will soon be replaced by yfR.  More
    ##   details about the change is available at github
    ##   <<www.github.com/msperlin/yfR> You can install yfR by executing:
    ## 
    ## remotes::install_github('msperlin/yfR')

    ## 
    ## Running BatchGetSymbols for:
    ##    tickers =TSLA
    ##    Downloading data for benchmark ticker
    ## ^GSPC | yahoo (1|1) | Not Cached | Saving cache
    ## TSLA | yahoo (1|1) | Not Cached | Saving cache - Got 100% of valid prices | Youre doing good!

``` r
price <- l.out[[2]]$price.close
day <- 1:length(price)
# Sample every 5th day
sampledday <- day[seq(from = 1, to = length(price), by = 5)]
sampledprice <- price[seq(from = 1, to = length(price), by = 5)]
# Fit interpolating polynomial
interpolatedprice <- lagrangeInterp(sampledday, sampledprice, day[1:max(sampledday)])
ymin <- min(price) - 0.1*(max(price) - min(price))
ymax <- max(price) + 0.1*(max(price) - min(price))
plot(day[1:max(sampledday)], interpolatedprice, col = "red", type = "l", ylim = c(ymin,ymax))
points(day, price, col = "blue")
```

![](coursebook_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

``` r
# Fit splines
TSLAspline <- splinefun(sampledday, sampledprice, method = "natural")
plot(day, TSLAspline(day), col = "blue", type = "l", ylim = c(ymin,ymax))
points(day, price, col = "blue")
```

![](coursebook_files/figure-gfm/unnamed-chunk-41-2.png)<!-- -->

# Least Squares I

## Big Picture

All of our study of interpolation has been based on the idea that the
model (a polynomial or polynomnial spline) should precisely pass through
every data point. But what if we loosen this restriction? What if the
data is thought to have error or noise that we don’t want to represent
in our model? Or what if we just don’t care about precise interpolation,
and decide approximation even at the data points themselves is good
enough?

## Goals

- Explain how least squares arises as a model fitting problem
- Explain how the solution of least squares involves vector projection
- Set up and solve least squares problems
- Define, identify, and calculate crucial quantities such as basis,
  target, normal equations, least squares solution, residual,
  pseudoinverse, projection, projection operator
- Recognize when a least squares approach is appropriate
- Explain how least squares allows data compression

## Model Fitting

By way of motivation, let’s examine a pedagogical data set. Suppose $a$
represents the amount of money (in \$1,000’s) a company spent on
advertising during different quarters, and $s$ represents money the
company earned on sales that quarter. We can plot the data to explore
it.

``` r
a <- c(3,4,5,6)
s <- c(105,117,141,152)
plot(a, s, xlab = "advertising", ylab = "sales")
```

![](coursebook_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

The company would like to model this data so they can predict sales for
other levels of advertising. The data looks roughly linear, and we have
no reason to expect a complicated relationship, so let’s try modeling
the data with a line, $s = x_0 + x_1 a$ where $x_{0,1}$ are unknown
coefficients. Plugging in to the model, we find

$$ \begin{aligned}
x_0 + 3x_1 &= 105\\ x_0 + 4x_1 &= 117\\ x_0 + 5x_1 &= 141\\ x_0 + 6x_1 &= 152
\end{aligned} $$

which we can write in vector form as

$$
x_0 \begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix} + x_1 \begin{pmatrix} 3 \\ 4 \\ 5 \\ 6 \end{pmatrix} = \begin{pmatrix} 105 \\ 117 \\ 141 \\ 152 \end{pmatrix}.
$$

By writing it this way, we can remember one interpretation of linear
systems. In this case we should imagine a four-dimensional space in
which we are trying to reach a particular target (the right hand side
vector) by taking linear combinations of the two vecors on the left.
Unfortunately, those two vectors only span a two dimensional subspace,
so the chance that we can make it to our target vector are pretty slim.
Stated differently: the sysem is overdetermined. But we would still like
to find a good model, so what should we do?

## Projection onto a Vector

To examine the details, let’s start with an even more fundamental
example: a single vector in the plane. Suppose I hand you the vector
$\mathbf{A} = (2,1)^T$ and tell you to use it to reach the target vector
$\mathbf{b}=(6,8)^T$. Well, you can’t do it exactly because there is no
scalar $x$ such that $x\mathbf{A} = \mathbf{b}$. So let us do the next
best thing: let’s find the value of $x$ such that $x\mathbf{A}$ is as
close as possible to $\mathbf{b}$. We can draw a picture to solve this
problem.

![](coursebook_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

Where should we stop on the dotted line? When we are perpendicular to
the end of $\mathbf{b}$. This results in the following picture.

![](coursebook_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

From this picture, two relationships arise:

$$
\begin{align}
x\mathbf{A} + \mathbf{r}&=\mathbf{b}\\
\mathbf{A} \cdot \mathbf{r} &=0.
\end{align}
$$

Let us dot the first equation with the vector $\mathbf{A}$. Then we have

$$
\mathbf{A} \cdot \mathbf{A}x + \mathbf{A} \cdot \mathbf{r} = \mathbf{A} \cdot \mathbf{b}.
$$

From the second equation above, we can eliminate the second term on the
left hand side to write

$$
\mathbf{A}^T \mathbf{A} x = \mathbf{A}^T \mathbf{b}
$$

where we have used the fact that
$\mathbf{y}\cdot\mathbf{z}=\mathbf{y}^T\mathbf{z}$. We can solve for
$\mathbf{x}$ by writing
$x = \left(\mathbf{A}^T \mathbf{A}\right)^{-1}\mathbf{A}^T \mathbf{b}$.
We can also calculate the vector that was as close as possible to
$\mathbf{b}$. We will call it $\widehat{\mathbf{b}}$ and it is

$$
\begin{align}
\widehat{\mathbf{b}} &= \mathbf{A}x \\
&=\mathbf{A}\left(\mathbf{A}^T \mathbf{A}\right)^{-1}\mathbf{A}^T \mathbf{b}\\
& = \mathbf{P} \mathbf{b}
\end{align}
$$

where the last equation defines a new quantity that we call
$\mathbf{P}$.

Let us know revisit what we have done and emphasize/introduce some
vocabulary. We started with a vector $\mathbf{A}$ that we used as a
**basis** to try to reach the **target** $\mathbf{b}$. We couldn’t do it
exactly, so we calcualted the closest we could come to $\mathbf{b}$,
which turned out to be $\widehat{\mathbf{b}}$. This is called the
**projection** of $\mathbf{b}$ into the subspace spanned by
$\mathbf{A}$. We found $\widehat{\mathbf{b}} = \mathbf{A}x$, where $x$
is called the **least squares solution**, which solved the **normal
equations** $\mathbf{A}^T \mathbf{A} x = \mathbf{A}^T \mathbf{b}$. We
can summarize the calculation of $x$ by remembering
$x = \left(\mathbf{A}^T \mathbf{A}\right)^{-1}\mathbf{A}^T \mathbf{b}$
where $\left(\mathbf{A}^T \mathbf{A}\right)^{-1}\mathbf{A}^T$ is called
the **pseudoinverse** of $\mathbf{A}$. Also, we can summarize the
calculation of $\widehat{\mathbf{b}}$ as
$\widehat{\mathbf{b}} = \mathbf{P} \mathbf{b}$ where
$\mathbf{P} = \mathbf{A} (\mathbf{A}^T\mathbf{A})^{-1} \mathbf{A}^T$ is
what we call a **projection operator** or a **projection matrix**. Since
we didn’t succeed in reaching $\mathbf{b}$, there is some error, and we
call this the **residual**,
$\mathbf{r} = \mathbf{b}-\widehat{\mathbf{b}}$.

What are the words/ideas you should make sure you understand in the
narrative above?

- basis
- target
- normal equations
- least squares solution
- pseudoinverse
- projection
- projection operator
- residual

Let’s calculate this concretely in `R`.

``` r
a <- c(2,1)
b <- c(6,8)
pseudoinv <- solve(t(a) %*% a) %*% t(a)
x <- pseudoinv %*% b
bhat <- a %*% x
r <- b - bhat
print(x)
```

    ##      [,1]
    ## [1,]    4

``` r
print(bhat)
```

    ##      [,1]
    ## [1,]    8
    ## [2,]    4

``` r
t(a) %*% r # Dot product
```

    ##      [,1]
    ## [1,]    0

## Projection onto a Plane

Let’s apply these same conceps to our original problem of predicting
sales from advertising. We had

$$
x_0 \begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix} + x_1 \begin{pmatrix} 3 \\ 4 \\ 5 \\ 6 \end{pmatrix} = \begin{pmatrix} 105 \\ 117 \\ 141 \\ 152 \end{pmatrix}
$$

which we’ll write symbolically as

$$
x_0 \mathbf{A}_0 + x_1 \mathbf{A}_1 = \mathbf{b}.
$$

We’d like to reach $\mathbf{b}$ using the basis vectors
$\mathbf{A}_{0,1}$, but we can’t, so let’s consider getting as close as
possible. The picture looks something like this.

![](PlaneProjection.png)

This picture gives rise to the equation

$$
\begin{align}
x_0 \mathbf{A}_0 + x_1 \mathbf{a_1} + \mathbf{r} &= \mathbf{b} \\
a_0 \cdot \mathbf{r} &= 0 \\
a_1 \cdot \mathbf{r} &= 0.
\end{align}
$$

Let’s take the dot product of the first equaton with each basis vector.
We find

$$
\begin{align}
x_0 \mathbf{A}_0\cdot\mathbf{A}_0 + x_1 \mathbf{A}_0\cdot\mathbf{A}_1 &= \mathbf{A}_0 \cdot \mathbf{b}\\
x_0 \mathbf{A}_1\cdot\mathbf{A}_0 + x_1 \mathbf{A}_1\cdot\mathbf{A}_1 &= \mathbf{A}_1 \cdot \mathbf{b}
\end{align}
$$

where we have elminated terms that turn out to be zero thanks to the
second and third equations previously. Note that aa matrix way to write
this. We can write

$$
\begin{pmatrix} \mathbf{A}_0^T \\ \mathbf{A}_1^T \end{pmatrix} \begin{pmatrix} \mathbf{A}_0 & \mathbf{A}_1 \end{pmatrix} \begin{pmatrix} x_0 \\ x_1 \end{pmatrix} = \begin{pmatrix} \mathbf{A}_0^T \\ \mathbf{A}_1^T \end{pmatrix} \mathbf{b}.
$$

If we let $\mathbf{A}$ represent the matrix with columns
$\mathbf{A}_{0,1}$ and if we let $\mathbf{x}=(x_0,x_1)^T$ then we can
write the last equation as

$$
\mathbf{A}^T\,\mathbf{A}\,\mathbf{x} = \mathbf{A}^T\,\mathbf{b}.
$$

These are the normal equations. For us, concretely, it looks like

$$
\begin{align}
\begin{pmatrix} 1 & 1 & 1 & 1\\ 3 & 4 & 5 & 6\end{pmatrix} \begin{pmatrix} 1 & 3\\1 & 4\\ 1 & 5\\ 1& 6\end{pmatrix}\begin{pmatrix}x_1 \\ x_2 \end{pmatrix} &= \begin{pmatrix} 1 & 1 & 1 & 1\\ 3 & 4 & 5 & 6\end{pmatrix} \begin{pmatrix} 105 \\ 117 \\ 141 \\ 152\end{pmatrix}\\
\begin{pmatrix} 4 & 18 \\ 18 & 86 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} &= \begin{pmatrix} 515 \\ 2400 \end{pmatrix}
\end{align}
$$

which has solution

$$
\begin{pmatrix} x_1 \\ x_2  \end{pmatrix} = \begin{pmatrix}  54.5 \\ 16.5 \end{pmatrix}.
$$

Let’s calculuate this in `R`.

``` r
a0 <- c(1,1,1,1)
a1 <- c(3,4,5,6)
A <- cbind(a0,a1)
b <- c(105,117,141,152)
pseudoinv <- solve(t(A) %*% A) %*% t(A)
x <- pseudoinv %*% b
bhat <- A %*% x
r <- b - bhat
print(x)
```

    ##    [,1]
    ## a0 54.5
    ## a1 16.5

``` r
print(bhat)
```

    ##       [,1]
    ## [1,] 104.0
    ## [2,] 120.5
    ## [3,] 137.0
    ## [4,] 153.5

``` r
t(r) %*% a0
```

    ##                   [,1]
    ## [1,] 1.84741111298e-13

``` r
t(r) %*% a1
```

    ##                   [,1]
    ## [1,] 8.81072992343e-13

``` r
plot(a1, b, xlab = "advertising", ylab = "sales")
xx <- seq(from = 0, to = 6, length = 200)
lines(xx, horner(rev(x), xx)$y) # Note we need to reverse order of coefficients
```

![](coursebook_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->

Symbolically, we calculated the least squares solution
$\mathbf{x} = \left(\mathbf{A}^T \mathbf{A}\right)^{-1}\mathbf{A}^T \mathbf{b}$
where $\left(\mathbf{A}^T \mathbf{A}\right)^{-1}\mathbf{A}^T$ is the
pseudoinverse. The projection is
$\widehat{\mathbf{b}} = \mathbf{A} \mathbf{x} = \mathbf{A} \left(\mathbf{A}^T \mathbf{A}\right)^{-1}\mathbf{A}^T \mathbf{b}$
where $\mathbf{A} \left(\mathbf{A}^T \mathbf{A}\right)^{-1}\mathbf{A}^T$
is the projection matrix.

Though our example here has used merely two basis vectors, the ideas
extend to any number.

## Model Fitting

The technique we have developed here of solving a least squares problem
by vector projection works in a model fitting context whenever the
unknowns $\mathbf{x}$ appear linearly. It doesn’t require that other
quantities appear linearly.

For instance, suppose that we wanted to fit a parabola to the data
$(0,6)$, $(1,5)$, $(2,2)$, and $(3,2)$. We choose the model
$y = c_0 + c_1 x + c_2 x^2$ where our unknown is the vector
$\mathbf{c}=(c_0,c_1,c_2)^T$. Though $x$ appears nonlinearly, the
coefficients $c_{0,1,2}$ don’t, so we are fine! Let’s solve this problem
using linear algebra. We want to solve

$$
\begin{align}
c_0 + c_1 \cdot 0 + c_2 \cdot 0^2 &= 6 \\
c_0 + c_1 \cdot 1 + c_2 \cdot 1^2 &= 5 \\
c_0 + c_1 \cdot 2 + c_2 \cdot 2^2 &= 2 \\
c_0 + c_1 \cdot 3 + c_2 \cdot 3^2 &= 2
\end{align}
$$

So, we calculuate in `R`:

``` r
x <- c(0,1,2,3)
A <- vander(x)[,2:4]
b <- c(6,5,2,12)
pseudoinv <- solve(t(A) %*% A) %*% t(A)
c <- pseudoinv %*% b
bhat <- A %*% c
r <- b - bhat
print(c)
```

    ##       [,1]
    ## [1,]  2.75
    ## [2,] -6.75
    ## [3,]  6.75

``` r
plot(x,b)
xx <- seq(from = 0, to = 3, length = 200)
lines(xx, horner(c, xx)$y)
```

![](coursebook_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->

``` r
print(r)
```

    ##       [,1]
    ## [1,] -0.75
    ## [2,]  2.25
    ## [3,] -2.25
    ## [4,]  0.75

If model parameters don’t appear in a linear fashion, sometimes you can
transform the equation so that they do. For example, if you wanted to
fit data to $y = C \mathrm{e}^{k\,x}$, you could take the log of both
sides to obtain $\ln y = \ln C + k\,x$. By considering the data
$(x,\,\ln y)$ you could take a least squares approach to find $\ln C$
and $k$.

Overall, here’s the process one would follow for model fitting in cases
where the model parameters appear in a lienar fashion.

1.  Look at data.
2.  Propose a model.
3.  Force the model, resulting in $\mathbf{A}\,\mathbf{x}=\mathbf{b}$.
4.  Solve the normal equations
    $\mathbf{A}^T\,\mathbf{A}\,\mathbf{x}=\mathbf{A}^T\,\mathbf{b}$.
5.  Assess the fit of the model visually and/or using the residual
    vector.

A statistics class would provide much more sophisticated ways of
analyzing the error.

## Least Squares and Data Compression

Adopting a least squares approach allows, potentially, massive
compression of data. Suppose we had 100 points that looked like this

``` r
x <- seq(from = 0, to = 1, length = 100)
y <- 3*x + 0.4*(2*runif(100) - 1)
plot(x, y)
```

![](coursebook_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->

If we decided to represent this data with a line, we’d go down from
having 100 pieces of information (the original data points) to merely 2
(a slope and an intercept).

## An Optimization Viewpoint on Least Squares

So far, we’ve taken a geometric approach to solving least-squares
problems, but there’s another way to get the same result: optimization.
We can start directly with the square of the norm of the residual vector
(a factor of 1/2 is included for algebraic convenience, but it doesn’t
change the result):

$$
\frac{1}{2}||\mathbf{r}||^2 = \frac{1}{2}||\mathbf{A}\mathbf{x}-\mathbf{b}||^2.
$$

More formally, we can define the scalar **objective function**

$$
f(\mathbf{x}) = \frac{1}{2} ||\mathbf{A}\,\mathbf{x}-\mathbf{b}||^2
$$

and define our least squares solution $\mathbf{x}_{L\,S}$as the value of
$\mathbf{x}$ that minimizes this objective function, that is,

$$
\mathbf{x}_{L\,S} = \mathop{\mathrm{arg\,min}}_{\mathbf{x}} f(\mathbf{x}).
$$

Now you know where the terminology least squares comes from! We can
apply calculus to minimize this expression and the normal equations will
result. This is a problem on your in-class exercise and your homework,
and it is critical that you work through the details.

# Least Squares II

## Big picture

We’ve talked about how to fit a model to data using least squares, and
we’ve examined the theoretical aspect of this process, but thus far
we’ve ignored computational issues. It turns out that solving the normal
equations can be difficult on a computer, and a technique called **QR**
factorization provides a potentially better way to solve a least squares
problem.

## Goals

- Define orthonormal matrices
- Perform Graham-Schmidt Orthogonalization
- Solve least squares problems with **QR** decomposition

## Orthogonality is nice

The eventual goal of this lesson is to show you how to solve the least
squares problem

$$
\mathbf{x}_{L\,S} = \mathop{\mathrm{arg\,min}}_x \frac{1}{2}||\mathbf{A}\,\mathbf{x}-\mathbf{b}||^2
$$

by writing the matrix $\mathbf{A}$ in a convenient way. But to build up
to that, we need to introduce a number of ideas. First off:
orthogonality.

You might remember that two vectors $\mathbf{v}$ and $\mathbf{w}$ are
**orthogonal** if
$\mathbf{v} \cdot \mathbf{w} = \mathbf{v}^T \mathbf{w} = 0$. We say that
two vectors are **orthonormal** if they are orthogonal and each have
norm of one, that is
$\mathbf{v} \cdot \mathbf{v} = \mathbf{w} \cdot \mathbf{w} = 1$.

A matrix is called **orthogonal** or **orthonormal** (they are used
interchangeably for matrices, sometimes) if its columns are orthonormal
vectors. A cool property arises from this. Suppose these orthonormal
columns are $\mathbf{v}_1,\ldots,\mathbf{v}_n$. Then we can consider the
quantity $\mathbf{A}^T \mathbf{A}$:

$$
\begin{align}
\mathbf{A}^T \mathbf{A} &= \begin{pmatrix} \mathbf{v}_1^T\\ \mathbf{v}_2^T \\ \vdots \\ \mathbf{v}_n^T \end{pmatrix} \begin{pmatrix} \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n \end{pmatrix} \\
&= \begin{pmatrix} \mathbf{v}_1^T \mathbf{v}_1 & \mathbf{v}_1^T \mathbf{v}_2 & \cdots & \mathbf{v}_1^T \mathbf{v}_n \\ \mathbf{v}_2^T \mathbf{v}_1 & \mathbf{v}_2^T \mathbf{v}_2 & \cdots &  \mathbf{v}_2^T \mathbf{v}_n \\
\vdots & \vdots & \cdots & \vdots \\
\mathbf{v}_n^T \mathbf{v}_1 & \mathbf{v}_n^T \mathbf{v}_2 & \cdots & \mathbf{v}_n^T \mathbf{v}_n
\end{pmatrix}\\
&= \begin{pmatrix} 1 & 0 & \cdots & 0 \\ 0 & 1 & \cdots & 0\\ \vdots & \vdots & \cdots & \vdots \\ 0 & 0 & \cdots & 1\end{pmatrix}\\
&= \mathbf{I}_n.
\end{align}
$$

If $\mathbf{A}$ is square and orthonormal, since
$\mathbf{A}^T\mathbf{A} = \mathbf{I}$, then by definition of inverse,
$\mathbf{A}^T = \mathbf{A}^{-1}$. This is a pretty great way to
calculate an inverse! Remember, though, that in our least squares
context , $\mathbf{A}$ is generally *not* a square matrix.

Another helpful property of orthogonal matrices is that they preserve
length, meaning that have norm of one. That is, suppose $\mathbf{Q}$ is
orthogonal. Then

$$
\begin{align}
||\mathbf{Q} \mathbf{v}||^2 &= \mathbf{Q} \mathbf{v} \cdot \mathbf{Q} \mathbf{v} \\
& = (\mathbf{Q}\mathbf{v})^T (\mathbf{Q}\mathbf{v})\\
& = \mathbf{v}^T \mathbf{Q}^T \mathbf{Q} \mathbf{v}\\
& = \mathbf{v}^T \mathbf{Q}^{-1} \mathbf{Q} \mathbf{v}\\
& = \mathbf{v}^T \mathbf{v}\\
& = ||\mathbf{v}||^2.
\end{align}
$$

Since $||\mathbf{Q}\,\mathbf{v}||^2 = ||\mathbf{v}||^2$, we know
$||\mathbf{Q}\,\mathbf{v}|| = ||\mathbf{v}||$ and therefore
$||\mathbf{Q}||=1$. A similar calculuation shows that
$||\mathbf{Q}^{-1}||=1$. And then, by definition of condition number,
$\kappa(\mathbf{Q}) = 1$. This means that orthogonal matrices are
incredibly well-conditioned.

Finally, let’s consider projecting a vector into a subspace spanned by
orthonormal vectors $\mathbf{q}_i$, $i=1,\ldots,n$. Define $\mathbf{A}$
as

$$
\mathbf{A} = \begin{pmatrix} | & | &  & | \\ \mathbf{q}_1 & \mathbf{q}_2 & \cdots & \  \mathbf{q}_n \\| & | &  & | \\ \end{pmatrix}.
$$

Let’s project a vector $\mathbf{w}$ onto the subspace spanned by the
columns of $\mathbf{A}$. By definition of the projection operator, and
using the fact that $\mathbf{A}$ is orthogonal,

$$
\begin{align}
\mathbf{P}\mathbf{w} &= \mathbf{A} (\mathbf{A}^T \mathbf{A})^{-1} \mathbf{A}^T \mathbf{w}\\
&= \mathbf{A} \mathbf{A}^{-1} \mathbf{A} \mathbf{A}^T \mathbf{w} \\
&= \mathbf{A} \mathbf{A}^T \mathbf{w} \\
&= \mathbf{q}_1 \mathbf{q}_1^T \mathbf{w} + \mathbf{q}_2 \mathbf{q}_2^T \mathbf{w} + \cdots + \mathbf{q}_n \mathbf{q}_n^T \mathbf{w}.
\end{align}
$$

So, when projecting onto the space spanned by orthonormal vectors, you
can just project onto each vector separately and add up the results. In
the calculation above, remember that $\mathbf{q}_i \mathbf{q}_i^T$ is a
matrix.

As an example, take

$$
\mathbf{q}_1 = \left(\frac{1}{3},\frac{2}{3},\frac{2}{3}\right)^T, \quad \mathbf{q}_2 = \left(\frac{2}{15},\frac{2}{3},-\frac{11}{15}\right)^T
$$

so that

$$
\mathbf{A} = \begin{pmatrix} \frac{1}{3} & \frac{2}{15} \\ \frac{2}{3} & \frac{2}{3} \\ \frac{2}{3} & -\frac{11}{15} \end{pmatrix}.
$$

Then

$$
\begin{align}
\mathbf{A} \mathbf{A}^T &= \begin{pmatrix} \frac{1}{3} & \frac{2}{15} \\ \frac{2}{3} & \frac{2}{3} \\ \frac{2}{3} & -\frac{11}{15} \end{pmatrix}
\begin{pmatrix} \frac{1}{3} &  \frac{2}{3} & \frac{2}{3} \\ \frac{2}{15} & \frac{2}{3} & -\frac{11}{15} \end{pmatrix} \\
&= \frac{1}{225}
\begin{pmatrix}
29 & 70 & 28 \\
70 & 200 & -10 \\
28 & -10 & 221
\end{pmatrix}\\
& = \underbrace{
\begin{pmatrix}
\frac{1}{9} & \frac{2}{9} & \frac{2}{9} \\
\frac{2}{9} & \frac{4}{9} & \frac{4}{9} \\
\frac{2}{9} & \frac{4}{9} & \frac{4}{9}
\end{pmatrix}}_{\mathbf{q}_1 \mathbf{q}_1^T}  +  \underbrace{
\begin{pmatrix}
\frac{4}{225} & \frac{4}{45} & -\frac{22}{225} \\
\frac{4}{45} & \frac{4}{9} & -\frac{22}{45} \\
-\frac{22}{225} & -\frac{22}{45} & \frac{121}{225}
\end{pmatrix}}_{\mathbf{q}_2 \mathbf{q}_2^T}.
\end{align}
$$

## Gram-Schmidt orthogonalization

Now that I’ve convinced you that it’s nice to have an orthonormal basis,
let’s talk about how to get one. Gram-Schmidt orthogonalization is an
iterative way of creating an orthonormal basis from a set of vectors.

Let’s see how it works via an example. Suppose we have
$\mathbf{v}_1 = (3,4,0)^T$ and $\mathbf{v}_2 = (1,1,0)^T$.

``` r
v1 <- c(1,2,2)
v2 <- c(2,1,-2)
v3 <- c(1,1,0)
```

Step 1. Take the first vector and turn it into a unit vector.

``` r
y1 <- v1
r11 <- Norm(y1, 2)
q1 <- y1 / r11
```

Step 2. Think of $\mathbf{v}_2$ as made up of stuff in the subspace
spanned by $\mathbf{q}_1$ and stuff orthogonal to it. Throw away stuff
in the span of $\mathbf{q}_1$ since we have it covered already. Take
what’s left of $\mathbf{v}_2$ and turn it unto a unit vector.

``` r
y2 <- v2 - q1 %*% t(q1) %*% v2
r22 <- Norm(y2, 2)
q2 <- y2 / r22
```

Step 3. Think of $\mathbf{v}_3$ as made up of stuff in the subspace
spanned by $\mathbf{q}_{1,2}$ and stuff orthogonal to it. Through away
stuff in the span of $\mathbf{q}_{1,2}$ since we have it covered
already. Take what’s left of $\mathbf{v}_3$ and turn it unto a unit
vector.

``` r
y3 <- v3 - q1 %*% t(q1) %*% v3 - q2 %*% t(q2) %*% v3
r33 <- Norm(y3, 2)
```

Let’s check the result and see if the procedure worked.

``` r
q1
```

    ## [1] 0.333333333333 0.666666666667 0.666666666667

``` r
q2
```

    ##                 [,1]
    ## [1,]  0.666666666667
    ## [2,]  0.333333333333
    ## [3,] -0.666666666667

``` r
Norm(q1, 2)
```

    ## [1] 1

``` r
Norm(q2, 2)
```

    ## [1] 1

``` r
t(q1) %*% q2
```

    ##      [,1]
    ## [1,]    0

## QR decomposition

Just like **LU** decomposition is a way of using matrices to encode the
process of Gaussian elimination, **QR** decomposition is a way of using
matrices to encode the process of Gram-Schmidt orthogonalization. Let’s
take the equations we implemented above:

$$
\begin{align}
\mathbf{y}_1 &= \mathbf{v}_1\\
\mathbf{q}_1 &= \frac{1}{||\mathbf{y}_1||} \mathbf{y}_1 \\
\mathbf{y}_2 &= \mathbf{v}_2 - \mathbf{q}_1 \mathbf{q}_1^T \mathbf{v}_2 \\
\mathbf{q}_2 &= \frac{1}{||\mathbf{y}_2||} \mathbf{y}_2\\
\mathbf{y}_3 &= \mathbf{v}_3 - ( \mathbf{q}_1 \mathbf{q}_1^T \mathbf{v}_3 + \mathbf{q}_2 \mathbf{q}_2^T \mathbf{v}_3)\\
\mathbf{q}_3 &= \frac{1}{||\mathbf{y}_3||} \mathbf{y}_3.
\end{align}
$$

Let’s rewrite these equations, solving for the $\mathbf{v}_i$:

$$
\begin{align}
\mathbf{v}_1 &= ||\mathbf{y}_1|| \mathbf{q}_1 \\
\mathbf{v}_2 &= (\mathbf{q}_1\cdot \mathbf{v}_2) \mathbf{q}_1 + ||\mathbf{y}_2|| \mathbf{q}_2 \\
\mathbf{v}_3 &= (\mathbf{q}_1\cdot \mathbf{v}_3) \mathbf{q}_1 + (\mathbf{q}_2\cdot \mathbf{v}_3) \mathbf{q}_2  + ||\mathbf{y}_3|| \mathbf{q}_3 
\end{align}
$$

or, as matrices,

$$
\begin{pmatrix} | & | &   | \\ \mathbf{v}_1 & \mathbf{v}_2  & \mathbf{v}_3 \\| & | &   | \\ \end{pmatrix} = 
\underbrace{ \begin{pmatrix} | & | &   | \\ \mathbf{q}_1 & \mathbf{q}_2  & \mathbf{q}_3 \\| & | &   | \\ \end{pmatrix}}_Q   \underbrace{\begin{pmatrix} r_{1,1} & r_{1,2} &   r_{1,3} \\ 0 &r_{2,2}  & r_{2,3} \\0 & 0 &   r_{3,3} \\ \end{pmatrix} }_R
$$

with

$$
r_{i,i} = ||\mathbf{y}_i||, \qquad r_{i,j} = \mathbf{q}_i \cdot \mathbf{v}_j.
$$

The **QR** decomposition is a matrix decomposition that writes a
$m \times n$ matrix $\mathbf{A}$ as a product
$\mathbf{A} = \mathbf{Q}\,\mathbf{R}$ where:

- $\mathbf{Q}$ is an $m \times r$ matrix with orthonormal columns, where
  $r$ is the number of linearly independent columns of $\mathbf{A}$.
- $\mathbf{R}$ is an $r \times n$ matrix which is upper triangular if
  $r=n$, or the top portion of an upper triangular matrix if $r<n$.
- The columns of $\mathbf{Q}$ span the same space as the columns of
  $\mathbf{A}$.
- The matrix $\mathbf{R}$ gives the change of basis between the vectors
  in $\mathbf{Q}$ and the vectors in $\mathbf{A}$.
- The decompisition is unique up to some sign changes, so if we require
  $R_{ii}\geq 0$, it is unique.
- If the columns of $\mathbf{A}$ are independent, then $R_{ii}\neq 0$.
- On the other hand, if column $j$ of $\mathbf{A}$ can be written as a
  linear combination of columns to the left, then $R_{jj}=0$.

The **QR** decomposition we’ve done so far is actually called the
partial **QR** decomposition We distinguish this from the **full** or
**complete** **QR** decomposition. In the latter, we include vectors
than span parts of the space not spanned by $\mathbf{A}$ itself. Of
course, these contribute nothing to the matrix $\mathbf{Q}$, so it
results in a bunch of 0’s in $\mathbf{R}$. Below is the key picture to
understand. Here, $\mathbf{Q}$ and $\mathbf{R}$ are the (partial) **QR**
decomposition and $\overline{Q}$ and $\overline{R}$ are the full
version.

![](fullQR.png)

Let’s make this clear with numerous examples!

First example. Let’s start with a $3 \times 3$ matrix and perform
Gram-Schmidt orthogonalization to obtain the **QR** decomposition.

``` r
# Define matrix
v1 <- c(2,3,6)
v2 <- c(10,8,9)
v3 <- c(19,-3,-13)
A <- cbind(v1,v2,v3)

# Step 1
y1 <- v1
r11 <- Norm(y1,2)
q1 <- y1/r11

# Step2
y2 <- v2 - (q1%*%t(q1))%*%v2
r22 <- Norm(y2, 2)
q2 <- y2/r22
r12 <- t(q1) %*% v2

# Step 3
y3 <- v3 - (q1%*%t(q1))%*%v3 - (q2%*%t(q2))%*%v3
r33 <- Norm(y3,2)
q3 <- y3/r33
r13 <- t(q1) %*% v3
r23 <- t(q2) %*% v3

# Assemble into Qbar and Rbar
Qbar <- cbind(q1,q2,q3)
Rbar <- rbind(c(r11,r12,r13), c(0,r22,r23), c(0,0,r33))
Qbar
```

    ##                  q1                                
    ## [1,] 0.285714285714  0.857142857143  0.428571428571
    ## [2,] 0.428571428571  0.285714285714 -0.857142857143
    ## [3,] 0.857142857143 -0.428571428571  0.285714285714

``` r
Rbar
```

    ##      [,1] [,2] [,3]
    ## [1,]    7   14   -7
    ## [2,]    0    7   21
    ## [3,]    0    0    7

In this example, notice that $\overline{\mathbf{Q}}$ has three columns.
This reflects the fact that the columns of $\mathbf{A}$ are linearly
independent. Hence, it spans all of $\mathbb{R}^3$ and the **QR** and
full **QR** decompositions are the same.

Let’s check our answer by hand and also check it against `R`’s built-in
capabilities.

``` r
Qbar%*%Rbar - A
```

    ##      v1 v2 v3
    ## [1,]  0  0  0
    ## [2,]  0  0  0
    ## [3,]  0  0  0

``` r
QRbarcheck <- qr(A)
Qbarcheck <- qr.Q(QRbarcheck, complete=TRUE)
Rbarcheck <- qr.R(QRbarcheck, complete=TRUE)
Qbar
```

    ##                  q1                                
    ## [1,] 0.285714285714  0.857142857143  0.428571428571
    ## [2,] 0.428571428571  0.285714285714 -0.857142857143
    ## [3,] 0.857142857143 -0.428571428571  0.285714285714

``` r
Qbarcheck
```

    ##                 [,1]            [,2]            [,3]
    ## [1,] -0.285714285714 -0.857142857143 -0.428571428571
    ## [2,] -0.428571428571 -0.285714285714  0.857142857143
    ## [3,] -0.857142857143  0.428571428571 -0.285714285714

``` r
S <- diag(c(-1,-1,-1))
Qbarcheck <- Qbarcheck%*%S
Rbarcheck <- S%*%Rbarcheck
Qbarcheck - Qbar
```

    ##                      q1                                     
    ## [1,] -1.11022302463e-16  0.00000000000e+00 1.66533453694e-16
    ## [2,]  0.00000000000e+00 -1.11022302463e-16 3.33066907388e-16
    ## [3,]  0.00000000000e+00 -1.11022302463e-16 6.10622663544e-16

``` r
Rbarcheck - Rbar
```

    ##      v1               v2                v3
    ## [1,]  0 1.7763568394e-15  0.0000000000e+00
    ## [2,]  0 8.8817841970e-16  0.0000000000e+00
    ## [3,]  0 0.0000000000e+00 -1.7763568394e-15

Let’s try an example we looked at earlier.

``` r
# Define matrix
v1 <- c(1,2,2)
v2 <- c(2,1,-2)
v3 <- c(1,1,0)
A <- cbind(v1,v2,v3)

# Step 1
y1 <- v1
r11 <- Norm(y1, 2)
q1 <- y1 / r11

# Step2
y2 <- v2 - (q1 %*% t(q1)) %*% v2
r22 <- Norm(y2, 2)
q2 <- y2 / r22
r12 <- t(q1) %*% v2

# Step 3
y3 <- v3 - (q1 %*% t(q1)) %*% v3 - (q2 %*% t(q2)) %*% v3
y3
```

    ##                   [,1]
    ## [1,] 1.11022302463e-16
    ## [2,] 5.55111512313e-17
    ## [3,] 0.00000000000e+00

Oh! it turns out there’s nothing left. If we want to compute the full QR
decomposition, we need to find something orthogonal to the span of
$q_{1,2}$. An easy way to do that is to pick a vector outside of their
span and continue the orthogonalization procedure. You can
algorithmically use your linear algebra skills to do this, but for
little cases like ours, you can eyeball it.

``` r
v3 <- c(1,1,1)
y3 <- v3 - (q1 %*% t(q1)) %*% v3 - (q2 %*% t(q2)) %*% v3
q3 <- y3 / Norm(y3, 2)
r13 <- dot(q1,v3)
r23 <- dot(q2,v3)
r33 <- 0

# Assemble into Qbar and Rbar
Qbar <- cbind(q1,q2,q3)
Rbar <- rbind(c(r11,r12,r13), c(0,r22,r23), c(0,0,r33))
Qbar
```

    ##                  q1                                
    ## [1,] 0.333333333333  0.666666666667  0.666666666667
    ## [2,] 0.666666666667  0.333333333333 -0.666666666667
    ## [3,] 0.666666666667 -0.666666666667  0.333333333333

``` r
Rbar
```

    ##      [,1] [,2]           [,3]
    ## [1,]    3    0 1.666666666667
    ## [2,]    0    3 0.333333333333
    ## [3,]    0    0 0.000000000000

To reiterate, we ended up with a row of zeros at the bottom of
$\overline{\mathbf{R}}$. That’s because the columns of $\mathbf{A}$ are
linearly dependent and don’t span $\mathbb{R}^3$. At any rate, let’s go
ahead and check our result.

``` r
Qbar %*% Rbar - A
```

    ##      v1 v2              v3
    ## [1,]  0  0 -0.222222222222
    ## [2,]  0  0  0.222222222222
    ## [3,]  0  0  0.888888888889

In this case, the full and reduced **QR** are not the same. Let’s see
how this works.

``` r
Q <- Qbar[1:3,1:2]
R <- Rbar[1:2,1:3]
Q
```

    ##                  q1                
    ## [1,] 0.333333333333  0.666666666667
    ## [2,] 0.666666666667  0.333333333333
    ## [3,] 0.666666666667 -0.666666666667

``` r
R
```

    ##      [,1] [,2]           [,3]
    ## [1,]    3    0 1.666666666667
    ## [2,]    0    3 0.333333333333

``` r
Q %*% R - A
```

    ##      v1 v2              v3
    ## [1,]  0  0 -0.222222222222
    ## [2,]  0  0  0.222222222222
    ## [3,]  0  0  0.888888888889

Ok, and let’s do one last example.

``` r
# Define matrix
v1 <- c(1,1,1,1)
v2 <- c(0,1,1,1)
v3 <- c(0,0,1,1)
A <- cbind(v1,v2,v3)
A
```

    ##      v1 v2 v3
    ## [1,]  1  0  0
    ## [2,]  1  1  0
    ## [3,]  1  1  1
    ## [4,]  1  1  1

``` r
# Step 1
y1 <- v1
r11 <- Norm(y1, 2)
q1 <- y1 / r11

# Step 2
y2 <- v2 - (q1 %*% t(q1)) %*% v2
r22 <- Norm(y2, 2)
q2 <- y2 / r22
r12 <- t(q1) %*% v2

# Step 3
y3 <- v3 - (q1 %*% t(q1)) %*% v3 - (q2 %*% t(q2)) %*% v3
r33 <- Norm(y3, 2)
q3 <- y3 / (r33)
r13 <- t(q1) %*% v3
r23 <- t(q2) %*% v3
```

But wait! We are living in $\mathbb{R}^4$ and we only have three vectors
to far, $\mathbf{q}_{1,2,3}$. If we want the full decomposition, we have
to find a basis for the orthogonal complement of $\mathbf{A}$.

``` r
# Choose a vector not in the span of q1, q2, q3
v4 <- c(1,2,3,4)
y4 <- v4 - (q1 %*% t(q1)) %*% v4 - (q2 %*% t(q2)) %*% v4  - (q3 %*% t(q3)) %*% v4
r44 <- Norm(y4, 2) 
q4 <- y4 / (r44)
r14 <- t(q1) %*% v4
r24 <- t(q2) %*% v4
r34 <- t(q3) %*% v4

# Assemble into Qbar and Rbar, check answer
Qbar <- cbind(q1,q2,q3,q4)
Rbar <- rbind(c(r11,r12,r13), c(0,r22,r23), c(0,0,r33), c(0,0,0))
Qbar
```

    ##       q1                                                      
    ## [1,] 0.5 -0.866025403784  1.35973995551e-16  3.92523114671e-16
    ## [2,] 0.5  0.288675134595 -8.16496580928e-01 -3.14018491737e-16
    ## [3,] 0.5  0.288675134595  4.08248290464e-01 -7.07106781187e-01
    ## [4,] 0.5  0.288675134595  4.08248290464e-01  7.07106781187e-01

``` r
Rbar
```

    ##      [,1]           [,2]           [,3]
    ## [1,]    2 1.500000000000 1.000000000000
    ## [2,]    0 0.866025403784 0.577350269190
    ## [3,]    0 0.000000000000 0.816496580928
    ## [4,]    0 0.000000000000 0.000000000000

``` r
Qbar %*% Rbar - A
```

    ##      v1 v2                 v3
    ## [1,]  0  0 -1.23259516441e-32
    ## [2,]  0  0  0.00000000000e+00
    ## [3,]  0  0  0.00000000000e+00
    ## [4,]  0  0  0.00000000000e+00

``` r
# Make reduced QR and check
Q <- Qbar[1:4,1:3]
R <- Rbar[1:3,1:3]
Q %*% R - A
```

    ##      v1 v2                 v3
    ## [1,]  0  0 -1.23259516441e-32
    ## [2,]  0  0  0.00000000000e+00
    ## [3,]  0  0  0.00000000000e+00
    ## [4,]  0  0  0.00000000000e+00

## Computational considerations

So far, we have calculated the **QR** decomposition using a methodology
that is based on Gram-Schmidt orthogonalization. This is because
Gram-Schmidt is the most conceptually straightforward method. For an
$m \times n$ matrix, it turns out that the computational cost of
factorization is $\mathcal{O}(2mn^2)$ multiplications and additions. An
alternative choice for factorization uses objects called Householder
reflections. I won’t go in to these here, but the method requires fewer
operations and is the one actually implemented in many software
packages.

## What is the point of QR factorization

Finally, we ask: why have we bothered to do all of this? It is really,
really convenient for least squares, and turns out to have very nice
numerical properties because of small condition numbers. You’ll work
with the numerical issue on your activities and/or homework, but for
now, here’s how least squares works when you use **QR** decomposition on
$\mathbf{A}\,\mathbf{x}=\mathbf{b}$. Note

$$
\begin{align}
||\mathbf{A}\mathbf{x}-\mathbf{b}||^2 &= ||\bar{\mathbf{Q}}\bar{\mathbf{R}} \mathbf{x} - \mathbf{b}||^2 \\
&= ||\bar{\mathbf{Q}}^T(\bar{\mathbf{Q}}\bar{\mathbf{R}} \mathbf{x} - \mathbf{b})||^2 \\
&= ||\bar{\mathbf{Q}}^T \bar{\mathbf{Q}}\bar{\mathbf{R}} \mathbf{x} - \bar{\mathbf{Q}}^T \mathbf{b}||^2 \\
&= ||\bar{\mathbf{R}} \mathbf{x} - \bar{\mathbf{Q}}^T \mathbf{b}||^2 \\
& = ||\mathbf{R} \mathbf{x}-\mathbf{Q}^T \mathbf{b}||^2+ ||\widehat{\mathbf{Q}}^T \mathbf{b}||^2. 
\end{align}
$$

Proceeding from the second to last line to the very last line is not
obvious, and takes you a little bit of writing out expressions more
explicitly. This is a good exercise for you to try.

In any case, to minimize this, note that the second term doesn’t even
include $\mathbf{x}$ so there’s nothing we can do about it (in fact, it
is the norm of the residual). To minimize the quantity, we can then just
force the first term to be zero. This is fine because it is a square
system! There’s a command that does this all automatically for you
called `qr.solve`.

# Eigenvalues I

## Big picture

Along with solving the linear system
$\mathbf{A}\,\mathbf{x}=\mathbf{b}$, finding the eigenvalues and
eigenvectors of a matrix $\mathbf{A}$ is one of the most important
problems in linear algebra. Knowing the eigenpairs can help simplify a
problem and reveal important information about systems modeled with
linear algebra.

## Goals

- Define eigenvalues and eigenvectors, and calculate them by hand for
  small matrices
- Diagonalize matrices and state the algebraic and geometric
  multiplicity of eigenvalues
- Apply the power iteration technique to find an eigenvalue

## Eigenvalue fundamentals

For an $n \times n$ matrix $\mathbf{A}$, a scalar
$\lambda \in \mathbb{C}$, and vector $\mathbf{v} \in \mathbb{R}^n$,
$\mathbf{v} \neq \mathbf{0}$, then we say $\lambda$ is an **eigenvalue**
of $\mathbf{A}$ and $\mathbf{v}$ is an **eigenvector** of $\mathbf{A}$
if $\mathbf{A}\,\mathbf{v}=\lambda \mathbf{v}$.

Stated in words: an eigenvector and eigenvalue are the magical vector
$\mathbf{v}$ and scalar $\lambda$ such that if you hit $\mathbf{v}$ with
$\mathbf{A}$, you get back the same vector $\mathbf{v}$ but multiplied
by a constant $\lambda$.

How do we calculate them? Let’s take the definition
$\mathbf{A} \,\mathbf{v} = \lambda \mathbf{v}$ and rearrange it to write
$(\mathbf{A} - \mathbf{I} \lambda) \mathbf{v} = \mathbf{0}$. There are
only two ways this can happen. One choice is $\mathbf{v}=0$, but that’s
trivial because it works for any $\mathbf{A}$. The other choice, by the
Invertible Matrix Theorem, is that the matrix on the left is singular.
Also from the Invertible Matrix Theorem, if it is singular, then it has
determinant zero, that is $\det (\mathbf{A}- \mathbf{I} \lambda) = 0$.
This equation is a polynomial in $\lambda$ and is called the
**characteristic polynomial**. When calculating by hand, we find the
characteristic polynomial first and then solve it to find the
eigenvalues. To find eigenvectors, we remember that
$\mathbf{A}\,\mathbf{v} = \lambda \mathbf{v} \rightarrow (\mathbf{A}-\mathbf{I}\lambda)\mathbf{v}=\mathbf{0}$
and solve for $\mathbf{v}$.

To recap:

1.  Solve the characteristic equation
    $\det (\mathbf{A}- \mathbf{I} \lambda) = 0$ to find the $\lambda_i$.
2.  Solve $(\mathbf{A}-\mathbf{I}\lambda_i)\mathbf{v_i}=\mathbf{0}$ to
    find the $\mathbf{v}_i$.

For example, let’s find the eigenvalues and eigenvectors of

$$
\mathbf{A} = \begin{pmatrix} -3 & 2 \\ 2 & -3 \end{pmatrix}.
$$

Using the result above, we can write the characteristic polynomial:

$$
\begin{align}
\det (\mathbf{A} - \mathbf{I} \lambda) &= 0\\
\det \begin{pmatrix} -3 - \lambda & 2 \\ 2 & -3 - \lambda \end{pmatrix} &= 0\\
(\lambda+3)^2-4 &= 0\\
\lambda^2 + 6\lambda + 5 &= 0\\
(\lambda+5)(\lambda+1) &= 0
\end{align}
$$

and therefore $\lambda_{1,2} = -1, -5$. To find $\mathbf{v}_1$, we
solve:

$$
\begin{align}
(\mathbf{A}-\mathbf{I}\lambda_1)\mathbf{v}_1&=\mathbf{0}\\
\begin{pmatrix} -2 & 2 \\ 2 & -2 \end{pmatrix}\mathbf{v}_1 &= \mathbf{0} \\
\mathbf{v}_1 &= \begin{pmatrix} 1 \\ 1 \end{pmatrix}
\end{align}
$$

or any scalar multiple of this vector. Similarly, we find

$$
\mathbf{v}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}.
$$

## Algebraic multiplicity, geometric multiplicity, and diagonalization

Many applications of eigenvalues are intimiately tied up with the idea
of **diagonalization** of matrices. Suppose $\mathbf{A}$ has eigenpairs
$\lambda_i$, $\mathbf{v}_i$, $i = 1,\ldots,n$. Then we can write down
the definition of eigenpair for all pairs simultaneously:
$\mathbf{A}\,\mathbf{v}_i = \lambda_i \mathbf{v}_i$ implies

$$
\begin{align}
\mathbf{A} \underbrace{\begin{pmatrix}
\vert & \vert & \cdots & \vert \\
\mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n \\
\vert & \vert & \cdots & \vert 
\end{pmatrix}}_\mathbf{S} &= 
{\begin{pmatrix} \vert & \vert & \cdots & \vert \\
\lambda_1 \mathbf{v}_1 & \lambda_2 \mathbf{v}_2 & \cdots & \lambda_n \mathbf{v}_n \\
\vert & \vert & \cdots & \vert 
\end{pmatrix}} \\
&=\underbrace{\begin{pmatrix} \vert & \vert & \cdots & \vert \\
\mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n \\
\vert & \vert & \cdots & \vert 
\end{pmatrix}}_\mathbf{S}
\underbrace{\begin{pmatrix} \lambda_1 &&& \\ & \lambda_2 & & \\ && \ddots & \\ &&& \lambda_n   \end{pmatrix}}_\mathbf{\Lambda}.
\end{align}
$$

Since $\mathbf{A}\,\mathbf{S} = \mathbf{S}\, \mathbf{\Lambda}$, we can
write $\mathbf{A} = \mathbf{S} \mathbf{\Lambda} \mathbf{S}^{-1}$. If we
think of $\mathbf{S}$ as describing a change of basis, this equation
says that the action of $\mathbf{A}$ is like going into another basis,
multiplying by a diagonal matrix, and then changing back to the original
basis.

Let’s show that we can do this with our example matrix $\mathbf{A}$ from
before,

$$
\mathbf{A} = \begin{pmatrix} -3 & 2 \\ 2 & -3 \end{pmatrix}.
$$

``` r
A <- matrix(c(-3,2,2,-3), byrow = TRUE, nrow = 2)
e <- eigen(A)
e
```

    ## eigen() decomposition
    ## $values
    ## [1] -1 -5
    ## 
    ## $vectors
    ##                [,1]            [,2]
    ## [1,] 0.707106781187  0.707106781187
    ## [2,] 0.707106781187 -0.707106781187

``` r
Lambda <- diag(e$values)
S <- e$vectors
S %*% Lambda %*% solve(S)
```

    ##      [,1] [,2]
    ## [1,]   -3    2
    ## [2,]    2   -3

Why would you ever need/want to diagonalize a matrix? There are many
reasons related to data analysis, differential equations, graph theory,
and more. During the next class you will see some of these applications
of eigenvalues and eigenvectors.

We’ll do one classic example: the Fibonacci sequence. First though,
let’s do a quick warm-up problem. Suppose we have a sequence defined as

$$
A_0 = 2,\quad A_n = 5A_{n-1}.
$$

What is a formula for $A_n$? Well, notice that

$$
\begin{align}
A_1 &= 2 \cdot 5\\
A_2 &= 2 \cdot 5 \cdot 5\\
&\vdots\\
A_n &= 2 \cdot 5^n.
\end{align}
$$

We’ll need this result! Now recall the Fibonacci sequence,

$$
F_0 = 1,\quad F_1 = 1, \quad F_n = F_{n-1} + F_{n-2}, \quad n \geq 2.
$$

If I asked you for $F_{1000}$, then using the recursive definition
above, you’d need to compute the $999$ terms that come before it. But is
there a way to avoid doing this? That is, is there a way to come up with
an explicit (non-recursive) formula for $F_n$?

So that we can use tools of linear algebra, let’s transform this to a
problem involving a matrix. Define $G_n = F_{n-1}$. Then substituting
into the definition of $F_n$, we can write $F_n = F_{n-1} + G_{n-1}$.
Putting it all together,

$$
\begin{align}
F_n &= F_{n-1} + G_{n-1} \\
G_n &= F_{n-1} \\
\begin{pmatrix} F_n \\ G_n \end{pmatrix} &= \begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix} \begin{pmatrix} F_{n-1} \\ G_{n-1} \end{pmatrix} \\
\mathbf{F}_n &= \mathbf{A} \mathbf{F}_{n-1}.
\end{align}
$$

Let’s diagonalize $\mathbf{A}$, and momentarily, you’ll see why.

``` r
A <- matrix(c(1,1,1,0), byrow = TRUE, nrow = 2)
lambdap <- (1+sqrt(5)) / 2
lambdap
```

    ## [1] 1.61803398875

``` r
lambdam <- (1-sqrt(5)) / 2
lambdam
```

    ## [1] -0.61803398875

``` r
e <- eigen(A)
e$values
```

    ## [1]  1.61803398875 -0.61803398875

``` r
Lambda <- diag(e$values)
Lambda
```

    ##               [,1]           [,2]
    ## [1,] 1.61803398875  0.00000000000
    ## [2,] 0.00000000000 -0.61803398875

``` r
S <- e$vectors
v1 <- S[,1] / S[2,1]
v2 <- S[,2] / S[2,2]
S <- cbind(v1, v2)
S
```

    ##                 v1             v2
    ## [1,] 1.61803398875 -0.61803398875
    ## [2,] 1.00000000000  1.00000000000

``` r
S %*% Lambda %*% solve(S)
```

    ##      [,1]               [,2]
    ## [1,]    1  1.00000000000e+00
    ## [2,]    1 -1.11022302463e-16

So to recap, we have

$$
\mathbf{\Lambda} = \begin{pmatrix} \lambda_+ & 0 \\ 0 & \lambda_- \end{pmatrix}, \quad \mathbf{S} = \begin{pmatrix} \lambda_+ & \lambda_- \\ 1 & 1 \end{pmatrix}.
$$

Why would we do this? Let’s use the diagonalization to re-write our
problem:

$$
\begin{align}
\mathbf{F}_n &= \mathbf{A} \mathbf{F}_{n-1}\\
\mathbf{F}_n &= \mathbf{S} \mathbf{\Lambda} \mathbf{S}^{-1} \mathbf{F}_{n-1}\\
\mathbf{S}^{-1} \mathbf{F} &= \mathbf{\Lambda} \mathbf{S}^{-1} \mathbf{F}_{n-1}.
\end{align}
$$

Now define

$$
\mathbf{S}^{-1} \mathbf{F} = \mathbf{S}^{-1} \begin{pmatrix} F_n \\ G_n \end{pmatrix} \equiv  \begin{pmatrix} A_n \\ B_n \end{pmatrix}.
$$

Then our problem becomes

$$
\begin{align}
\begin{pmatrix} A_n \\ B_n \end{pmatrix} &= \mathbf{\Lambda} \begin{pmatrix} A_{n-1} \\ B_{n-1} \end{pmatrix} \\
\begin{pmatrix} A_n \\ B_n \end{pmatrix} &= \begin{pmatrix} \lambda_+ &  0 \\ 0 & \lambda_- \end{pmatrix} \begin{pmatrix} A_{n-1} \\ B_{n-1}. \end{pmatrix}
\end{align}
$$

Because this problem is diagonal, the top and bottom equations are
uncoupled, which makes them much easier to solve. In fact, they are
geometrtic, just like the warm-up problem we did. We can write

$$
\begin{pmatrix} A_n \\ B_n \end{pmatrix} = \begin{pmatrix}  A_0 \lambda_+^n \\ B_0 \lambda_-^n \end{pmatrix}.
$$

But we didn’t want to know $A_n$ and $B_n$. We wanted to know $F_n$ and
$G_n$. So we have to undo the transformation we did before. We have

$$
\begin{align}
\mathbf{S}^{-1} \mathbf{F} &= \begin{pmatrix} A_n \\ B_n \end{pmatrix} \\
\mathbf{S}^{-1} \mathbf{F} &= \begin{pmatrix}  A_0 \lambda_1^n \\ B_0 \lambda_2^n \end{pmatrix} \\
\mathbf{F} &= \mathbf{S} \begin{pmatrix}  A_0 \lambda_1^n \\ B_0 \lambda_2^n \end{pmatrix} \\
\mathbf{F} &= \begin{pmatrix} \lambda_1 & \lambda_2 \\ 1 & 1 \end{pmatrix} \begin{pmatrix}  A_0 \lambda_1^n \\ B_0 \lambda_2^n \end{pmatrix} \\
\mathbf{F} &= \begin{pmatrix} \lambda_1 A_0 \lambda_1^n + \lambda_2 B_0 \lambda^2n \\ \mathrm{doesn't\ matter}\end{pmatrix}.
\end{align}
$$

Therefore,

$$
F_n = A_0 \lambda_+^{n+1} + B_0 \lambda_-^{n+1}.
$$

To find the unknown constants, we plug in the initial conditions.

$$
\begin{align}
F_0 &= A_0 \lambda_+ + B_0 \lambda_- = 1\\
F_1 &= A_0 \lambda_+^2 + B_0 \lambda_-^2 = 1.
\end{align}
$$

Solving, we find $A_0 = 1/\sqrt{5}$, $B_0 = -1/\sqrt{5}$. The final
answer, then, is

$$
F_n = \frac{1}{\sqrt{5}}\lambda_+^{n+1} - \frac{1}{\sqrt{5}}\lambda_-^{n+1},\quad \lambda_\pm = \frac{1 \pm \sqrt{5}}{2}.
$$

Let’s test this out.

``` r
fib <- function(n){
  lambdap <- (1 + sqrt(5)) / 2
  lambdam <- (1 - sqrt(5)) / 2
  Fn <- 1/sqrt(5) * lambdap^(n+1) - 1/sqrt(5) * lambdam^(n+1)
  return(Fn)
}
fib(1:10)
```

    ##  [1]  1  2  3  5  8 13 21 34 55 89

``` r
fib(1000)
```

    ## [1] 7.03303677114e+208

One cool thing about this is that it lets us understand the behavior of
$F_n$ for large $n$ in a simpler way. Since $|\lambda_-| < 1$, after
many repeated iterations, the term involving $\lambda_-$ will die out
and we can write, for large $n$,

$$
F_n \approx \frac{1}{\sqrt{5}}\lambda_+^{n+1}.
$$

Let’s test this out!

``` r
fibapprox <- function(n){
  lambdap <- (1 + sqrt(5)) / 2
  Fn <- 1/sqrt(5) * lambdap^(n+1)
  return(Fn)
}
n <- 1:10
fib(n)
```

    ##  [1]  1  2  3  5  8 13 21 34 55 89

``` r
fibapprox(n)
```

    ##  [1]  1.17082039325  1.89442719100  3.06524758425  4.95967477525  8.02492235950
    ##  [6] 12.98459713475 21.00951949425 33.99411662900 55.00363612325 88.99775275225

``` r
n <- 1:20
error <- abs(fib(n) - fibapprox(n))
plot(n, log10(error))
```

![](coursebook_files/figure-gfm/unnamed-chunk-65-1.png)<!-- -->

It’s important to remember from linear algebra that not every matrix can
be diagonalized. For a matrix to be diagonalizable, you need for the
multiplicity of each eigenvalue to be the same as the dimension of the
eigenspace. Here are a couple of examples.

``` r
A <- matrix(c(5,-4,4,12,-11,12,4,-4,5), byrow = TRUE, nrow = 3)
e <- eigen(A)
e$values
```

    ## [1] -3  1  1

``` r
e$vectors
```

    ##                 [,1]           [,2]            [,3]
    ## [1,] -0.301511344578 0.534522483825 -0.588533689677
    ## [2,] -0.904534033733 0.801783725737 -0.784390372213
    ## [3,] -0.301511344578 0.267261241912 -0.195856682537

Here, the **algebraic multiplicity** of $-3$ is 1 because it is only an
eigenvalue once, and the **geometric multiplicity** is 1 because it only
has one eigenvector. The algebraic multiplicity of $1$ is $2$ because it
is an eigenvalue twice. Since it has two independent eigenvectors, the
geometric multiplicity is also 2. Therefore, this matrix is
diagonalizable. Stated differently: we need enough independent
eigenvectors to form the $\mathbf{S}$ matrix.

By way of counterexample, consider this problem.

``` r
A <- matrix(c(1,-1,0,1), byrow = TRUE, nrow = 2)
A
```

    ##      [,1] [,2]
    ## [1,]    1   -1
    ## [2,]    0    1

``` r
e <- eigen(A)
e$values
```

    ## [1] 1 1

``` r
e$vectors
```

    ##      [,1]              [,2]
    ## [1,]    1 1.00000000000e+00
    ## [2,]    0 2.22044604925e-16

The eigenvalue $1$ has algebraic multiplicity 2, but there is only one
eigenvector (the trivial eigenvector doesn’t count) so it has geometric
multiplicity 1. We don’t have enough eigenvectors to make $\mathbf{S}$,
so the matrix is not diagonalizable.

## Power iteration

Please find the eigenvalues of this matrix:

``` r
set.seed(123)
A <- matrix(sample(-100:100,64), nrow = 8)
print(A)
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
    ## [1,]   58   98   36  -23    8   62   -4   15
    ## [2,]   78   97   -2  -20   80   77   87   -7
    ## [3,]  -87   52  -29   93   99  -67  -63  -95
    ## [4,]   94  -11  -75    2  -27  -32  -80  -15
    ## [5,]   69  -10  -94   16  -78   82  -60   41
    ## [6,]  -51   88   96  -25   54   71   89  -62
    ## [7,]   17   84   86   42  -48  -38  -41   91
    ## [8,]  -58   -9   63  -69   34   40  -85   95

Haha, you (probably) can’t! That’s because writing down the
characteristic polynomial is messy, and also, because for any polynomial
of degree greater than 4, we don’t have a formula for the solution, so
you’d have to resort to numerical methods to find the roots.

Instead, we will develop numerical methods that can find an eigenvalue
more directly. Let’s do a numerical experiment.

``` r
A <- matrix(c(-13,170,240,19,-224,-320,-14,166,237), byrow = TRUE, nrow = 3)
A
```

    ##      [,1] [,2] [,3]
    ## [1,]  -13  170  240
    ## [2,]   19 -224 -320
    ## [3,]  -14  166  237

``` r
v <- c(1,1,1)
v <- A %*% v
v
```

    ##      [,1]
    ## [1,]  397
    ## [2,] -525
    ## [3,]  389

``` r
v <- v / Norm(v, 2)
v
```

    ##                 [,1]
    ## [1,]  0.519251568355
    ## [2,] -0.686667691150
    ## [3,]  0.508788060681

``` r
v <- A %*% v
v
```

    ##                 [,1]
    ## [1,] -1.374643320759
    ## [2,]  0.867163198538
    ## [3,] -0.673588306557

``` r
v <- v / Norm(v, 2)
for (i in 1:100){
  v <- A %*% v
  v <- v / Norm(v, 2)
}
(A %*%v )/ v
```

    ##      [,1]
    ## [1,]   -3
    ## [2,]   -3
    ## [3,]   -3

``` r
e <- eigen(A)
e$values
```

    ## [1] -3  2  1

Interesting. Using this iteration, we found the eigenvalue $-3$. Any
guesses why?

Let’s try another example.

``` r
A <- matrix(c(2,-520,8,-90,468,-360,40,-698,160), byrow = TRUE, nrow = 3)
v <- c(1,1,1)
for (i in 1:100){
  v <- A %*% v
  v <- v / Norm(v, 2)
}
(A %*% v)/v
```

    ##      [,1]
    ## [1,]  882
    ## [2,]  882
    ## [3,]  882

``` r
e <- eigen(A)
e$values
```

    ## [1]  8.82000000000e+02 -2.52000000000e+02 -4.63379334903e-14

This method is called **power iteration**. You might not realize it, but
you use power iteration nearly every day. Google searches are based on
an algorithm called PageRank, which is applies power iteration to a
matrix that encodes links between web pages. The matrix in question is
gigantic, which is why numerical methods are necessary.
