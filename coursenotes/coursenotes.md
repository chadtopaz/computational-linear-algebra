Computational Linear Algebra Course Notes
================
Chad M. Topaz

-   <a href="#r-bootcamp" id="toc-r-bootcamp">R Bootcamp</a>
-   <a href="#how-computers-store-numbers"
    id="toc-how-computers-store-numbers">How Computers Store Numbers</a>
-   <a href="#fundamentals-of-linear-systems"
    id="toc-fundamentals-of-linear-systems">Fundamentals of Linear
    Systems</a>
-   <a href="#lu-decomposition" id="toc-lu-decomposition">LU
    decomposition</a>
-   <a href="#iterative-methods-for-linear-systems"
    id="toc-iterative-methods-for-linear-systems">Iterative methods for
    linear systems</a>
-   <a href="#polynomial-interpolation"
    id="toc-polynomial-interpolation">Polynomial interpolation</a>
-   <a href="#interpolation-error-and-chebyshev-interpolation"
    id="toc-interpolation-error-and-chebyshev-interpolation">Interpolation
    error and Chebyshev interpolation</a>
-   <a href="#splines" id="toc-splines">Splines</a>
-   <a href="#least-squares" id="toc-least-squares">Least squares</a>
-   <a href="#qr-factorization" id="toc-qr-factorization">QR
    Factorization</a>
-   <a href="#eigenvalues" id="toc-eigenvalues">Eigenvalues</a>

# R Bootcamp

## Big Picture

It’s time to review some math and start getting experience with R
programming.

## Goals

-   Implement Horner’s method
-   Explain advantage of Horner’s method
-   Define different ways to multiply vectors
-   Compute Taylor polynomials and bound the error

## Horner’s Method

Evaluating a polynomial might seem boring, and admittedly, it sort of
is. But even a simple task like evaluating a polynomial can involve
complexities. Suppose that evaluating a polynomial on your comptuer
takes 0.005 seconds. That seems fast, but what if you are solving a
heavy-duty industrial problem that requires millions of evaluations?

Consider some methods for evaluating the polynomial
![2x^4+3x^3-3x^2+5x-1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2x%5E4%2B3x%5E3-3x%5E2%2B5x-1 "2x^4+3x^3-3x^2+5x-1"):

1.  Naively,
    ![2 \cdot x \cdot x \cdot x \cdot x + 3 \cdot x \cdot x \cdot x \ldots](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2%20%5Ccdot%20x%20%5Ccdot%20x%20%5Ccdot%20x%20%5Ccdot%20x%20%2B%203%20%5Ccdot%20x%20%5Ccdot%20x%20%5Ccdot%20x%20%5Cldots "2 \cdot x \cdot x \cdot x \cdot x + 3 \cdot x \cdot x \cdot x \ldots").

-   10 multiplications, 4 additions (includes subtractions)

2.  Store
    ![(1/2)(1/2)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%281%2F2%29%281%2F2%29 "(1/2)(1/2)"),
    bootstrap to
    ![(1/2)(1/2)(1/2)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%281%2F2%29%281%2F2%29%281%2F2%29 "(1/2)(1/2)(1/2)"),
    and so forth.

-   7 multiplications, 4 additions

3.  Horner’s method, factoring out successive powers of
    ![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x"),
    as in
    ![-1+x(5+x(-3+x(3+2x)))](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-1%2Bx%285%2Bx%28-3%2Bx%283%2B2x%29%29%29 "-1+x(5+x(-3+x(3+2x)))").

-   4 multiplications, 4 additions

We can test:

``` r
xvec <- runif(10^6)
t1 <- system.time(
  for (x in xvec){
    2*x^4 + 3*x^3 - 3*x^2 + 5*x - 1
  }
)[3]
t2 <- system.time(
  for (x in xvec){
    -1 + x*(5 + x*(-3 + x*(3 + 2*x)))
  }
)[3]
t1/t2
```

    ##  elapsed 
    ## 2.195652

## Inner and Outer Products

There are several different ways to “multiply” vectors
![\textbf{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D "\textbf{x}")
and
![\textbf{y}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7By%7D "\textbf{y}"):

1.  Element-wise. Just multiply each element of
    ![\textbf{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D "\textbf{x}")
    with the element in the corresponding poisiotn in
    ![\textbf{y}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7By%7D "\textbf{y}").
    Note that
    ![\textbf{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D "\textbf{x}")
    and
    ![\textbf{y}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7By%7D "\textbf{y}")
    must have the same number of elements, and the result is a vector
    with the same number of elements.

2.  Dot product, also known as the inner product. To compute, calculate
    ![\textbf{x} \cdot \textbf{y} \equiv \textbf{x}^T \textbf{y}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D%20%5Ccdot%20%5Ctextbf%7By%7D%20%5Cequiv%20%5Ctextbf%7Bx%7D%5ET%20%5Ctextbf%7By%7D "\textbf{x} \cdot \textbf{y} \equiv \textbf{x}^T \textbf{y}")
    where the superscript
    ![T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T "T")
    indicates transpose. Dot product is related to the angle
    ![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta")
    between the two vectors as
    ![\textbf{x} \cdot \textbf{y} = \|\textbf{x}\|\|\textbf{y}\| \cos \theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D%20%5Ccdot%20%5Ctextbf%7By%7D%20%3D%20%7C%5Ctextbf%7Bx%7D%7C%7C%5Ctextbf%7By%7D%7C%20%5Ccos%20%5Ctheta "\textbf{x} \cdot \textbf{y} = |\textbf{x}||\textbf{y}| \cos \theta").
    Rearranging this as
    ![\textbf{x} \cdot \textbf{y}/\|\textbf{y}\| = \|\textbf{x}\| \cos \theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D%20%5Ccdot%20%5Ctextbf%7By%7D%2F%7C%5Ctextbf%7By%7D%7C%20%3D%20%7C%5Ctextbf%7Bx%7D%7C%20%5Ccos%20%5Ctheta "\textbf{x} \cdot \textbf{y}/|\textbf{y}| = |\textbf{x}| \cos \theta")
    suggests the intuition of the dot product. It calculates a
    *projection* of one vector onto the other, or restated, it tells us
    how much of one vector is pointing in the direction of the other
    vector. Note that
    ![\textbf{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D "\textbf{x}")
    and
    ![\textbf{y}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7By%7D "\textbf{y}")
    must have the same number of elements, and the result is a scalar.

3.  Outer product. To compute, calculate
    ![\textbf{x} \otimes \textbf{y} \equiv \textbf{x}\textbf{y}^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D%20%5Cotimes%20%5Ctextbf%7By%7D%20%5Cequiv%20%5Ctextbf%7Bx%7D%5Ctextbf%7By%7D%5ET "\textbf{x} \otimes \textbf{y} \equiv \textbf{x}\textbf{y}^T").
    Note that
    ![\textbf{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D "\textbf{x}")
    and
    ![\textbf{y}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7By%7D "\textbf{y}")
    do not need to have the same number of elements, and the result is a
    matrix. If
    ![\textbf{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D "\textbf{x}")
    is
    ![m \times 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m%20%5Ctimes%201 "m \times 1")
    and
    ![\textbf{y}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7By%7D "\textbf{y}")
    is
    ![n \times 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%201 "n \times 1")
    then
    ![\textbf{x} \otimes \textbf{y}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D%20%5Cotimes%20%5Ctextbf%7By%7D "\textbf{x} \otimes \textbf{y}")
    is
    ![m \times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m%20%5Ctimes%20n "m \times n").

4.  Cross product. We are not going to worry about this one for now.

Let’s calculate some examples.

``` r
x <- c(1,2,3)
y <- c(4,5,6)
z <- c(7,8,9,10)
x*y # element-wise
```

    ## [1]  4 10 18

``` r
t(x) %*% y # dot product
```

    ##      [,1]
    ## [1,]   32

``` r
sum(x*y) # dot product
```

    ## [1] 32

``` r
x %*% t(z) # outer product
```

    ##      [,1] [,2] [,3] [,4]
    ## [1,]    7    8    9   10
    ## [2,]   14   16   18   20
    ## [3,]   21   24   27   30

``` r
x %o% z # outer product
```

    ##      [,1] [,2] [,3] [,4]
    ## [1,]    7    8    9   10
    ## [2,]   14   16   18   20
    ## [3,]   21   24   27   30

## Taylor’s Theorem

This theorem will be especially useful for error analysis of some
algorithms we use. The basic idea of Taylor’s theorem is that for many
functions, we can approximate the function around a point
![x_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_0 "x_0")
as a polynomial of finite degree, called the Taylor polynomial, plus an
error term that accounts for all the higher degree terms we ignored.

More formally, suppose
![f](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f "f")
is
![n+1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%2B1 "n+1")
times continuously differentiablw on the interval between
![x_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_0 "x_0")
and
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x").
Then

![f(x) = \left( \displaystyle \sum\_{k=0}^n \frac{f^{(k)}(x_0)}{k!}(x-x_0)^k \right) + \frac{f^{(n+1)}(c)}{(n+1)!}(x-x_0)^{n+1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28x%29%20%3D%20%5Cleft%28%20%5Cdisplaystyle%20%5Csum_%7Bk%3D0%7D%5En%20%5Cfrac%7Bf%5E%7B%28k%29%7D%28x_0%29%7D%7Bk%21%7D%28x-x_0%29%5Ek%20%5Cright%29%20%2B%20%5Cfrac%7Bf%5E%7B%28n%2B1%29%7D%28c%29%7D%7B%28n%2B1%29%21%7D%28x-x_0%29%5E%7Bn%2B1%7D "f(x) = \left( \displaystyle \sum_{k=0}^n \frac{f^{(k)}(x_0)}{k!}(x-x_0)^k \right) + \frac{f^{(n+1)}(c)}{(n+1)!}(x-x_0)^{n+1}")

where
![c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c "c")
is an (unknown) number between
![x_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_0 "x_0")
and
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x").

For example, suppose
![f(x)=\sin(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28x%29%3D%5Csin%28x%29 "f(x)=\sin(x)")
and
![x_0=0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_0%3D0 "x_0=0").
We can pre-compute some derivatives,
![f(0)=\sin(0)=0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%280%29%3D%5Csin%280%29%3D0 "f(0)=\sin(0)=0"),
![f^{(1)}(0)=\cos(0)=1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%5E%7B%281%29%7D%280%29%3D%5Ccos%280%29%3D1 "f^{(1)}(0)=\cos(0)=1"),
![f^{(2)}(0)=-\sin(0)=0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%5E%7B%282%29%7D%280%29%3D-%5Csin%280%29%3D0 "f^{(2)}(0)=-\sin(0)=0"),
![f^{(3)}(0)=-\cos(0)=-1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%5E%7B%283%29%7D%280%29%3D-%5Ccos%280%29%3D-1 "f^{(3)}(0)=-\cos(0)=-1"),
![f^{(4)}(0)=\sin(0)=0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%5E%7B%284%29%7D%280%29%3D%5Csin%280%29%3D0 "f^{(4)}(0)=\sin(0)=0"),
and so on. Then the Taylor polynomial to 4th order and the error term
are are

![0 + x + 0x^2 - \frac{x^3}{6} + 0x^4 + \frac{x^5}{120}\cos c.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0%20%2B%20x%20%2B%200x%5E2%20-%20%5Cfrac%7Bx%5E3%7D%7B6%7D%20%2B%200x%5E4%20%2B%20%5Cfrac%7Bx%5E5%7D%7B120%7D%5Ccos%20c. "0 + x + 0x^2 - \frac{x^3}{6} + 0x^4 + \frac{x^5}{120}\cos c.")

Here is a plot of the function and the subsequent Taylor approximations.

``` r
P1 <- function(x) {x}
P3 <- function(x) {x-x^3/6}
x = seq(from = -pi/2, to = pi/2, length = 200)
plot(x, sin(x), col = "black", type = "l")
lines(x, P1(x), col = "green")
lines(x, P3(x), col = "red")
```

![](coursenotes_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

We can also calculaute a bound on the error if we use the fourth degree
polynomial to approxiamte
![\sin(0.1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csin%280.1%29 "\sin(0.1)").
Note
![\sin](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csin "\sin")
is at most one in magnitude, so the error term is at most
![(0.1)^5/120\approx 8.333 \times 10^{-8}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%280.1%29%5E5%2F120%5Capprox%208.333%20%5Ctimes%2010%5E%7B-8%7D "(0.1)^5/120\approx 8.333 \times 10^{-8}").
The actual error achieved is

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

-   Convert between decimal and binary
-   Explain the form of double precision floating point numbers
-   Convert between real numbers and floating point numbers
-   Implement nearest neighbor rounding
-   Recognize loss of numerical significance
-   Implement strategies to reduce loss of significance
-   Calculate/estimate the error arising in numerical computations

## Binary Numbers

We usually work in the base 10 system, where each place in the number
represents a power of 10 with a coefficient of 0 through 9, *e.g.*,

![13.25 = 1 \times 10^1 + 3\times 10^0 + 2 \times 10^{-1} + 5 \times 10^{-2}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;13.25%20%3D%201%20%5Ctimes%2010%5E1%20%2B%203%5Ctimes%2010%5E0%20%2B%202%20%5Ctimes%2010%5E%7B-1%7D%20%2B%205%20%5Ctimes%2010%5E%7B-2%7D. "13.25 = 1 \times 10^1 + 3\times 10^0 + 2 \times 10^{-1} + 5 \times 10^{-2}.")

In the binary system, we instead use powers of 2 with coefficients of 0
through 1, *e.g.*,

![1101.01 = 1 \times 2^3 + 1\times 2^2 + 0 \times 2^1 + 1 \times 2^0 + 0 \times 2^{-1} + 1 \times 2^{-2}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1101.01%20%3D%201%20%5Ctimes%202%5E3%20%2B%201%5Ctimes%202%5E2%20%2B%200%20%5Ctimes%202%5E1%20%2B%201%20%5Ctimes%202%5E0%20%2B%200%20%5Ctimes%202%5E%7B-1%7D%20%2B%201%20%5Ctimes%202%5E%7B-2%7D. "1101.01 = 1 \times 2^3 + 1\times 2^2 + 0 \times 2^1 + 1 \times 2^0 + 0 \times 2^{-1} + 1 \times 2^{-2}.")

## Machine Numbers

Arithmetic performed on a machine can be different from what we do in
our heads. One fundamental difference is that computers can’t perfectly
store most numbers.

Most computers use a standard representation of real numbers called a
*floating point* number. Sometimes we just say *float* for short. We
represent the floating point version of a number
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
as
![\mathop{\mathrm{fl}}(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathop%7B%5Cmathrm%7Bfl%7D%7D%28x%29 "\mathop{\mathrm{fl}}(x)"),
and to reiterate the point above, it is not necessarily true that
![\mathop{\mathrm{fl}}(x)=x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathop%7B%5Cmathrm%7Bfl%7D%7D%28x%29%3Dx "\mathop{\mathrm{fl}}(x)=x").

We focus on double precision floats. A double precision float consists
of 64 bits (1’s and 0’s) and, expressed in base 10, is equal to

![(-1)^s 2^{c-1023}(1+f)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28-1%29%5Es%202%5E%7Bc-1023%7D%281%2Bf%29 "(-1)^s 2^{c-1023}(1+f)")

The parts of this number are:

| Part of number | Name     | Bits              | Notes                                             |
|----------------|----------|-------------------|---------------------------------------------------|
| s              | sign     | 1                 | Determines positive or negative                   |
| c              | exponent | 2 - 12 (11 bits)  | Subtract 1023 to correct exponential bias         |
| f              | mantissa | 13 - 64 (52 bits) | Interpreted as a binary fraction (after radix pt) |

For example, suppose our 64 bits are

![1 \| 10000000111 \| 1101000000000000000000000000000000000000000000000000](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%20%7C%2010000000111%20%7C%201101000000000000000000000000000000000000000000000000 "1 | 10000000111 | 1101000000000000000000000000000000000000000000000000")

Then expressed in base 10,
![s=1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s%3D1 "s=1"),
![c=1024+4+2+1=1031](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c%3D1024%2B4%2B2%2B1%3D1031 "c=1024+4+2+1=1031"),
and
![f=1/2+1/4+1/16=0.8125](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%3D1%2F2%2B1%2F4%2B1%2F16%3D0.8125 "f=1/2+1/4+1/16=0.8125").
So the floating point number is
![(-1)^1 \cdot 2^{1031-1023} \cdot 1.8125=-208](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28-1%29%5E1%20%5Ccdot%202%5E%7B1031-1023%7D%20%5Ccdot%201.8125%3D-208 "(-1)^1 \cdot 2^{1031-1023} \cdot 1.8125=-208").

Now let’s go the other way, namely converting from a decimal number to a
floating point number. In doing so, we’ll see that there can be error in
the computer storage of numbers. Let’s take the number
![59\\ 1/3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;59%5C%201%2F3 "59\ 1/3").

1.  Convert to binary. The number is
    ![111011.01010101\ldots](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;111011.01010101%5Cldots "111011.01010101\ldots").
2.  Write in normalized form, similar to scientific notation for
    base 10. In normalized form, there is a 1 before the radix point
    (nothing more, nothing less). This yelds
    ![1.11011\overline{01}\times 2^5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1.11011%5Coverline%7B01%7D%5Ctimes%202%5E5 "1.11011\overline{01}\times 2^5").
3.  Work backwords to figure out
    ![s](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s "s"),
    ![c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c "c"),
    and
    ![f](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f "f").
    Since the number is positive, the sign is
    ![s=0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s%3D0 "s=0").
    Since we have
    ![2^5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2%5E5 "2^5"),
    the exponent must equal
    ![5+1023 = 1028](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;5%2B1023%20%3D%201028 "5+1023 = 1028")
    in base 10, which is
    ![10000000100](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;10000000100 "10000000100").
    The mantissa
    ![f](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f "f")
    is
    ![11011\overline{01}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;11011%5Coverline%7B01%7D "11011\overline{01}").

The problem is that
![f](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f "f")
is repeating, but on a computer only holds 52 bits. The computer
truncates to 52 bits by nearest neigbor rounding, which works like this:

-   If the 53rd bit is a 0, drop the 53rd bit and everything after it.
-   If the 53rd bit is a 1 and there is any nonzero bit after it, add 1
    to the 52nd bit (and carry as necessary).
-   If the 53rd bit is a 1 and there are only zeros after it,
-   If the 52nd bit is 1, add 1 to it and carry as necessary
-   If the 52nd bit is a 0, do nothing (truncate after 52nd bit)

For our example above, in the mantissa
![f](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f "f")
we have
![11011](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;11011 "11011")
followed by 23 copies of the pattern
![01](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;01 "01")
followed by a
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0").
This means the 53rd bit is a
![1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1 "1"),
and there are nonzero bits after it, so we have to add
![1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1 "1")
to the 52nd. In base 10, we have

![f = 2^{-1} + 2^{-2} + 2^{-4} + 2^{-5} + 2^{-7} + 2^{-9} + + \ldots + 2^{-51} + 2^{-52} = 0.85416666666666674068.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%20%3D%202%5E%7B-1%7D%20%2B%202%5E%7B-2%7D%20%2B%202%5E%7B-4%7D%20%2B%202%5E%7B-5%7D%20%2B%202%5E%7B-7%7D%20%2B%202%5E%7B-9%7D%20%2B%20%2B%20%5Cldots%20%2B%202%5E%7B-51%7D%20%2B%202%5E%7B-52%7D%20%3D%200.85416666666666674068. "f = 2^{-1} + 2^{-2} + 2^{-4} + 2^{-5} + 2^{-7} + 2^{-9} + + \ldots + 2^{-51} + 2^{-52} = 0.85416666666666674068.")

Therefore, on the machine, our number is stored as
![2^5(1.85416666666666674068)=59.333333333333335702 \neq 59\\ 1/3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2%5E5%281.85416666666666674068%29%3D59.333333333333335702%20%5Cneq%2059%5C%201%2F3 "2^5(1.85416666666666674068)=59.333333333333335702 \neq 59\ 1/3").

The relative error in the storage of this number is

![\frac{\|59\\ 1/3 -  59.333333333333335702\|}{59\\ 1/3} \approx -4 \times 10^{-17}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B%7C59%5C%201%2F3%20-%20%2059.333333333333335702%7C%7D%7B59%5C%201%2F3%7D%20%5Capprox%20-4%20%5Ctimes%2010%5E%7B-17%7D. "\frac{|59\ 1/3 -  59.333333333333335702|}{59\ 1/3} \approx -4 \times 10^{-17}.")

An important number is *machine epsilon*, denoted
![\epsilon\_{mach}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7Bmach%7D "\epsilon_{mach}"),
which is the distance from the number 1 to the next smallest number that
can be represented exactly in floating point form. This distance is
![\epsilon\_{mach}=2^{-52}\approx 2.2 \times 10^{-16}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7Bmach%7D%3D2%5E%7B-52%7D%5Capprox%202.2%20%5Ctimes%2010%5E%7B-16%7D "\epsilon_{mach}=2^{-52}\approx 2.2 \times 10^{-16}").

## Machine Addition

Machine addition is defined as

![\mathop{\mathrm{fl}}(x+y) = \mathop{\mathrm{fl}}(\mathop{\mathrm{fl}}(x)+\mathop{\mathrm{fl}}(y)).](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathop%7B%5Cmathrm%7Bfl%7D%7D%28x%2By%29%20%3D%20%5Cmathop%7B%5Cmathrm%7Bfl%7D%7D%28%5Cmathop%7B%5Cmathrm%7Bfl%7D%7D%28x%29%2B%5Cmathop%7B%5Cmathrm%7Bfl%7D%7D%28y%29%29. "\mathop{\mathrm{fl}}(x+y) = \mathop{\mathrm{fl}}(\mathop{\mathrm{fl}}(x)+\mathop{\mathrm{fl}}(y)).")

That is, take
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
and
![y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y "y"),
convert each to a machine number, add them (exactly, since more
registers are available for this operation), and convert the result to a
machine number. Subtraction is just addition with a negative sign.

For example, let’s add
![1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1 "1")
and
![2^{-53}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2%5E%7B-53%7D "2^{-53}").
Well,
![\mathop{\mathrm{fl}}(1)=1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathop%7B%5Cmathrm%7Bfl%7D%7D%281%29%3D1 "\mathop{\mathrm{fl}}(1)=1")
and
![\mathop{\mathrm{fl}}(2^{-53})=2^{-53}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathop%7B%5Cmathrm%7Bfl%7D%7D%282%5E%7B-53%7D%29%3D2%5E%7B-53%7D "\mathop{\mathrm{fl}}(2^{-53})=2^{-53}").
The (exact) sum is
![1+2^{-53}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%2B2%5E%7B-53%7D "1+2^{-53}")
but due to the rounding rules on machines,
![\mathop{\mathrm{fl}}(1+2^{-53})=1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathop%7B%5Cmathrm%7Bfl%7D%7D%281%2B2%5E%7B-53%7D%29%3D1 "\mathop{\mathrm{fl}}(1+2^{-53})=1").
Therefore, on a machine, the sum is
![1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1 "1").
Let’s try:

``` r
1 + 2^(-53)
```

    ## [1] 1

## Loss of Significance

We’ve seen that computer storage of numbers can have error, and
therefore arithmetic can have error. This error is sometimes called
*loss of significance* and the most dangerous operation is subtraction
of nearly equal numbers. Consider the expression
![(1-\cos x)/\sin^2 x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%281-%5Ccos%20x%29%2F%5Csin%5E2%20x "(1-\cos x)/\sin^2 x"),
which can also be written as
![1/(1+\cos x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%2F%281%2B%5Ccos%20x%29 "1/(1+\cos x)").
We compute this both ways for
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
decreasing from
![1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1 "1").

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

The result becomes drastically wrong for small
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
because of subtraction of nearly equal numbers. In such problems, be
aware and seek alternative ways to represent the necessary computation,
as we did with the second option above.

# Fundamentals of Linear Systems

## Big Picture

We begin considering the solution of systems of linear equations,
![\textbf{A}\textbf{x}=\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%3D%5Ctextbf%7Bb%7D "\textbf{A}\textbf{x}=\textbf{b}").
Linear systems arise in analysis of different equations (modeling
macromolecules, electromagnetics, heat flow, wave motion, structural
engineering, a million other examples), in curve fitting, in
optimization, and many other applications. To understand solution of
linear systems, it’s helpful to recall some fundamental ideas of linear
algebra. Additionally, as we consider the solution of problems on a
computer, we have to think about the effect of small errors on the
solution of the problem. This is called conditioning. In the linear
algebra setting, conditioning is intimately related to matrix and vector
norms.

## Goals

-   Interpret linear systems geometrically
-   Define linear algebra terms including inverse, determinant,
    eigenvalues, nullspace, linear independence, span, and image
-   Give intuitive explanations for the “big theorem” of linear algebra
-   Perform the steps of Gaussian elimination
-   Establish the computational complexity of Gaussian elimination
-   Define forward and backward error
-   Define vector norms
-   Define matrix norms
-   Define condition number for solution of
    ![\textbf{A}\textbf{x}=\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%3D%5Ctextbf%7Bb%7D "\textbf{A}\textbf{x}=\textbf{b}")

## Linear Systems

For concreteness let’s work in three dimensions and let’s consider the
system

$$
=
. $$

Let’s interpret this system two ways.

**Intersecting hyperplanes.** Carrying out the matrix multiplication, we
write the equations

![\begin{eqnarray}
2x_1 + 4x_2 -2x_3 & = & 8\\\\
x_1 + 4x_2 -3x_3 & = & 8\\\\
-2x_1 -6x_2 + 7x_3 & = & -3
\end{eqnarray}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Beqnarray%7D%0A2x_1%20%2B%204x_2%20-2x_3%20%26%20%3D%20%26%208%5C%5C%0Ax_1%20%2B%204x_2%20-3x_3%20%26%20%3D%20%26%208%5C%5C%0A-2x_1%20-6x_2%20%2B%207x_3%20%26%20%3D%20%26%20-3%0A%5Cend%7Beqnarray%7D "\begin{eqnarray}
2x_1 + 4x_2 -2x_3 & = & 8\\
x_1 + 4x_2 -3x_3 & = & 8\\
-2x_1 -6x_2 + 7x_3 & = & -3
\end{eqnarray}")

This form suggests thinking of the set of points that are at the
intersection of these three planes, which could be the empty set, or a
point, or a line, or a plane.

**Vector spans.** Write the equation in terms of its columns, as

![x_1 \begin{pmatrix} 2 \\\\ 1 \\\\ -2 \end{pmatrix} + x_2 \begin{pmatrix} 4 \\\\ 4 \\\\ -6 \end{pmatrix} + x_3 \begin{pmatrix} -2 \\\\ -3 \\\\ 7 \end{pmatrix} = \begin{pmatrix} 8 \\\\ 8 \\\\ -3 \end{pmatrix}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_1%20%5Cbegin%7Bpmatrix%7D%202%20%5C%5C%201%20%5C%5C%20-2%20%5Cend%7Bpmatrix%7D%20%2B%20x_2%20%5Cbegin%7Bpmatrix%7D%204%20%5C%5C%204%20%5C%5C%20-6%20%5Cend%7Bpmatrix%7D%20%2B%20x_3%20%5Cbegin%7Bpmatrix%7D%20-2%20%5C%5C%20-3%20%5C%5C%207%20%5Cend%7Bpmatrix%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%208%20%5C%5C%208%20%5C%5C%20-3%20%5Cend%7Bpmatrix%7D "x_1 \begin{pmatrix} 2 \\ 1 \\ -2 \end{pmatrix} + x_2 \begin{pmatrix} 4 \\ 4 \\ -6 \end{pmatrix} + x_3 \begin{pmatrix} -2 \\ -3 \\ 7 \end{pmatrix} = \begin{pmatrix} 8 \\ 8 \\ -3 \end{pmatrix}")

This form suggests thinking of the linear combination of three basis
vectors necessary to reach a particular target vector. There could be 0,
1, or inifinity depending on the arrangement of those vectors.

## Linear Algebra Review and the Big Theorem

Let’s consider an
![n \times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%20n "n \times n")
matrix
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}").

1.  **Invertible** means that the inverse
    ![\textbf{A}^{-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5E%7B-1%7D "\textbf{A}^{-1}")
    exists. This matrix satisfies
    ![\textbf{A}\textbf{A}^{-1} = \textbf{I}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Ctextbf%7BA%7D%5E%7B-1%7D%20%3D%20%5Ctextbf%7BI%7D "\textbf{A}\textbf{A}^{-1} = \textbf{I}"),
    where
    ![\textbf{I}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BI%7D "\textbf{I}")
    is the
    ![n \times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%20n "n \times n")
    identity matrix.

2.  The notation
    ![\det](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cdet "\det")
    means **determinant**. Think of it as a scaling factor for the
    transformation defined by a matrix. That is, multiplication by
    ![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
    can cause a region to contract
    (![\|\det \textbf{A}\| \< 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%5Cdet%20%5Ctextbf%7BA%7D%7C%20%3C%201 "|\det \textbf{A}| < 1"))
    or expand
    (![\|\det \textbf{A}\| \> 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%5Cdet%20%5Ctextbf%7BA%7D%7C%20%3E%201 "|\det \textbf{A}| > 1"))
    and/or reflect
    (![\det \textbf{A} \< 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cdet%20%5Ctextbf%7BA%7D%20%3C%200 "\det \textbf{A} < 0")).
    As an example, let

![\textbf{A}=\begin{pmatrix}a & b\\\c & d\end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%3D%5Cbegin%7Bpmatrix%7Da%20%26%20b%5C%5Cc%20%26%20d%5Cend%7Bpmatrix%7D. "\textbf{A}=\begin{pmatrix}a & b\\c & d\end{pmatrix}.")

For this example,

![\textbf{A}^{-1} = \frac{1}{ad-bc}\begin{pmatrix}d & -b\\\\-c & a\end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5E%7B-1%7D%20%3D%20%5Cfrac%7B1%7D%7Bad-bc%7D%5Cbegin%7Bpmatrix%7Dd%20%26%20-b%5C%5C-c%20%26%20a%5Cend%7Bpmatrix%7D. "\textbf{A}^{-1} = \frac{1}{ad-bc}\begin{pmatrix}d & -b\\-c & a\end{pmatrix}.")

So
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
invertible
![\iff \det \textbf{A} \neq 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ciff%20%5Cdet%20%5Ctextbf%7BA%7D%20%5Cneq%200 "\iff \det \textbf{A} \neq 0").

3.  The **eigenvalues**
    ![\lambda_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_i "\lambda_i")
    of
    ![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
    satisfy

![\textbf{A}\mathbf{v_i}=\lambda_i \mathbf{v_i},](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Cmathbf%7Bv_i%7D%3D%5Clambda_i%20%5Cmathbf%7Bv_i%7D%2C "\textbf{A}\mathbf{v_i}=\lambda_i \mathbf{v_i},")

where
![\mathbf{v_i}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv_i%7D "\mathbf{v_i}")
are the **eigenvectors**. You can prove that

![\prod_i \lambda_i = \det \textbf{A},](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cprod_i%20%5Clambda_i%20%3D%20%5Cdet%20%5Ctextbf%7BA%7D%2C "\prod_i \lambda_i = \det \textbf{A},")

so no
![\lambda_i = 0 \iff \det \textbf{A} \neq 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_i%20%3D%200%20%5Ciff%20%5Cdet%20%5Ctextbf%7BA%7D%20%5Cneq%200 "\lambda_i = 0 \iff \det \textbf{A} \neq 0").

4.  ![\textbf{A} \mathbf{z} \neq 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%5Cmathbf%7Bz%7D%20%5Cneq%200 "\textbf{A} \mathbf{z} \neq 0")
    for all
    ![\mathbf{z} \in \mathbb{R}^n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bz%7D%20%5Cin%20%5Cmathbb%7BR%7D%5En "\mathbf{z} \in \mathbb{R}^n")
    except
    ![\mathbf{z}=0 \iff \textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bz%7D%3D0%20%5Ciff%20%5Ctextbf%7BA%7D "\mathbf{z}=0 \iff \textbf{A}")
    is invertible. Why? If
    ![\textbf{A} \mathbf{z} = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%5Cmathbf%7Bz%7D%20%3D%200 "\textbf{A} \mathbf{z} = 0")
    for
    ![\mathbf{z} \neq 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bz%7D%20%5Cneq%200 "\mathbf{z} \neq 0"),
    then
    ![\textbf{A} \mathbf{z} = 0 \mathbf{z}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%5Cmathbf%7Bz%7D%20%3D%200%20%5Cmathbf%7Bz%7D "\textbf{A} \mathbf{z} = 0 \mathbf{z}"),
    so 0 is an eigenvalue. But for
    ![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
    to be invertible, we know
    ![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0")
    can’t be an eigenvalue. (If you are proof oriented you might notice
    that the implication needs to be shown both ways, but I am trying
    purposely not to prove here – just to give you some intuition.)

5.  **Nullspace** just means the parts of
    ![\mathbb{R}^n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BR%7D%5En "\mathbb{R}^n")
    that gets mapped to
    ![\mathbf{0}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7B0%7D "\mathbf{0}")
    by
    ![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}").
    Another name for nullspace is **kernel**. Mathematically, the
    nullspace is all the vectors
    ![\mathbf{v}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D "\mathbf{v}")
    for which
    ![\textbf{A} \mathbf{v}=0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%5Cmathbf%7Bv%7D%3D0 "\textbf{A} \mathbf{v}=0").
    So for an invertible matrix
    ![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}"),
    the nullspace is
    ![\mathbf{0}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7B0%7D "\mathbf{0}").
    This is basically assigning a definition to the previous point
    (above).

6.  If vectors are **linearly independent**, it means that none of the
    vectors can be written as a linear combination of the others. If
    ![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
    is invertible, then its columns are linearly independent. Why? If
    the columns were linearly dependent, you could take a linear
    combination of them to reach
    ![\mathbf{0}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7B0%7D "\mathbf{0}")
    nontrivially, for instance,

![\begin{pmatrix}1 & -2\\\2 & -4\end{pmatrix}\begin{pmatrix}x_1 \\\\ x_2 \end{pmatrix}=\begin{pmatrix}0 \\\\ 0\end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Bpmatrix%7D1%20%26%20-2%5C%5C2%20%26%20-4%5Cend%7Bpmatrix%7D%5Cbegin%7Bpmatrix%7Dx_1%20%5C%5C%20x_2%20%5Cend%7Bpmatrix%7D%3D%5Cbegin%7Bpmatrix%7D0%20%5C%5C%200%5Cend%7Bpmatrix%7D. "\begin{pmatrix}1 & -2\\2 & -4\end{pmatrix}\begin{pmatrix}x_1 \\ x_2 \end{pmatrix}=\begin{pmatrix}0 \\ 0\end{pmatrix}.")

This violates our previous condition (above) about the nullspace only
being
![\mathbf{0}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7B0%7D "\mathbf{0}").

7.  **Span** means the set of points reachable by taking linear
    combinations of a set of vectors. If you have
    ![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
    linearly independent vectors in
    ![\mathbb{R}^n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BR%7D%5En "\mathbb{R}^n"),
    they span
    ![\mathbb{R}^n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BR%7D%5En "\mathbb{R}^n").
    The **rank** of
    ![\mathbb{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BA%7D "\mathbb{A}")
    is just the dimension of the space spanned by the column vectors.

8.  The **image** of
    ![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
    means all the points that
    ![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
    can map to. This is synonymous with the point above: it’s the span
    of the columns of
    ![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}").
    If the rank of
    ![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
    is
    ![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n"),
    then the column vectors are linearly independent so they span
    ![\mathbb{R}^n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BR%7D%5En "\mathbb{R}^n"),
    so the image of
    ![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
    is
    ![\mathbb{R}^n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BR%7D%5En "\mathbb{R}^n").

9.  ![\textbf{A}\textbf{x}=\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%3D%5Ctextbf%7Bb%7D "\textbf{A}\textbf{x}=\textbf{b}")
    has a unique solution for all
    ![\textbf{b} \iff](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bb%7D%20%5Ciff "\textbf{b} \iff")
    ![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
    is invertible, since you can left multiply by
    ![\textbf{A}^{-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5E%7B-1%7D "\textbf{A}^{-1}").

## Gaussian Elimination

First, we define **row echelon form**. A matrix is in row echelon for if

-   all nonzero rows (rows with at least one nonzero element) are above
    any rows of all zeroes
-   the leading coefficient (the first nonzero number from the left,
    also called the **pivot**) of a nonzero row is always strictly to
    the right of the leading coefficient of the row above it

Gaussian elimination refers to transforming a matrix to row echelon form
by applying the following operations, which do not change the solution
set:

-   swap two rows
-   multiply a row by a nonzero scalar
-   Add one row to a scalar multiple of another

To solve a linear system
![\textbf{A}\textbf{x}=\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%3D%5Ctextbf%7Bb%7D "\textbf{A}\textbf{x}=\textbf{b}"),
write it as an augmented matrix, reduce it to row echelon form, and then
use back substitution to solve.

For example, take

![\textbf{A} = 
\begin{pmatrix}
1 & 3 & 1\\\\
1 & 1 & -1 \\\\
3 & 11 & 5
\end{pmatrix}, \quad
\textbf{b} = \begin{pmatrix}
9 \\\\ 1 \\\\ 35
\end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%3D%20%0A%5Cbegin%7Bpmatrix%7D%0A1%20%26%203%20%26%201%5C%5C%0A1%20%26%201%20%26%20-1%20%5C%5C%0A3%20%26%2011%20%26%205%0A%5Cend%7Bpmatrix%7D%2C%20%5Cquad%0A%5Ctextbf%7Bb%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%0A9%20%5C%5C%201%20%5C%5C%2035%0A%5Cend%7Bpmatrix%7D. "\textbf{A} = 
\begin{pmatrix}
1 & 3 & 1\\
1 & 1 & -1 \\
3 & 11 & 5
\end{pmatrix}, \quad
\textbf{b} = \begin{pmatrix}
9 \\ 1 \\ 35
\end{pmatrix}.")

Write the augmented matrix

![\mathbf{A_a} = \begin{pmatrix}
1 & 3 & 1 & 9 \\\\
1 & 1 & -1 & 1\\\\
3 & 11 & 5 & 35
\end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA_a%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%0A1%20%26%203%20%26%201%20%26%209%20%5C%5C%0A1%20%26%201%20%26%20-1%20%26%201%5C%5C%0A3%20%26%2011%20%26%205%20%26%2035%0A%5Cend%7Bpmatrix%7D. "\mathbf{A_a} = \begin{pmatrix}
1 & 3 & 1 & 9 \\
1 & 1 & -1 & 1\\
3 & 11 & 5 & 35
\end{pmatrix}.")

Apply
![II \leftarrow II - I](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;II%20%5Cleftarrow%20II%20-%20I "II \leftarrow II - I")
and
![III \leftarrow III - 3 I](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;III%20%5Cleftarrow%20III%20-%203%20I "III \leftarrow III - 3 I")
to obtain:

![\mathbf{A_a} = \begin{pmatrix}
1 & 3 & 1 & 9 \\\\
0 & -2 & -2 & -8\\\\
0 & 2 & 2 & 8
\end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA_a%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%0A1%20%26%203%20%26%201%20%26%209%20%5C%5C%0A0%20%26%20-2%20%26%20-2%20%26%20-8%5C%5C%0A0%20%26%202%20%26%202%20%26%208%0A%5Cend%7Bpmatrix%7D. "\mathbf{A_a} = \begin{pmatrix}
1 & 3 & 1 & 9 \\
0 & -2 & -2 & -8\\
0 & 2 & 2 & 8
\end{pmatrix}.")

Apply
![III \leftarrow III + II](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;III%20%5Cleftarrow%20III%20%2B%20II "III \leftarrow III + II")
to obtain:

![\mathbf{A_a} = \begin{pmatrix}
1 & 3 & 1 & 9 \\\\
0 & -2 & -2 & -8\\\\
0 & 0 & 0 & 0
\end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA_a%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%0A1%20%26%203%20%26%201%20%26%209%20%5C%5C%0A0%20%26%20-2%20%26%20-2%20%26%20-8%5C%5C%0A0%20%26%200%20%26%200%20%26%200%0A%5Cend%7Bpmatrix%7D. "\mathbf{A_a} = \begin{pmatrix}
1 & 3 & 1 & 9 \\
0 & -2 & -2 & -8\\
0 & 0 & 0 & 0
\end{pmatrix}.")

The bottom row tells us nothing. The second row tells us there is a free
variable, which we take to be
![x_3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_3 "x_3").
So we solve this equation for
![x_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_2 "x_2"),
finding
![x_2 = 4 -x_3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_2%20%3D%204%20-x_3 "x_2 = 4 -x_3").
This is called back substitution. Then we do back substitution on the
top row, from which we find
![x_1 = 9 - x_3 - 3 x_2 = 9 - x_3 - 3(4 - x_3) = -3 + 2x_3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_1%20%3D%209%20-%20x_3%20-%203%20x_2%20%3D%209%20-%20x_3%20-%203%284%20-%20x_3%29%20%3D%20-3%20%2B%202x_3 "x_1 = 9 - x_3 - 3 x_2 = 9 - x_3 - 3(4 - x_3) = -3 + 2x_3").
Therefore, the solution is

![\textbf{x} = \begin{pmatrix}  -3 + 2x_3 \\\\ 4 - x_3 \\\\ x_3 \end{pmatrix} = \begin{pmatrix}  -3 \\\\ 4\\\\ 0 \end{pmatrix} + x_3 \begin{pmatrix}  2 \\\\ -1 \\\\ 1 \end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20%20-3%20%2B%202x_3%20%5C%5C%204%20-%20x_3%20%5C%5C%20x_3%20%5Cend%7Bpmatrix%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20%20-3%20%5C%5C%204%5C%5C%200%20%5Cend%7Bpmatrix%7D%20%2B%20x_3%20%5Cbegin%7Bpmatrix%7D%20%202%20%5C%5C%20-1%20%5C%5C%201%20%5Cend%7Bpmatrix%7D. "\textbf{x} = \begin{pmatrix}  -3 + 2x_3 \\ 4 - x_3 \\ x_3 \end{pmatrix} = \begin{pmatrix}  -3 \\ 4\\ 0 \end{pmatrix} + x_3 \begin{pmatrix}  2 \\ -1 \\ 1 \end{pmatrix}.")

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
the method. For solving
![\textbf{A}\textbf{x}=\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%3D%5Ctextbf%7Bb%7D "\textbf{A}\textbf{x}=\textbf{b}"),
with
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
an
![n \times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%20n "n \times n")
matrix, we write the complexity in terms of
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n").
Then, we are usually concerned with the behavior of the operation count
for
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
large, so we might retain just the leading term in
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
as an approximation, or even ignore the coefficient in front that
leading term.

For Gaussian elimination, we have to compute the complexity of the two
stages.

1.  Reduce to echelon form. This takes

![\frac{2}{3}n^3 + \frac{1}{2}n^2-\frac{7}{6}n = \mathcal{O}(n^3)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B2%7D%7B3%7Dn%5E3%20%2B%20%5Cfrac%7B1%7D%7B2%7Dn%5E2-%5Cfrac%7B7%7D%7B6%7Dn%20%3D%20%5Cmathcal%7BO%7D%28n%5E3%29 "\frac{2}{3}n^3 + \frac{1}{2}n^2-\frac{7}{6}n = \mathcal{O}(n^3)")

operations.

2.  Back substitute. This takes
    ![n^2 = \mathcal{O}(n^2)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%5E2%20%3D%20%5Cmathcal%7BO%7D%28n%5E2%29 "n^2 = \mathcal{O}(n^2)")
    operations.

Back substitution is comptutationally cheap compared to row reduction.
For large enough
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n"),
the back substitution step is negligible since
![n^3 \gg n^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%5E3%20%5Cgg%20n%5E2 "n^3 \gg n^2").

We can use these operation counts to make estimates of how long
calculations should take.

For example, suppose row reduction on a
![500 \times 500](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;500%20%5Ctimes%20500 "500 \times 500")
matrix takes 1 second. How long does back substitution take? Well, we
just use leading terms in the complexity. Since
![n=500](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%3D500 "n=500"),
we have
![(2/3)\*500^3 t =1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%282%2F3%29%2A500%5E3%20t%20%3D1 "(2/3)*500^3 t =1")
where
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t")
is time per operation. Solving,
![t = 1.2 \times 10^{-8}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t%20%3D%201.2%20%5Ctimes%2010%5E%7B-8%7D "t = 1.2 \times 10^{-8}").
Back substitution takes
![n^2=500^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%5E2%3D500%5E2 "n^2=500^2")
operations, so the total time is
![n^2t=500^2\times 1.2 \times 10^{-8}=0.003](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%5E2t%3D500%5E2%5Ctimes%201.2%20%5Ctimes%2010%5E%7B-8%7D%3D0.003 "n^2t=500^2\times 1.2 \times 10^{-8}=0.003")
seconds.

Let’s test scaling of the reduction step on Chad’s machine for a random
matrix.

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
    ## 4.79054054054

## Forward and Backward Error

In the solution of a computational problem, **forward error** is the
difference between the exact and computed solution, and **backwards
error** is the difference between the original problem and the so-called
modified problem that the approximate solution satisfies. This probably
sounds abstract, so let’s make it concrete in the cases of a
root-finding problem and a linear algebra problem.

Suppose we want to solve
![\textbf{A}\textbf{x}=\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%3D%5Ctextbf%7Bb%7D "\textbf{A}\textbf{x}=\textbf{b}").
The true solution is
![\textbf{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D "\textbf{x}")
but our computational method finds an approximate solution
![\textbf{x}\_a](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D_a "\textbf{x}_a").
The forward error is the distance between the two solutions, that is,
![\|\|\textbf{x}-\mathbf{x_a}\|\|](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%7C%5Ctextbf%7Bx%7D-%5Cmathbf%7Bx_a%7D%7C%7C "||\textbf{x}-\mathbf{x_a}||").
The backward error is the distance between what the matrix outputs when
applied to those solutions, that is,
![\|\|\textbf{A}\textbf{x}-\textbf{A}\textbf{x}\_a\|\|=\|\|\textbf{b}-\textbf{A}\textbf{x}\_a\|\|](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%7C%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D-%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D_a%7C%7C%3D%7C%7C%5Ctextbf%7Bb%7D-%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D_a%7C%7C "||\textbf{A}\textbf{x}-\textbf{A}\textbf{x}_a||=||\textbf{b}-\textbf{A}\textbf{x}_a||").
Distance here is the length of the difference between two quantities.

Notice that we haven’t specified what distance means! This is why we
need to define vector and matrix norms.

## Vector Norms

A vector norm is a rule that assigns a real number to every vector.
Intuitively, it measures length. There are a bunch of requirements that
this rule must satisfy in order to be a norm, but rather than stating
those requirements, I’m going to just tell you some practicalities.

The vector norm we’ll work with is called the
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")-norm.
The
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")-norm
for
![1 \le p \le \infty](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%20%5Cle%20p%20%5Cle%20%5Cinfty "1 \le p \le \infty")
is defined as

![\|\| \textbf{x} \|\|\_p = \left({\|x_1\|^p + \|x_2\|^p + \cdots + \|x_n\|^p} \right)^{1/p}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%7C%20%5Ctextbf%7Bx%7D%20%7C%7C_p%20%3D%20%5Cleft%28%7B%7Cx_1%7C%5Ep%20%2B%20%7Cx_2%7C%5Ep%20%2B%20%5Ccdots%20%2B%20%7Cx_n%7C%5Ep%7D%20%5Cright%29%5E%7B1%2Fp%7D. "|| \textbf{x} ||_p = \left({|x_1|^p + |x_2|^p + \cdots + |x_n|^p} \right)^{1/p}.")

The three most common
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")-norms
are
![p = 1, 2, \infty](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%20%3D%201%2C%202%2C%20%5Cinfty "p = 1, 2, \infty")
since they are the easiest to compute with and in some sense are the
most natural:

-   ![p=1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%3D1 "p=1")
    (the Manhattan or taxicab norm) $$ \|\| \|\|\_1 =\|x_1\| +
    \|x_2\| + + \|x_n\|

![-   $p =2$ (the Euclidean norm)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-%20%20%20%24p%20%3D2%24%20%28the%20Euclidean%20norm%29 "-   $p =2$ (the Euclidean norm)")

\|\| \|\|\_2 = = $$

-   ![p=\infty](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%3D%5Cinfty "p=\infty")

![\|\| \textbf{x} \|\|\_\infty = \max{\left(\| x_1\|, \| x_2\|,  \ldots, \|x_n\|\right)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%7C%20%5Ctextbf%7Bx%7D%20%7C%7C_%5Cinfty%20%3D%20%5Cmax%7B%5Cleft%28%7C%20x_1%7C%2C%20%7C%20x_2%7C%2C%20%20%5Cldots%2C%20%7Cx_n%7C%5Cright%29%7D "|| \textbf{x} ||_\infty = \max{\left(| x_1|, | x_2|,  \ldots, |x_n|\right)}")

For
![p = \infty](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%20%3D%20%5Cinfty "p = \infty")
it takes a little analysis to show why the computational definition is
what it is, but a numerical study is usually convincing. We can use the
`Norm` command

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
![\mathbb{R}^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BR%7D%5E2 "\mathbb{R}^2")
for various values of
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p"),
that is, the set of points that are unit distance from the origin.
![](UnitCircleGrid.png)

## Matrix Norms

The matrix
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")-norm
is closely related to the vector
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")-norm,
and is given by

![\|\|\textbf{A}\|\|\_p = \max\_{\textbf{x} \not = \mathbf{ 0}} \frac{ \|\| \textbf{A} \textbf{x}\|\|\_p} { \|\|\textbf{x}\|\|\_p} = \max\_{\|\|\textbf{x}\|\|\_p  = 1}  \|\| \textbf{A} \textbf{x}\|\|\_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%7C%5Ctextbf%7BA%7D%7C%7C_p%20%3D%20%5Cmax_%7B%5Ctextbf%7Bx%7D%20%5Cnot%20%3D%20%5Cmathbf%7B%200%7D%7D%20%5Cfrac%7B%20%7C%7C%20%5Ctextbf%7BA%7D%20%5Ctextbf%7Bx%7D%7C%7C_p%7D%20%7B%20%7C%7C%5Ctextbf%7Bx%7D%7C%7C_p%7D%20%3D%20%5Cmax_%7B%7C%7C%5Ctextbf%7Bx%7D%7C%7C_p%20%20%3D%201%7D%20%20%7C%7C%20%5Ctextbf%7BA%7D%20%5Ctextbf%7Bx%7D%7C%7C_p "||\textbf{A}||_p = \max_{\textbf{x} \not = \mathbf{ 0}} \frac{ || \textbf{A} \textbf{x}||_p} { ||\textbf{x}||_p} = \max_{||\textbf{x}||_p  = 1}  || \textbf{A} \textbf{x}||_p")

The matrix
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")-norm
says: apply
![A](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A "A")
to the unit sphere, and
![\|\|\textbf{A}\|\|\_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%7C%5Ctextbf%7BA%7D%7C%7C_p "||\textbf{A}||_p")
is the length of the vector that is farthest from the origin. This is
not trivial to calculate! You have an infinite number of vectors to
consider.

Fortunately, just like for the vector case, the matrix
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")-norm
has a few special values of
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")
for which it is easy to compute. We have:

-   ![p=1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%3D1 "p=1")

![\|\|\textbf{A} \|\|\_1 = \displaystyle{\max\_{1 \le j \le n} \sum\_{i=1}^n \|a\_{ij}\|} = \text{maximum absolute column sum}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%7C%5Ctextbf%7BA%7D%20%7C%7C_1%20%3D%20%5Cdisplaystyle%7B%5Cmax_%7B1%20%5Cle%20j%20%5Cle%20n%7D%20%5Csum_%7Bi%3D1%7D%5En%20%7Ca_%7Bij%7D%7C%7D%20%3D%20%5Ctext%7Bmaximum%20absolute%20column%20sum%7D "||\textbf{A} ||_1 = \displaystyle{\max_{1 \le j \le n} \sum_{i=1}^n |a_{ij}|} = \text{maximum absolute column sum}")

-   ![p=2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%3D2 "p=2")

![\|\| \textbf{A} \|\|\_2 = \sqrt{\max\\{ \text{eigenvalue}(A^TA) \\} }](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%7C%20%5Ctextbf%7BA%7D%20%7C%7C_2%20%3D%20%5Csqrt%7B%5Cmax%5C%7B%20%5Ctext%7Beigenvalue%7D%28A%5ETA%29%20%5C%7D%20%7D "|| \textbf{A} ||_2 = \sqrt{\max\{ \text{eigenvalue}(A^TA) \} }")

-   ![p=\infty](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%3D%5Cinfty "p=\infty")

![\|\| \textbf{A} \|\|\_\infty = \displaystyle{\max\_{1 \le i \le n} \sum\_{j=1}^n \|a\_{ij}\|} = \text{maximum absolute row sum}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%7C%20%5Ctextbf%7BA%7D%20%7C%7C_%5Cinfty%20%3D%20%5Cdisplaystyle%7B%5Cmax_%7B1%20%5Cle%20i%20%5Cle%20n%7D%20%5Csum_%7Bj%3D1%7D%5En%20%7Ca_%7Bij%7D%7C%7D%20%3D%20%5Ctext%7Bmaximum%20absolute%20row%20sum%7D "|| \textbf{A} ||_\infty = \displaystyle{\max_{1 \le i \le n} \sum_{j=1}^n |a_{ij}|} = \text{maximum absolute row sum}")

To see why these definitions are true requires some analysis. If you are
interested, I am happy to point you to proofs.

You can calculate the
![1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1 "1"),
![2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2 "2"),
and
![\infty](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cinfty "\infty")
matrix norms using the R command `norm`.

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

![\|\|\textbf{A}\textbf{x}\|\|\_p \leq \|\|\textbf{A}\|\|\_p \|\|\textbf{x}\|\|\_p.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%7C%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%7C%7C_p%20%5Cleq%20%7C%7C%5Ctextbf%7BA%7D%7C%7C_p%20%7C%7C%5Ctextbf%7Bx%7D%7C%7C_p. "||\textbf{A}\textbf{x}||_p \leq ||\textbf{A}||_p ||\textbf{x}||_p.")

To see this, we start with the right hand side and note

![\|\|\textbf{A}\|\|\_p \|\|\textbf{x}\|\|\_p = \left( \max\_{\mathbf{y} \not = \mathbf{ 0}} \frac{ \|\| \textbf{A} \mathbf{y}\|\|\_p} { \|\|\mathbf{y}\|\|\_p} \right) \|\|\textbf{x}\|\|\_p \geq \frac{ \|\| \textbf{A} \textbf{x}\|\|\_p} { \|\|\textbf{x}\|\|\_p} \|\|\textbf{x}\|\|\_p = \|\| \textbf{A} \textbf{x}\|\|\_p.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%7C%5Ctextbf%7BA%7D%7C%7C_p%20%7C%7C%5Ctextbf%7Bx%7D%7C%7C_p%20%3D%20%5Cleft%28%20%5Cmax_%7B%5Cmathbf%7By%7D%20%5Cnot%20%3D%20%5Cmathbf%7B%200%7D%7D%20%5Cfrac%7B%20%7C%7C%20%5Ctextbf%7BA%7D%20%5Cmathbf%7By%7D%7C%7C_p%7D%20%7B%20%7C%7C%5Cmathbf%7By%7D%7C%7C_p%7D%20%5Cright%29%20%7C%7C%5Ctextbf%7Bx%7D%7C%7C_p%20%5Cgeq%20%5Cfrac%7B%20%7C%7C%20%5Ctextbf%7BA%7D%20%5Ctextbf%7Bx%7D%7C%7C_p%7D%20%7B%20%7C%7C%5Ctextbf%7Bx%7D%7C%7C_p%7D%20%7C%7C%5Ctextbf%7Bx%7D%7C%7C_p%20%3D%20%7C%7C%20%5Ctextbf%7BA%7D%20%5Ctextbf%7Bx%7D%7C%7C_p. "||\textbf{A}||_p ||\textbf{x}||_p = \left( \max_{\mathbf{y} \not = \mathbf{ 0}} \frac{ || \textbf{A} \mathbf{y}||_p} { ||\mathbf{y}||_p} \right) ||\textbf{x}||_p \geq \frac{ || \textbf{A} \textbf{x}||_p} { ||\textbf{x}||_p} ||\textbf{x}||_p = || \textbf{A} \textbf{x}||_p.")

## Condition Number for ![\textbf{A}\textbf{x} = \textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%20%3D%20%5Ctextbf%7Bb%7D "\textbf{A}\textbf{x} = \textbf{b}")

Let us consider solving
![\textbf{A}\textbf{x}=\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%3D%5Ctextbf%7Bb%7D "\textbf{A}\textbf{x}=\textbf{b}").
Suppose we find approximate solution
![\textbf{x}\_a](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D_a "\textbf{x}_a").
The **relative forward error** is

![\frac{\|\|\textbf{x}-\textbf{x}\_a\|\|}{\|\|\textbf{x}\|\|}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B%7C%7C%5Ctextbf%7Bx%7D-%5Ctextbf%7Bx%7D_a%7C%7C%7D%7B%7C%7C%5Ctextbf%7Bx%7D%7C%7C%7D "\frac{||\textbf{x}-\textbf{x}_a||}{||\textbf{x}||}")

and the **relative backward error** is

![\frac{\|\|\textbf{A}\textbf{x}-\textbf{A}\textbf{x}\_a\|\|}{\|\|\textbf{A}\textbf{x}\|\|}=\frac{\|\|\textbf{b}-\textbf{A}\textbf{x}\_a\|\|}{\|\|\textbf{b}\|\|}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B%7C%7C%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D-%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D_a%7C%7C%7D%7B%7C%7C%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%7C%7C%7D%3D%5Cfrac%7B%7C%7C%5Ctextbf%7Bb%7D-%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D_a%7C%7C%7D%7B%7C%7C%5Ctextbf%7Bb%7D%7C%7C%7D. "\frac{||\textbf{A}\textbf{x}-\textbf{A}\textbf{x}_a||}{||\textbf{A}\textbf{x}||}=\frac{||\textbf{b}-\textbf{A}\textbf{x}_a||}{||\textbf{b}||}.")

We define **error magnification** as the ratio

![\frac{\text{relative forward error}}{\text{relative backward error}}=\frac{\frac{\|\|\textbf{x}-\textbf{x}\_a\|\|}{\|\|\textbf{x}\|\|}}{\frac{\|\|\textbf{b}-\textbf{A}\textbf{x}\_a\|\|}{\|\|\textbf{b}\|\|}}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B%5Ctext%7Brelative%20forward%20error%7D%7D%7B%5Ctext%7Brelative%20backward%20error%7D%7D%3D%5Cfrac%7B%5Cfrac%7B%7C%7C%5Ctextbf%7Bx%7D-%5Ctextbf%7Bx%7D_a%7C%7C%7D%7B%7C%7C%5Ctextbf%7Bx%7D%7C%7C%7D%7D%7B%5Cfrac%7B%7C%7C%5Ctextbf%7Bb%7D-%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D_a%7C%7C%7D%7B%7C%7C%5Ctextbf%7Bb%7D%7C%7C%7D%7D. "\frac{\text{relative forward error}}{\text{relative backward error}}=\frac{\frac{||\textbf{x}-\textbf{x}_a||}{||\textbf{x}||}}{\frac{||\textbf{b}-\textbf{A}\textbf{x}_a||}{||\textbf{b}||}}.")

The **condition number**
![\kappa(\textbf{A})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ckappa%28%5Ctextbf%7BA%7D%29 "\kappa(\textbf{A})")
is the largest possible error magnification (over all possible
![\textbf{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D "\textbf{x}")).
Or restated, it’s the worst possible ratio of relative forward error to
relative backward error.

Why do we care about this? It tells us whether we expect a small
residual to imply a small error in the solution or not. Let’s make this
concrete with an example.

Consider:

-   Let
    ![\textbf{A} = \begin{pmatrix}0.913 & 0.659 \\\\ 0.457 & 0.330 \end{pmatrix}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%3D%20%5Cbegin%7Bpmatrix%7D0.913%20%26%200.659%20%5C%5C%200.457%20%26%200.330%20%5Cend%7Bpmatrix%7D "\textbf{A} = \begin{pmatrix}0.913 & 0.659 \\ 0.457 & 0.330 \end{pmatrix}")
-   Then
    ![\kappa_2(A) = 1.25\times10^4](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ckappa_2%28A%29%20%3D%201.25%5Ctimes10%5E4 "\kappa_2(A) = 1.25\times10^4")
-   Let
    ![\textbf{b} = \begin{pmatrix} 0.254 \\\\ 0.127 \end{pmatrix}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bb%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%200.254%20%5C%5C%200.127%20%5Cend%7Bpmatrix%7D "\textbf{b} = \begin{pmatrix} 0.254 \\ 0.127 \end{pmatrix}")
-   Then
    ![\textbf{x} =(1,-1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D%20%3D%281%2C-1%29 "\textbf{x} =(1,-1)").
-   Consider two approximate solutions
    ![\textbf{x}\_{1,2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D_%7B1%2C2%7D "\textbf{x}_{1,2}")

![\begin{array}{lcl}
\textbf{x}\_1 = (-0.0827,0.5) && \textbf{x}\_2 = (0.999,-1.001) \\\\
\triangle \textbf{x}\_1 = (1.0827, -1.5) && \triangle \textbf{x}\_2 = (0.001,0.001) \\\\
\|\|\triangle \textbf{x}\_1 \|\| = 1.85&& \|\|\triangle \textbf{x}\_2\|\| = .0014\\\\
\|\|\triangle \textbf{x}\_1 \|\|/\|\|\textbf{x}\|\| = 1.308&& \|\|\triangle \textbf{x}\_2\|\|/\|\|\textbf{x}\|\| = .001\\\\
\textbf{b}\_1 = (0.2539949, 0.1272061) &&  \textbf{b}\_2 = (0.252428, 0.126213) \\\\
\triangle\textbf{b}\_1 = (0.0000051,- 0.0002061) &&  \triangle\textbf{b}\_2 = (0.001572, 0.000787) \\\\
\|\|\triangle \textbf{b}\_1 \|\| = 0.000206&& \|\|\triangle \textbf{b}\_2\|\| = .00176\\\\
\|\|\triangle \textbf{b}\_1 \|\|/\|\|\textbf{b}\|\| = 0.000726&& \|\|\triangle \textbf{b}\_2\|\|/\|\|\textbf{b}\|\| = .0062\\\\
mag = 1.8 \times 10^3 && mag = 1.6 \times 10^1
\end{array}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Barray%7D%7Blcl%7D%0A%5Ctextbf%7Bx%7D_1%20%3D%20%28-0.0827%2C0.5%29%20%26%26%20%5Ctextbf%7Bx%7D_2%20%3D%20%280.999%2C-1.001%29%20%5C%5C%0A%5Ctriangle%20%5Ctextbf%7Bx%7D_1%20%3D%20%281.0827%2C%20-1.5%29%20%26%26%20%5Ctriangle%20%5Ctextbf%7Bx%7D_2%20%3D%20%280.001%2C0.001%29%20%5C%5C%0A%7C%7C%5Ctriangle%20%5Ctextbf%7Bx%7D_1%20%7C%7C%20%3D%201.85%26%26%20%7C%7C%5Ctriangle%20%5Ctextbf%7Bx%7D_2%7C%7C%20%3D%20.0014%5C%5C%0A%7C%7C%5Ctriangle%20%5Ctextbf%7Bx%7D_1%20%7C%7C%2F%7C%7C%5Ctextbf%7Bx%7D%7C%7C%20%3D%201.308%26%26%20%7C%7C%5Ctriangle%20%5Ctextbf%7Bx%7D_2%7C%7C%2F%7C%7C%5Ctextbf%7Bx%7D%7C%7C%20%3D%20.001%5C%5C%0A%5Ctextbf%7Bb%7D_1%20%3D%20%280.2539949%2C%200.1272061%29%20%26%26%20%20%5Ctextbf%7Bb%7D_2%20%3D%20%280.252428%2C%200.126213%29%20%5C%5C%0A%5Ctriangle%5Ctextbf%7Bb%7D_1%20%3D%20%280.0000051%2C-%200.0002061%29%20%26%26%20%20%5Ctriangle%5Ctextbf%7Bb%7D_2%20%3D%20%280.001572%2C%200.000787%29%20%5C%5C%0A%7C%7C%5Ctriangle%20%5Ctextbf%7Bb%7D_1%20%7C%7C%20%3D%200.000206%26%26%20%7C%7C%5Ctriangle%20%5Ctextbf%7Bb%7D_2%7C%7C%20%3D%20.00176%5C%5C%0A%7C%7C%5Ctriangle%20%5Ctextbf%7Bb%7D_1%20%7C%7C%2F%7C%7C%5Ctextbf%7Bb%7D%7C%7C%20%3D%200.000726%26%26%20%7C%7C%5Ctriangle%20%5Ctextbf%7Bb%7D_2%7C%7C%2F%7C%7C%5Ctextbf%7Bb%7D%7C%7C%20%3D%20.0062%5C%5C%0Amag%20%3D%201.8%20%5Ctimes%2010%5E3%20%26%26%20mag%20%3D%201.6%20%5Ctimes%2010%5E1%0A%5Cend%7Barray%7D "\begin{array}{lcl}
\textbf{x}_1 = (-0.0827,0.5) && \textbf{x}_2 = (0.999,-1.001) \\
\triangle \textbf{x}_1 = (1.0827, -1.5) && \triangle \textbf{x}_2 = (0.001,0.001) \\
||\triangle \textbf{x}_1 || = 1.85&& ||\triangle \textbf{x}_2|| = .0014\\
||\triangle \textbf{x}_1 ||/||\textbf{x}|| = 1.308&& ||\triangle \textbf{x}_2||/||\textbf{x}|| = .001\\
\textbf{b}_1 = (0.2539949, 0.1272061) &&  \textbf{b}_2 = (0.252428, 0.126213) \\
\triangle\textbf{b}_1 = (0.0000051,- 0.0002061) &&  \triangle\textbf{b}_2 = (0.001572, 0.000787) \\
||\triangle \textbf{b}_1 || = 0.000206&& ||\triangle \textbf{b}_2|| = .00176\\
||\triangle \textbf{b}_1 ||/||\textbf{b}|| = 0.000726&& ||\triangle \textbf{b}_2||/||\textbf{b}|| = .0062\\
mag = 1.8 \times 10^3 && mag = 1.6 \times 10^1
\end{array}")

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
it’s the maximum possible error magnification. Computing
![\kappa](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ckappa "\kappa")
exactly using this definition is impossible because there are an
infinite number vectors one must consider
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
acting on.

Fortunately, there’s another way to calculate condition number:

![\kappa_p(\textbf{A})=\|\|\textbf{A}\|\|\_p\|\|\textbf{A}^{-1}\|\|\_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ckappa_p%28%5Ctextbf%7BA%7D%29%3D%7C%7C%5Ctextbf%7BA%7D%7C%7C_p%7C%7C%5Ctextbf%7BA%7D%5E%7B-1%7D%7C%7C_p "\kappa_p(\textbf{A})=||\textbf{A}||_p||\textbf{A}^{-1}||_p")

The derivation of this identity is about 10 to 20 lines of linear
algebra that I am happy to show you if you are interested. We can check
it numerically for now.

``` r
set.seed(123)
N <- 10
A <- matrix(runif(N^2), nrow = N)
ans1 <- norm(A, "2") * norm(solve(A), "2")
ans2 <- kappa(A, norm = "2", exact = TRUE)
```

# LU decomposition

## Big picture

When solving systems of linear equations, depending on the context,
solution by Gaussian elimination can be computationally costly.
Sometimes it is better to decompose (factor) the matrix. There are a
number of useful ways to do this but the one we will focus on is called
LU decomposision.

## Goals

-   Define and implement LU decomposition and state the potential
    advantages

## LU decomposition

Recall that when solving an
![n \times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%20n "n \times n")
system
![\textbf{A}\textbf{x} = \textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%20%3D%20%5Ctextbf%7Bb%7D "\textbf{A}\textbf{x} = \textbf{b}")
with Gaussian elimination, the elimination step is
![\mathcal{O}(n^3)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathcal%7BO%7D%28n%5E3%29 "\mathcal{O}(n^3)")
and back substitution is
![\mathcal{O}(n^2)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathcal%7BO%7D%28n%5E2%29 "\mathcal{O}(n^2)").
In some applications, it is necessary to solve

![\textbf{A}\textbf{x} = \textbf{b}\_1, \quad \textbf{A}\textbf{x} = \textbf{b}\_2, \quad \textbf{A}\textbf{x} = \textbf{b}\_3, \quad \ldots, \quad \textbf{A}\textbf{x} = \textbf{b}\_M](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%20%3D%20%5Ctextbf%7Bb%7D_1%2C%20%5Cquad%20%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%20%3D%20%5Ctextbf%7Bb%7D_2%2C%20%5Cquad%20%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%20%3D%20%5Ctextbf%7Bb%7D_3%2C%20%5Cquad%20%5Cldots%2C%20%5Cquad%20%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%20%3D%20%5Ctextbf%7Bb%7D_M "\textbf{A}\textbf{x} = \textbf{b}_1, \quad \textbf{A}\textbf{x} = \textbf{b}_2, \quad \textbf{A}\textbf{x} = \textbf{b}_3, \quad \ldots, \quad \textbf{A}\textbf{x} = \textbf{b}_M")

where
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
is the same each time and
![M](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;M "M")
is large.
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
itself needs the same row reductions each time. Only the augmented part
![\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bb%7D "\textbf{b}")
changes. It would be a waste of computation to run Gaussian elimination
![M](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;M "M")
times.

LU decomposition is a way of storing the Gaussian elimination steps in
matrix form so that they can be applied to many
![\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bb%7D "\textbf{b}").
We take
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
and decompose (or factorize) it as the product

![A  = 
\underbrace{\begin{array}{\|cccccc\|}
\hline
1 &&&&&\\\\
\ast & 1 &&&&\\\\
\ast & \ast  & 1 &&&\\\\
\ast & \ast & \ast & 1 &&\\\\
\ast & \ast & \ast & \ast & 1  &\\\\
\ast & \ast & \ast & \ast & \ast & 1 \\\\ \hline
\end{array}}\_{\mathbf{L}=\text{Lower unit triangular}} \\
\underbrace{\begin{array}{\|cccccc\|}
\hline
\ast & \ast & \ast & \ast & \ast & \ast \\\\ 
&\ast & \ast & \ast & \ast & \ast  \\\\
&&\ast & \ast & \ast & \ast \\\\
&&&\ast & \ast  & \ast \\\\
&&&&\ast & \ast \\\\
&&&&&\ast \\\\
\hline
\end{array}}\_{\mathbf{U}=\text{Upper triangular}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A%20%20%3D%20%0A%5Cunderbrace%7B%5Cbegin%7Barray%7D%7B%7Ccccccc%7C%7D%0A%5Chline%0A1%20%26%26%26%26%26%5C%5C%0A%5Cast%20%26%201%20%26%26%26%26%5C%5C%0A%5Cast%20%26%20%5Cast%20%20%26%201%20%26%26%26%5C%5C%0A%5Cast%20%26%20%5Cast%20%26%20%5Cast%20%26%201%20%26%26%5C%5C%0A%5Cast%20%26%20%5Cast%20%26%20%5Cast%20%26%20%5Cast%20%26%201%20%20%26%5C%5C%0A%5Cast%20%26%20%5Cast%20%26%20%5Cast%20%26%20%5Cast%20%26%20%5Cast%20%26%201%20%5C%5C%20%5Chline%0A%5Cend%7Barray%7D%7D_%7B%5Cmathbf%7BL%7D%3D%5Ctext%7BLower%20unit%20triangular%7D%7D%20%5C%0A%5Cunderbrace%7B%5Cbegin%7Barray%7D%7B%7Ccccccc%7C%7D%0A%5Chline%0A%5Cast%20%26%20%5Cast%20%26%20%5Cast%20%26%20%5Cast%20%26%20%5Cast%20%26%20%5Cast%20%5C%5C%20%0A%26%5Cast%20%26%20%5Cast%20%26%20%5Cast%20%26%20%5Cast%20%26%20%5Cast%20%20%5C%5C%0A%26%26%5Cast%20%26%20%5Cast%20%26%20%5Cast%20%26%20%5Cast%20%5C%5C%0A%26%26%26%5Cast%20%26%20%5Cast%20%20%26%20%5Cast%20%5C%5C%0A%26%26%26%26%5Cast%20%26%20%5Cast%20%5C%5C%0A%26%26%26%26%26%5Cast%20%5C%5C%0A%5Chline%0A%5Cend%7Barray%7D%7D_%7B%5Cmathbf%7BU%7D%3D%5Ctext%7BUpper%20triangular%7D%7D "A  = 
\underbrace{\begin{array}{|cccccc|}
\hline
1 &&&&&\\
\ast & 1 &&&&\\
\ast & \ast  & 1 &&&\\
\ast & \ast & \ast & 1 &&\\
\ast & \ast & \ast & \ast & 1  &\\
\ast & \ast & \ast & \ast & \ast & 1 \\ \hline
\end{array}}_{\mathbf{L}=\text{Lower unit triangular}} \
\underbrace{\begin{array}{|cccccc|}
\hline
\ast & \ast & \ast & \ast & \ast & \ast \\ 
&\ast & \ast & \ast & \ast & \ast  \\
&&\ast & \ast & \ast & \ast \\
&&&\ast & \ast  & \ast \\
&&&&\ast & \ast \\
&&&&&\ast \\
\hline
\end{array}}_{\mathbf{U}=\text{Upper triangular}}")

The
![\mathbf{L}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BL%7D "\mathbf{L}")
matrix encodes the multipliers used to eliminate elements during
Gaussian elimination and the
![\mathbf{U}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BU%7D "\mathbf{U}")
matrix is the result of the elimination process. Therefore, putting
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
into its LU factorization takes one application of Gaussian elimination,
or approximately
![\frac{2}{3} n^3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B2%7D%7B3%7D%20n%5E3 "\frac{2}{3} n^3")
operations. Solving
![\mathbf{LU}\textbf{x} = \textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BLU%7D%5Ctextbf%7Bx%7D%20%3D%20%5Ctextbf%7Bb%7D "\mathbf{LU}\textbf{x} = \textbf{b}")
requires 2 back substitutions, namely one to solve
![\mathbf{L}\mathbf{y} = \textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BL%7D%5Cmathbf%7By%7D%20%3D%20%5Ctextbf%7Bb%7D "\mathbf{L}\mathbf{y} = \textbf{b}")
for
![\mathbf{y}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7By%7D "\mathbf{y}")
and one to solve
![\mathbf{U}\textbf{x}=\mathbf{y}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BU%7D%5Ctextbf%7Bx%7D%3D%5Cmathbf%7By%7D "\mathbf{U}\textbf{x}=\mathbf{y}")
for
![\textbf{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D "\textbf{x}").
This takes
![2n^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2n%5E2 "2n^2")
operations. So, to solve
![\textbf{A}\textbf{x} = \textbf{b}\_1, \ldots, \textbf{A}\textbf{x} = \textbf{b}\_M](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%20%3D%20%5Ctextbf%7Bb%7D_1%2C%20%5Cldots%2C%20%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%20%3D%20%5Ctextbf%7Bb%7D_M "\textbf{A}\textbf{x} = \textbf{b}_1, \ldots, \textbf{A}\textbf{x} = \textbf{b}_M")
takes approximately
![\frac{2}{3}n^3 + 2 M n^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B2%7D%7B3%7Dn%5E3%20%2B%202%20M%20n%5E2 "\frac{2}{3}n^3 + 2 M n^2")
operations, in contrast to
![\frac{2}{3}Mn^3 + Mn^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B2%7D%7B3%7DMn%5E3%20%2B%20Mn%5E2 "\frac{2}{3}Mn^3 + Mn^2")
for Gaussian elimination.

The LU decomposition exists if and only if the upper-left sub-blocks
![\textbf{A}\_{1:k,1:k}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D_%7B1%3Ak%2C1%3Ak%7D "\textbf{A}_{1:k,1:k}")
are non-singular for all
![1\leq k \leq n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%5Cleq%20k%20%5Cleq%20n "1\leq k \leq n")
(not proven here). If the decomposition exists, it is unique.

Let’s do an example of how LU decomposition works. Take

![\textbf{A} = 
\begin{pmatrix}
1 & 3 & 1\\\\
1 & 1 & -1 \\\\
3 & 11 & 5
\end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%3D%20%0A%5Cbegin%7Bpmatrix%7D%0A1%20%26%203%20%26%201%5C%5C%0A1%20%26%201%20%26%20-1%20%5C%5C%0A3%20%26%2011%20%26%205%0A%5Cend%7Bpmatrix%7D. "\textbf{A} = 
\begin{pmatrix}
1 & 3 & 1\\
1 & 1 & -1 \\
3 & 11 & 5
\end{pmatrix}.")

Start by defining
![\mathbf{U} = \textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BU%7D%20%3D%20%5Ctextbf%7BA%7D "\mathbf{U} = \textbf{A}")
(it’s not upper triangluar yet, but I am still calling it
![\mathbf{U}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BU%7D "\mathbf{U}"))
and

![\mathbf{L} = \mathbf{I}\_3 = 
\begin{pmatrix}
1 & 0 & 0\\\\
0 & 1 & 0 \\\\
0 & 0 & 1
\end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BL%7D%20%3D%20%5Cmathbf%7BI%7D_3%20%3D%20%0A%5Cbegin%7Bpmatrix%7D%0A1%20%26%200%20%26%200%5C%5C%0A0%20%26%201%20%26%200%20%5C%5C%0A0%20%26%200%20%26%201%0A%5Cend%7Bpmatrix%7D. "\mathbf{L} = \mathbf{I}_3 = 
\begin{pmatrix}
1 & 0 & 0\\
0 & 1 & 0 \\
0 & 0 & 1
\end{pmatrix}.")

Apply
![II \leftarrow II - 1\cdot I](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;II%20%5Cleftarrow%20II%20-%201%5Ccdot%20I "II \leftarrow II - 1\cdot I")
and
![III \leftarrow III - 3\cdot I](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;III%20%5Cleftarrow%20III%20-%203%5Ccdot%20I "III \leftarrow III - 3\cdot I")
to
![\mathbf{U}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BU%7D "\mathbf{U}"),
so

![\mathbf{U} = \begin{pmatrix}
1 & 3 & 1  \\\\
0 & -2 & -2\\\\
0 & 2 & 2
\end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BU%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%0A1%20%26%203%20%26%201%20%20%5C%5C%0A0%20%26%20-2%20%26%20-2%5C%5C%0A0%20%26%202%20%26%202%0A%5Cend%7Bpmatrix%7D. "\mathbf{U} = \begin{pmatrix}
1 & 3 & 1  \\
0 & -2 & -2\\
0 & 2 & 2
\end{pmatrix}.")

Also, store the multupliers
![1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1 "1")
and
![3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;3 "3")
in their corresponding rows in the first column of
![\mathbf{L}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BL%7D "\mathbf{L}"),
so that

![\mathbf{L} =
\begin{pmatrix}
1 & 0 & 0\\\\
1 & 1 & 0 \\\\
3 & 0 & 1
\end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BL%7D%20%3D%0A%5Cbegin%7Bpmatrix%7D%0A1%20%26%200%20%26%200%5C%5C%0A1%20%26%201%20%26%200%20%5C%5C%0A3%20%26%200%20%26%201%0A%5Cend%7Bpmatrix%7D. "\mathbf{L} =
\begin{pmatrix}
1 & 0 & 0\\
1 & 1 & 0 \\
3 & 0 & 1
\end{pmatrix}.")

Now apply
![III \leftarrow III + II](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;III%20%5Cleftarrow%20III%20%2B%20II "III \leftarrow III + II"),
which yields

![\mathbf{U} = \begin{pmatrix}
1 & 3 & 1  \\\\
0 & -2 & -2\\\\
0 & 0 & 0
\end{pmatrix}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BU%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%0A1%20%26%203%20%26%201%20%20%5C%5C%0A0%20%26%20-2%20%26%20-2%5C%5C%0A0%20%26%200%20%26%200%0A%5Cend%7Bpmatrix%7D "\mathbf{U} = \begin{pmatrix}
1 & 3 & 1  \\
0 & -2 & -2\\
0 & 0 & 0
\end{pmatrix}")

and

![\mathbf{L} =
\begin{pmatrix}
1 & 0 & 0\\\\
1 & 1 & 0 \\\\
3 & -1 & 1
\end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BL%7D%20%3D%0A%5Cbegin%7Bpmatrix%7D%0A1%20%26%200%20%26%200%5C%5C%0A1%20%26%201%20%26%200%20%5C%5C%0A3%20%26%20-1%20%26%201%0A%5Cend%7Bpmatrix%7D. "\mathbf{L} =
\begin{pmatrix}
1 & 0 & 0\\
1 & 1 & 0 \\
3 & -1 & 1
\end{pmatrix}.")

Since
![\mathbf{U}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BU%7D "\mathbf{U}")
is in echelon form, we are done! We can check that our decomposition
worked.

``` r
A <- matrix(c(1,3,1,1,1,-1,3,11,5),nrow=3,byrow=TRUE)
L <- matrix(c(1,0,0,1,1,0,3,-1,1),nrow=3,byrow=TRUE)
U <- matrix(c(1,3,1,0,-2,-2,0,0,0),nrow=3,byrow=TRUE)
A - L%*%U
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

Let’s use Gaussian elimination with the `echelon` command and LU
decomposition with the to compare the time for solving a
![100 \times 100](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;100%20%5Ctimes%20100 "100 \times 100")
system for
![100](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;100 "100")
different right hand sides.

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

    ## [1] 6.33333333333

# Iterative methods for linear systems

## Big picture

So far, all of the methods you have seen for solving linear systems
provide exact solutions (excluding numerical error). They are called
direct methods. However, they all involve at least one step that is
![\mathcal{O}(n^3)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathcal%7BO%7D%28n%5E3%29 "\mathcal{O}(n^3)").
If this computational cost is prohibitive, consider a potentially faster
iterative method at the expense of giving up having an exact solution.

## Goals

-   Explain how fixed point iteration relates to solution of an equation
-   Derive and implement Jacobi’s method, and explain advantages and
    limitations
-   Use convergence criteria for Jacobi’s method
-   State other iterative solution methods

## Fixed point iteration

Sometimes you can solve a problem by a method called **fixed point
iteration** whereby you just keep plugging into an expression until the
output equals the input. For example, suppose you want to solve
![(x-3)(x+1)=x^2-2x-3=0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28x-3%29%28x%2B1%29%3Dx%5E2-2x-3%3D0 "(x-3)(x+1)=x^2-2x-3=0").
Pretend you don’t know where the roots are but you think there is one
near x = -2, so you start out with that guess. You also notice you can
write
![x^2 - 2x - 3 =0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%5E2%20-%202x%20-%203%20%3D0 "x^2 - 2x - 3 =0")
as
![x = 3/(x-2)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%20%3D%203%2F%28x-2%29 "x = 3/(x-2)").
So you define an iteration
![x\_{i+1} = 3/(x_i-2)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_%7Bi%2B1%7D%20%3D%203%2F%28x_i-2%29 "x_{i+1} = 3/(x_i-2)").

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

## Jacobi iteration

Let’s take this idea and apply it to solving
![\textbf{A}\textbf{x}=\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%3D%5Ctextbf%7Bb%7D "\textbf{A}\textbf{x}=\textbf{b}").
Let
![\textbf{A} = \mathbf{D} + \mathbf{R}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%3D%20%5Cmathbf%7BD%7D%20%2B%20%5Cmathbf%7BR%7D "\textbf{A} = \mathbf{D} + \mathbf{R}")
where
![\mathbf{D}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BD%7D "\mathbf{D}")
contains the diagonal elements of
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
and
![\mathbf{R}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BR%7D "\mathbf{R}")
contains everything else. Then we can write

![\begin{align}
\textbf{A}\textbf{x} &= \textbf{b} \\\\
(\mathbf{D}+\mathbf{R})\textbf{x} &= \textbf{b} \\\\
\mathbf{D}\textbf{x} + \mathbf{R} \textbf{x} &= \textbf{b}\\\\
\mathbf{D} \textbf{x} &= \textbf{b} - \mathbf{R} \textbf{x}\\\\
\textbf{x} &= \mathbf{D}^{-1} (\textbf{b}-\mathbf{R}\textbf{x})
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0A%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%20%26%3D%20%5Ctextbf%7Bb%7D%20%5C%5C%0A%28%5Cmathbf%7BD%7D%2B%5Cmathbf%7BR%7D%29%5Ctextbf%7Bx%7D%20%26%3D%20%5Ctextbf%7Bb%7D%20%5C%5C%0A%5Cmathbf%7BD%7D%5Ctextbf%7Bx%7D%20%2B%20%5Cmathbf%7BR%7D%20%5Ctextbf%7Bx%7D%20%26%3D%20%5Ctextbf%7Bb%7D%5C%5C%0A%5Cmathbf%7BD%7D%20%5Ctextbf%7Bx%7D%20%26%3D%20%5Ctextbf%7Bb%7D%20-%20%5Cmathbf%7BR%7D%20%5Ctextbf%7Bx%7D%5C%5C%0A%5Ctextbf%7Bx%7D%20%26%3D%20%5Cmathbf%7BD%7D%5E%7B-1%7D%20%28%5Ctextbf%7Bb%7D-%5Cmathbf%7BR%7D%5Ctextbf%7Bx%7D%29%0A%5Cend%7Balign%7D "\begin{align}
\textbf{A}\textbf{x} &= \textbf{b} \\
(\mathbf{D}+\mathbf{R})\textbf{x} &= \textbf{b} \\
\mathbf{D}\textbf{x} + \mathbf{R} \textbf{x} &= \textbf{b}\\
\mathbf{D} \textbf{x} &= \textbf{b} - \mathbf{R} \textbf{x}\\
\textbf{x} &= \mathbf{D}^{-1} (\textbf{b}-\mathbf{R}\textbf{x})
\end{align}")

We can consider this a fixed point iteration,

![\textbf{x}\_{i+1} = \mathbf{D}^{-1} (\textbf{b}-\mathbf{R}\textbf{x}\_i).](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D_%7Bi%2B1%7D%20%3D%20%5Cmathbf%7BD%7D%5E%7B-1%7D%20%28%5Ctextbf%7Bb%7D-%5Cmathbf%7BR%7D%5Ctextbf%7Bx%7D_i%29. "\textbf{x}_{i+1} = \mathbf{D}^{-1} (\textbf{b}-\mathbf{R}\textbf{x}_i).")

Here are some computational advantages of this method.

1.  The matrix
    ![\mathbf{D}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BD%7D "\mathbf{D}")
    is very cheap to invert, because it is a diagonal matrix. The
    inverse is simply the diagonal matrix with the reciprocals of the
    original elements.
2.  Each iteration is only
    ![\mathcal{O}(n^2)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathcal%7BO%7D%28n%5E2%29 "\mathcal{O}(n^2)"),
    and if the matrix is sparse, it is even cheaper.

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
    ## 2.41025641026

## Convergence of Jacobi’s method

There’s no reason at all to expect that Jacobi’s method converges.
There’s a really useful theorem that says it it will converge if and
only if the eigenvalues of
![-\mathbf{D}^{-1}\mathbf{R}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-%5Cmathbf%7BD%7D%5E%7B-1%7D%5Cmathbf%7BR%7D "-\mathbf{D}^{-1}\mathbf{R}")
are all less than one in magnitude. If this criterion is met, then the
closer to one the eigenvalues are in magnitude, the slower convergence
will be. This takes about half a page to prove and is a worthwhile
exercise to understand, so feel free to ask me for the proof.

However, calculating the eigenvalues of that matrix could be really hard
and costly! There’s a sufficient condition for convergence that is much
easier to check computationally, namely that
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
is **strictly diagonally dominant**. This means

![\|a\_{ii}\| \> \sum\_{j \not=i} \|a\_{ij}\| \quad\text{in each row $i$.}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7Ca_%7Bii%7D%7C%20%3E%20%5Csum_%7Bj%20%5Cnot%3Di%7D%20%7Ca_%7Bij%7D%7C%20%5Cquad%5Ctext%7Bin%20each%20row%20%24i%24.%7D "|a_{ii}| > \sum_{j \not=i} |a_{ij}| \quad\text{in each row $i$.}")

## Other iterative methods

There are other iterative solution methods for linear systems that all
are inspired by Jacobi’s method. Some of these include Gauss-Seidel
iteration and Successive Over-Relaxation. Details of these appear in
book and you are welcome to discuss them with me.

# Polynomial interpolation

## Big picture

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

-   Explain the advantage of using polynomials to describe data
-   Implement Vandermonde interpolation and explain its pros/cons
-   Implement Lagrange inteprolation and explain its pros/cons
-   Explain and implement data compression via interpolation

## Why polynomials?

Suppose we have incomplete data and we’d like to estimate a piece of
information that we don’t have.

``` r
# I created some mystery data and hid it from you
plot(x,y)
```

![](coursenotes_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

We might try to use polynomials to describe the data, and then glean
information from the polynomial. Polynomials are convenient for several
reasons.

First, for any function that is defined and continuous on a closed,
bounded interval, there exists a polynomial that is as “close” to the
given function as desired. Stated more precisely, this is the
**Weierstrass Approximation Theorem**. Suppose
![f](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f "f")
is defined and continuous on
![\[a,b\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ba%2Cb%5D "[a,b]").
For each
![\epsilon \> 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon%20%3E%200 "\epsilon > 0"),
there exists a polynomial
![P(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;P%28x%29 "P(x)")
with the property that
![\|f(x)-P(x)\| \< \epsilon](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7Cf%28x%29-P%28x%29%7C%20%3C%20%5Cepsilon "|f(x)-P(x)| < \epsilon")
for all
![x \in \[a,b\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%20%5Cin%20%5Ba%2Cb%5D "x \in [a,b]").

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
Consider
![f(x) = 1/x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28x%29%20%3D%201%2Fx "f(x) = 1/x")
on the interval
![\[1,3\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5B1%2C3%5D "[1,3]").
The second degree Taylor polynomal through
![x=1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%3D1 "x=1")
(as an example) is
![T(x)=3-3x+x^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T%28x%29%3D3-3x%2Bx%5E2 "T(x)=3-3x+x^2").
The interpolating polynomial using three points from sampled equally
across the interval, namely
![(1,1),(2,1/2),(3,1/3)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%281%2C1%29%2C%282%2C1%2F2%29%2C%283%2C1%2F3%29 "(1,1),(2,1/2),(3,1/3)"),
is
![11/6 - x + (1/6)x^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;11%2F6%20-%20x%20%2B%20%281%2F6%29x%5E2 "11/6 - x + (1/6)x^2").

``` r
f <- function(x){1/x}
tee <- function(x){3-3*x+x^2}
p <- function(x){11/6 - x + 1/6*x^2}
xdata <- seq(from=1,to=3,length=3)
x <- seq(from=1,to=3,length=200)
plot(x,f(x),type="l",lwd=2)
lines(x,tee(x),col="red",lwd=2)
points(xdata,f(xdata),cex=2)
lines(x,p(x),col="green",lwd=2)
```

![](coursenotes_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

A few important points to know about the interpolating polynomial
through
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
points with distinct
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
coordinates include:

-   It exists
-   It is unique
-   It has degree at most
    ![n-1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n-1 "n-1")

Let’s examine some different ways to compute the interpolating
polynomial.

## Vandermonde matrix

I also call this method of interpolation “brute force.” Let’s start with
an example. Suppose we have three data points
![(x_1,y_1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28x_1%2Cy_1%29 "(x_1,y_1)"),
![(x_2,y_2)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28x_2%2Cy_2%29 "(x_2,y_2)"),
![(x_3,y_3)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28x_3%2Cy_3%29 "(x_3,y_3)").
The interpolating polynomial is
![P(x) = c_0 + c_1 x + c_2 x^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;P%28x%29%20%3D%20c_0%20%2B%20c_1%20x%20%2B%20c_2%20x%5E2 "P(x) = c_0 + c_1 x + c_2 x^2"),
where we need to determine the coefficients. Let’s determine them simply
by plugging in. We have the equations

![\begin{align\*}
c_0 + c_1 x_1 + c_2 x_1^2 &= y_1 \\\\
c_0 + c_1 x_2 + c_2 x_2^2 &= y_2 \\\\
c_0 + c_1 x_3 + c_2 x_3^2 &= y_3
\end{align\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%2A%7D%0Ac_0%20%2B%20c_1%20x_1%20%2B%20c_2%20x_1%5E2%20%26%3D%20y_1%20%5C%5C%0Ac_0%20%2B%20c_1%20x_2%20%2B%20c_2%20x_2%5E2%20%26%3D%20y_2%20%5C%5C%0Ac_0%20%2B%20c_1%20x_3%20%2B%20c_2%20x_3%5E2%20%26%3D%20y_3%0A%5Cend%7Balign%2A%7D "\begin{align*}
c_0 + c_1 x_1 + c_2 x_1^2 &= y_1 \\
c_0 + c_1 x_2 + c_2 x_2^2 &= y_2 \\
c_0 + c_1 x_3 + c_2 x_3^2 &= y_3
\end{align*}")

We can write this in matrix form as

![\begin{pmatrix}
1 & x_1 & x_1^2 \\\\
1 & x_2 & x_2^2 \\\\
1 & x_3 & x_3^2
\end{pmatrix}
\begin{pmatrix}
c_0 \\\\ c_1 \\\\ c_2
\end{pmatrix}
=
\begin{pmatrix}
y_1 \\\\ y_2 \\\\ y_3
\end{pmatrix}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Bpmatrix%7D%0A1%20%26%20x_1%20%26%20x_1%5E2%20%5C%5C%0A1%20%26%20x_2%20%26%20x_2%5E2%20%5C%5C%0A1%20%26%20x_3%20%26%20x_3%5E2%0A%5Cend%7Bpmatrix%7D%0A%5Cbegin%7Bpmatrix%7D%0Ac_0%20%5C%5C%20c_1%20%5C%5C%20c_2%0A%5Cend%7Bpmatrix%7D%0A%3D%0A%5Cbegin%7Bpmatrix%7D%0Ay_1%20%5C%5C%20y_2%20%5C%5C%20y_3%0A%5Cend%7Bpmatrix%7D "\begin{pmatrix}
1 & x_1 & x_1^2 \\
1 & x_2 & x_2^2 \\
1 & x_3 & x_3^2
\end{pmatrix}
\begin{pmatrix}
c_0 \\ c_1 \\ c_2
\end{pmatrix}
=
\begin{pmatrix}
y_1 \\ y_2 \\ y_3
\end{pmatrix}")

and solve the linear system to find the coefficients
![c_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c_i "c_i").
More generally, for
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
data points, the problem is

![\begin{pmatrix}
1 & x_1 & x_1^2 &  & x_1^{n-1} \\\\ 
1 & x_2 & x_2^2 & \cdots & x_2^{n-1} \\\\ 
1 & x_3 & x_3^2 &  & x_3^{n-1} \\\\ 
&  \vdots &  &  & \vdots \\\\ 
1 & x_n & x_n^2 & \cdots & x_n^{n-1} \\\\ 
\end{pmatrix}
\begin{pmatrix} c_0 \\\\ c_1 \\\\ c_2 \\\\ \vdots \\\\ c\_{n-1} \end{pmatrix}
= 
\begin{pmatrix} y_1 \\\\ y_2 \\\\ y_3 \\\\ \vdots \\\\ y_n \end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Bpmatrix%7D%0A1%20%26%20x_1%20%26%20x_1%5E2%20%26%20%20%26%20x_1%5E%7Bn-1%7D%20%5C%5C%20%0A1%20%26%20x_2%20%26%20x_2%5E2%20%26%20%5Ccdots%20%26%20x_2%5E%7Bn-1%7D%20%5C%5C%20%0A1%20%26%20x_3%20%26%20x_3%5E2%20%26%20%20%26%20x_3%5E%7Bn-1%7D%20%5C%5C%20%0A%26%20%20%5Cvdots%20%26%20%20%26%20%20%26%20%5Cvdots%20%5C%5C%20%0A1%20%26%20x_n%20%26%20x_n%5E2%20%26%20%5Ccdots%20%26%20x_n%5E%7Bn-1%7D%20%5C%5C%20%0A%5Cend%7Bpmatrix%7D%0A%5Cbegin%7Bpmatrix%7D%20c_0%20%5C%5C%20c_1%20%5C%5C%20c_2%20%5C%5C%20%5Cvdots%20%5C%5C%20c_%7Bn-1%7D%20%5Cend%7Bpmatrix%7D%0A%3D%20%0A%5Cbegin%7Bpmatrix%7D%20y_1%20%5C%5C%20y_2%20%5C%5C%20y_3%20%5C%5C%20%5Cvdots%20%5C%5C%20y_n%20%5Cend%7Bpmatrix%7D. "\begin{pmatrix}
1 & x_1 & x_1^2 &  & x_1^{n-1} \\ 
1 & x_2 & x_2^2 & \cdots & x_2^{n-1} \\ 
1 & x_3 & x_3^2 &  & x_3^{n-1} \\ 
&  \vdots &  &  & \vdots \\ 
1 & x_n & x_n^2 & \cdots & x_n^{n-1} \\ 
\end{pmatrix}
\begin{pmatrix} c_0 \\ c_1 \\ c_2 \\ \vdots \\ c_{n-1} \end{pmatrix}
= 
\begin{pmatrix} y_1 \\ y_2 \\ y_3 \\ \vdots \\ y_n \end{pmatrix}.")

It can be proven that if the
![x_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_i "x_i")
are distinct, the matrix has nonzero determinant, and hence the system
is solvable with a unique solution.

In solving problems using the Vandermonde matrix, we are using a basis
for the interpolating polynomial that is
![\\{1,x,x^2,\ldots\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7B1%2Cx%2Cx%5E2%2C%5Cldots%5C%7D "\{1,x,x^2,\ldots\}").
This seems very natural since this is how we usually think of
polynomials! The problems is that to find the coefficients in this
basis, we have to solve a linear problem whose matrix is very
ill-conditioned. Let’s see what happens if we sample more and more
points from our function and construct the interpolating polynomial.
We’ll look at
![\kappa](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ckappa "\kappa")
for the Vandermonde matrix.

``` r
Vandermonde <- function(x){
  n <- length(x)
  V <- outer(x, 0:(n-1), "^")
  return(V)
}
nvals <- 2^(1:8)
kappavals <- NULL
for (n in nvals){
  x <- seq(from=1,to=3,length=n)
  kappavals <- c(kappavals,kappa(Vandermonde(x)))
}
kable(cbind(nvals,kappavals))
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
7.50000000000e+00
</td>
</tr>
<tr>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
2.49186893872e+03
</td>
</tr>
<tr>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
5.75363020789e+08
</td>
</tr>
<tr>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
4.55564575643e+19
</td>
</tr>
<tr>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
4.91326454030e+28
</td>
</tr>
<tr>
<td style="text-align:right;">
64
</td>
<td style="text-align:right;">
5.17340173411e+44
</td>
</tr>
<tr>
<td style="text-align:right;">
128
</td>
<td style="text-align:right;">
3.07806605038e+76
</td>
</tr>
<tr>
<td style="text-align:right;">
256
</td>
<td style="text-align:right;">
5.91626505144e+136
</td>
</tr>
</tbody>
</table>

## Lagrange interpolating polynomial

Let’s seek an alternative method that gets around these issues. There’s
actually a way to simply write down the interpolating polynomial without
doing any calculation of coefficients at all. Let’s do an example.

For points
![(1,10), (2, 6), (3,4), (4,10)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%281%2C10%29%2C%20%282%2C%206%29%2C%20%283%2C4%29%2C%20%284%2C10%29 "(1,10), (2, 6), (3,4), (4,10)"),
consider the polynomial

![\begin{align\*}
p(x) &= 10 \frac{(x-2)(x-3)(x-4)}{(1-2)(1-3)(1-4)}
+ 6 \frac{(x-1)(x-3)(x-4)}{(2-1)(2-3)(2-4)} \\\\
&+ 4 \frac{(x-1)(x-2)(x-4)}{(3-1)(3-2)(3-4)} + 10 \frac{(x-1)(x-2)(x-3)}{(4-1)(4-2)(4-3)}
\end{align\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%2A%7D%0Ap%28x%29%20%26%3D%2010%20%5Cfrac%7B%28x-2%29%28x-3%29%28x-4%29%7D%7B%281-2%29%281-3%29%281-4%29%7D%0A%2B%206%20%5Cfrac%7B%28x-1%29%28x-3%29%28x-4%29%7D%7B%282-1%29%282-3%29%282-4%29%7D%20%5C%5C%0A%26%2B%204%20%5Cfrac%7B%28x-1%29%28x-2%29%28x-4%29%7D%7B%283-1%29%283-2%29%283-4%29%7D%20%2B%2010%20%5Cfrac%7B%28x-1%29%28x-2%29%28x-3%29%7D%7B%284-1%29%284-2%29%284-3%29%7D%0A%5Cend%7Balign%2A%7D "\begin{align*}
p(x) &= 10 \frac{(x-2)(x-3)(x-4)}{(1-2)(1-3)(1-4)}
+ 6 \frac{(x-1)(x-3)(x-4)}{(2-1)(2-3)(2-4)} \\
&+ 4 \frac{(x-1)(x-2)(x-4)}{(3-1)(3-2)(3-4)} + 10 \frac{(x-1)(x-2)(x-3)}{(4-1)(4-2)(4-3)}
\end{align*}")

First, notice the pattern of how each term is constructed: there’s a
coefficient times a polynomial. The coefficient is the
![y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y "y")-value
of an interpolation point. The polynomial is one that is equal to
![1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1 "1")
at that the
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")-coordinate
of that point, and
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0")
at the
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")-coordinate
of every other point. So now we ask…

-   Is the expression a polynomial?
-   Is it degree at most
    ![n-1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n-1 "n-1")?
-   Does it pass through each data point?

Since the answer to all these questions is yes, it is the interpolating
polynomial.

If we expand out our Lagrange polynomial we find
![p(x) = x^3 - 5 x^2 + 4 x+10](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%28x%29%20%3D%20x%5E3%20-%205%20x%5E2%20%2B%204%20x%2B10 "p(x) = x^3 - 5 x^2 + 4 x+10").
We can use `R` to verify this result by solving the Vandermonde problem.

``` r
x <- c(1,2,3,4)
y <- c(10,6,4,10)
c <- solve(Vandermonde(x),y)
c
```

    ## [1] 10  4 -5  1

Following the pattern we established above, the Langrange polynomial for
points
![(x_1,y_1),\ldots,(x_n,y_n)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28x_1%2Cy_1%29%2C%5Cldots%2C%28x_n%2Cy_n%29 "(x_1,y_1),\ldots,(x_n,y_n)")
is

![p(x) = \sum\_{i = 1}^n y_i \prod\_{j \not = i} \frac{(x - x_j)}{(x_i - x_j)}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%28x%29%20%3D%20%5Csum_%7Bi%20%3D%201%7D%5En%20y_i%20%5Cprod_%7Bj%20%5Cnot%20%3D%20i%7D%20%5Cfrac%7B%28x%20-%20x_j%29%7D%7B%28x_i%20-%20x_j%29%7D. "p(x) = \sum_{i = 1}^n y_i \prod_{j \not = i} \frac{(x - x_j)}{(x_i - x_j)}.")

The advantage of this method is that it doesn’t require any numerical
solution… just evaluation. Let’s try a comparison: Vandermonde
vs. Lagrange. Here, I’ll use R’s `baryalg` function. This function has
some strengths and some weaknesses.

``` r
set.seed(123)
n <- 10
x <- 1:n
y <- runif(n)
x0 <- 1.5
c <- solve(vander(x),y)
horner(c, x0)
```

    ## $y
    ## [1] -0.157060132098
    ## 
    ## $dy
    ## [1] 2.41438899384

``` r
barylag(x,y,x0)
```

    ## [1] -0.157060132133

We can also do a speed comparison test.

``` r
set.seed(123)
numTrials <- 1000
n <- 10
x <- 1:n
t1 <- system.time(
  for (i in 1:numTrials){
    y <- runif(n)
    c <- echelon(vander(x),y)[, n+1]
    horner(c,x0)
  }
)[3]
t2 <- system.time(
  for (i in 1:numTrials){
    y <- runif(n)
    barylag(x,y,x0)
  }
)[3]
as.numeric(t1/t2)
```

    ## [1] 24.3461538462

## Data compression

One of the powerful things interpolation can do is compress data. Let’s
do an example. Suppose that we need to know values for the function
![\sin x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csin%20x "\sin x").
A computer doesn’t magically know this function, so it has to have some
way to compute/evaluate it. One option would be to store a giant look-up
table. There are an infinite number of numbers to store, though, even
for the interval
![\[0,2\pi)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5B0%2C2%5Cpi%29 "[0,2\pi)").
Another option is to fit a polynomial based on a finite number of
points, store the coefficients, and evaluate the polynomial as needed.
Let’s do this using 5 points to begin with. We’ll write a function that
takes a specified number of points, samples them from the function,
constructs the interpolating polynomial, plots the function and the
polynomial, and calculates the maximum error. The inputs are your
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
data, your
![y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y "y")
data, and the
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
values at which you’d like interpolated values.

``` r
interperror <- function(n,plotflag=FALSE){
  x <- seq(from=0,to=2*pi,length=n)
  y <- sin(x)
  xexact <- seq(from=0,to=2*pi,length=1000)
  yexact <- sin(xexact)
  yinterp <- barylag(x,y,xexact)
  if (plotflag==TRUE){
    plot(xexact,yexact,type="l")
    points(x,y)
    lines(xexact,yinterp,col="blue")
  }
  error <- max(abs(yexact-yinterp))
  return(error)
}
interperror(5,plotflag=TRUE)
```

![](coursenotes_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

    ## [1] 0.180757796555

Not bad for just 5 points. Let’s examine how the error changes as a
function of
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n").

``` r
nvec <- 2:20
errorvec <- NULL
for (n in nvec){
  errorvec <- c(errorvec,interperror(n))
}
orderofmag <- round(log10(errorvec))
plot(nvec,orderofmag)
```

![](coursenotes_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

This means that we can represent the sine function with
![10^{-13}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;10%5E%7B-13%7D "10^{-13}")
error using only 20 pieces of information, instead of storing a huge
lookup table.

# Interpolation error and Chebyshev interpolation

## Big picture

We have been thinking about data and talking about polynomial
interpolation as a way of estimating it, compressing it, and
representing it conveniently to do other mathematical operations on it.
So far, we have only examined error numerically. Now it is time to look
at the error in more detail, including finding out when it is
potentially large and thinking about how we can reduce it.

## Goals

-   Describe and recognize Runge’s phenomenon
-   State the error term for polynomial interpolation and bound it
-   Explain the advantages of Chebyshev integration and implement the
    technique
-   Compare approaches to interpolation

## Runge’s phenomenon

Before doing polynomial interpolation, let’s start out with an example
that is about Taylor polynomials. We’ll construct Taylor polynomials of
increasing degree to estimate the function
![f(x) = \\\cos(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28x%29%20%3D%20%5C%5Ccos%28x%29 "f(x) = \\cos(x)")
around the point
![x_0=0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_0%3D0 "x_0=0")
on
![\[0,2\\\pi\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5B0%2C2%5C%5Cpi%5D "[0,2\\pi]").
The
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")th
degree taylor polynomial is

![\\\sum\\\_{i=1}^n (-1)^{i/2} \\\frac{x^i}{i!}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%5Csum%5C_%7Bi%3D1%7D%5En%20%28-1%29%5E%7Bi%2F2%7D%20%5C%5Cfrac%7Bx%5Ei%7D%7Bi%21%7D. "\\sum\_{i=1}^n (-1)^{i/2} \\frac{x^i}{i!}.")

``` r
mytaylor1 <- function(x,n){
  ans <- 0
  for (i in seq(from=0,to=n,by=2)){
    ans <- ans + (-1)^(i/2)*x^i/factorial(i)
  }
  return(ans)
}
x <- seq(from=0,to=2*pi,length=1000)
plot(x,cos(x),type="l",col="red",lwd=5,xlim=c(0,2*pi),ylim=c(-1.5,1.5),xlab="x",ylab="y")
for (n in seq(from=0,to=14,by=2)){
  lines(x,mytaylor1(x,n))
}
```

![](coursenotes_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

So more terms are better, right? Let’s try again with the function
![f(x) = 1/x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28x%29%20%3D%201%2Fx "f(x) = 1/x")
around the point
![x_0=1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_0%3D1 "x_0=1").
The
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")th
degree Taylor polynomial is

![\\\sum\\\_{i=1}^n (-1)^i (x-1)^i.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%5Csum%5C_%7Bi%3D1%7D%5En%20%28-1%29%5Ei%20%28x-1%29%5Ei. "\\sum\_{i=1}^n (-1)^i (x-1)^i.")

``` r
mytaylor2 <- function(x,n){
  ans <- 0
  for (i in seq(from=0,to=n,by=1)){
    ans <- ans + (-1)^(i)*(x-1)^(i)
  }
  return(ans)
}
x <- seq(from=1,to=2.5,length=1000)
plot(x,1/x,type="l",col="red",lwd=5,xlim=c(1,2.5),ylim=c(0,1.5),xlab="x",ylab="y")
for (n in seq(from=0,to=40,by=4)){
  lines(x,mytaylor2(x,n))
}
```

![](coursenotes_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

Oh! I guess that more isn’t always better.

Now let’s think now about interpolating polynomials. Let’s consider
![\\\cos(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%5Ccos%28x%29 "\\cos(x)")
with
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
equally sampled points across
![\[0,2\\\pi\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5B0%2C2%5C%5Cpi%5D "[0,2\\pi]")
for different values of
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n").

``` r
x <- seq(from=0,to=2*pi,length=1000)
y <- cos(x)
plot(x,y,type="l",col="red",lwd=5,xlim=c(0,2*pi),ylim=c(-1.1,1.1))
nvec <- 2:20
error <- NULL
for (n in nvec){
  xdata <- seq(from=0,to=2*pi,length=n)
  ydata <- cos(xdata)
  yinterp <- barylag(xdata,ydata,x)
  lines(x,yinterp)
  error <- c(error,max(abs(y-yinterp)))
}
```

![](coursenotes_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

``` r
plot(nvec,log10(error),xlab="n",ylab="log10 of error")
```

![](coursenotes_files/figure-gfm/unnamed-chunk-32-2.png)<!-- -->

Looks good! Let’s try again with a different function,
![f(x) = (1+x^2)^{-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28x%29%20%3D%20%281%2Bx%5E2%29%5E%7B-1%7D "f(x) = (1+x^2)^{-1}")
on
![\[-1,1\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5B-1%2C1%5D "[-1,1]").

``` r
x <- seq(from=-5,to=5,length=1000)
y <- 1/(1+x^2)
plot(x,y,type="l",col="red",lwd=5,xlim=c(-5,5),ylim=c(-3,3))
nvec <- seq(from=2,to=32,by=6)
error <- NULL
for (n in nvec){
  xdata <- seq(from=-5,to=5,length=n)
  ydata <- 1/(1+xdata^2)
  yinterp <- barylag(xdata,ydata,x)
  lines(x,yinterp)
  error <- c(error,max(abs(y-yinterp)))
}
```

![](coursenotes_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

``` r
plot(nvec,log10(error),xlab="n",ylab="log10 of error")
```

![](coursenotes_files/figure-gfm/unnamed-chunk-33-2.png)<!-- -->

Not good! The error goes up as we take more and more points.
Equally-spaced nodes are very natural in many applications (scientific
measurements, audio/video signals, etc.). But sadly, it turns out that
in some circumstances, inteprolating with polynomials through
equally-spaced nodes leads to very undesirable oscillations like those
above, called \*\*Runge’s phenomenon\*\*.

## Interpolation error

You might remember that for an
![n-1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n-1 "n-1")st
degree Taylor series (that is,
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
coefficients) of
![f(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28x%29 "f(x)")
centered around
![x=x_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%3Dx_0 "x=x_0"),
the magnitude of the error term is

![\\\left\\\|\\\frac{f^{(n)}(c)}{n!}(x-x_0)^{n}\\\right\\\|](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%5Cleft%5C%7C%5C%5Cfrac%7Bf%5E%7B%28n%29%7D%28c%29%7D%7Bn%21%7D%28x-x_0%29%5E%7Bn%7D%5C%5Cright%5C%7C "\\left\|\\frac{f^{(n)}(c)}{n!}(x-x_0)^{n}\\right\|")

where
![c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c "c")
is a number between
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
and
![x_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_0 "x_0").
For an interpolating polynomial constructed from
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
points on an interval
![\[x_1,x_n\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Bx_1%2Cx_n%5D "[x_1,x_n]"),
the interpolation error is

![\\\left\\\|\\\frac{f^{(n)}(c)}{n!}(x-x_1)(x-x_2)\\\cdots (x-x_n)\\\right\\\|](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%5Cleft%5C%7C%5C%5Cfrac%7Bf%5E%7B%28n%29%7D%28c%29%7D%7Bn%21%7D%28x-x_1%29%28x-x_2%29%5C%5Ccdots%20%28x-x_n%29%5C%5Cright%5C%7C "\\left\|\\frac{f^{(n)}(c)}{n!}(x-x_1)(x-x_2)\\cdots (x-x_n)\\right\|")

for some
![c \\\in \[x_1,x_n\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c%20%5C%5Cin%20%5Bx_1%2Cx_n%5D "c \\in [x_1,x_n]").
If you are interested in the proof, it appears in your textbook and/or I
am happy to go over it with you in office hours. Honestly, though, I
find the proof to be totally nonintuitive and unenlightening.

Let’s use the error term to do an example, where we again ask the
question HOW DO WE EVEN KNOW ANYTHING? For instance, suppose you want to
calculate the value of an exponential function. [This turns out to be
hard](https://math.stackexchange.com/questions/1239352/how-do-pocket-calculators-calculate-exponents)!
Let’s restrict the problem a bit. Suppose we want to calculate
![\\\mathrm{e}^{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%5Cmathrm%7Be%7D%5E%7Bx%7D "\\mathrm{e}^{x}")
on
![\[0,1\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5B0%2C1%5D "[0,1]")
with 5 digits of accuracy. One way to do this would be to choose
equally-spaced points on
![\[0,1\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5B0%2C1%5D "[0,1]")
and make a table of
![\\\mathrm{e}^{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%5Cmathrm%7Be%7D%5E%7Bx%7D "\\mathrm{e}^{x}")
for those values using some very accurate method. Then, we let the user
choose a value of
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x").
We find the two values that surround it in our table, construct the
linear interpolating polynomial between them, and use it to estimate our
answer. For instance, if we chose
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
values of
![0.1, 0.2, \\\ldots, 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.1%2C%200.2%2C%20%5C%5Cldots%2C%201 "0.1, 0.2, \\ldots, 1")
and the user input
![x=0.3268](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%3D0.3268 "x=0.3268"),
we would construct the inteprolating polynomial through
![(0.3,\\\mathrm{e}^{0.3})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%280.3%2C%5C%5Cmathrm%7Be%7D%5E%7B0.3%7D%29 "(0.3,\\mathrm{e}^{0.3})")
and
![(0.4,\\\mathrm{e}^{0.4})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%280.4%2C%5C%5Cmathrm%7Be%7D%5E%7B0.4%7D%29 "(0.4,\\mathrm{e}^{0.4})")
and then plug in our
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x").
Let’s ask the question: to achieve 5-digit accuracy with this method,
how many points must we take on
![\[0,1\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5B0%2C1%5D "[0,1]")?
(Note: crucially, this is not at all the same as creating a single
interpolating polynomial that goes through all the points.) We have

![\\\begin{align\\\*}
\\\|f(x) - P_1(x)\\\| &= \\\left\\\|\\\frac{f^{(n)}(c)}{n!}(x-x_1)(x-x_2)\\\cdots (x-x_n)\\\right\\\|\\
&= \\\left\\\|\\\frac{\\\left(\\\mathrm{e}^{x}\\\right)^{''}\\\|\\\_c}{2!}(x-x_1)(x-x_2)\\\right\\\| \\
& = \\\left\\\|\\\frac{\\\mathrm{e}^{c}}{2} (x-x_1)(x-x_1-h)\\\right\\\|
\\\end{align\\\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%5Cbegin%7Balign%5C%2A%7D%0A%5C%7Cf%28x%29%20-%20P_1%28x%29%5C%7C%20%26%3D%20%5C%5Cleft%5C%7C%5C%5Cfrac%7Bf%5E%7B%28n%29%7D%28c%29%7D%7Bn%21%7D%28x-x_1%29%28x-x_2%29%5C%5Ccdots%20%28x-x_n%29%5C%5Cright%5C%7C%5C%0A%26%3D%20%5C%5Cleft%5C%7C%5C%5Cfrac%7B%5C%5Cleft%28%5C%5Cmathrm%7Be%7D%5E%7Bx%7D%5C%5Cright%29%5E%7B%27%27%7D%5C%7C%5C_c%7D%7B2%21%7D%28x-x_1%29%28x-x_2%29%5C%5Cright%5C%7C%20%5C%0A%26%20%3D%20%5C%5Cleft%5C%7C%5C%5Cfrac%7B%5C%5Cmathrm%7Be%7D%5E%7Bc%7D%7D%7B2%7D%20%28x-x_1%29%28x-x_1-h%29%5C%5Cright%5C%7C%0A%5C%5Cend%7Balign%5C%2A%7D "\\begin{align\*}
\|f(x) - P_1(x)\| &= \\left\|\\frac{f^{(n)}(c)}{n!}(x-x_1)(x-x_2)\\cdots (x-x_n)\\right\|\
&= \\left\|\\frac{\\left(\\mathrm{e}^{x}\\right)^{''}\|\_c}{2!}(x-x_1)(x-x_2)\\right\| \
& = \\left\|\\frac{\\mathrm{e}^{c}}{2} (x-x_1)(x-x_1-h)\\right\|
\\end{align\*}")

Here, I’ve called the spacing between the two points
![h](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;h "h"),
that is
![x_2 = x_1 + h](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_2%20%3D%20x_1%20%2B%20h "x_2 = x_1 + h").
Now, we can ask what is the worst (biggest) that the term
![\left\|(x-x_1)(x-x_1-h)\right\|](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%7C%28x-x_1%29%28x-x_1-h%29%5Cright%7C "\left|(x-x_1)(x-x_1-h)\right|")
can be, and via calculus, we can show it is
![h^2/4](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;h%5E2%2F4 "h^2/4").
Also, on our interval of interest, we know that
![\mathrm{e}^{c} \leq \mathrm{e}^{1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathrm%7Be%7D%5E%7Bc%7D%20%5Cleq%20%5Cmathrm%7Be%7D%5E%7B1%7D "\mathrm{e}^{c} \leq \mathrm{e}^{1}").
Therefore, we can write

![\\\|f(x) - P_1(x)\\\| \leq \\\frac{\\\mathrm{e}h^2}{8}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7Cf%28x%29%20-%20P_1%28x%29%5C%7C%20%5Cleq%20%5C%5Cfrac%7B%5C%5Cmathrm%7Be%7Dh%5E2%7D%7B8%7D. "\|f(x) - P_1(x)\| \leq \\frac{\\mathrm{e}h^2}{8}.")

For our desired level of accuracy, we need

![\frac{\mathrm{e}h^2}{8} \leq 0.5 \times 10^{-5}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B%5Cmathrm%7Be%7Dh%5E2%7D%7B8%7D%20%5Cleq%200.5%20%5Ctimes%2010%5E%7B-5%7D. "\frac{\mathrm{e}h^2}{8} \leq 0.5 \times 10^{-5}.")

Solving for
![h](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;h "h"),
we find
![h \lessapprox 3.8 \times 10^{-3}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;h%20%5Clessapprox%203.8%20%5Ctimes%2010%5E%7B-3%7D "h \lessapprox 3.8 \times 10^{-3}").
That is the space between points, so on
![\[0,1\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5B0%2C1%5D "[0,1]")
this corresponds to just over 260 points.

## Chebyshev nodes

Look at the error expression again. There are three parts of it.
Factorial is the same for every function. If the derivative gets smaller
with more
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
then we’re good. If it gets bigger, the error term could potentially
grow overall with
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n").
Since there’s nothing we can do about the derivative (we can’t choose
it!) we should try to do what we can to minimize what we can control,
namely the remaining part of the error term.

Let’s compare two approaches in Mathematica. Without loss of generality,
we’ll live on the interval
![\[-1,1\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5B-1%2C1%5D "[-1,1]").
Please see the demos in [this
notebook](https://drive.google.com/file/d/0B3Www1z6Tm8xNFlXOVhZZDFITlU/view?usp=sharing).

Now let’s summarize the comparison of equally-spaced nodes vs. Chebyshev
nodes.

``` r
nvec <- 1:30
x <- seq(from=-1,to=1,length=5000)
equallyspaced <- NULL
for (n in nvec){
  nodes <- seq(from=-1,to=1,length=n)
  prod <- 1
  for (i in 1:n){
    prod <- prod*(x-nodes[i])
  }
  equallyspaced <- c(equallyspaced,unique(max(abs(prod))))
}
chebyshev <- 1/2^(nvec-1)
plot(nvec,log10(chebyshev),col="green",pch=16,ylim=c(-10,0),xlab="points",ylab="bound on portion of error")
points(nvec,log10(equallyspaced),col="red",pch=16)
```

![](coursenotes_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

Chebyshev is much better! We won’t prove the result, but I hope my
numerical experiment convinced you that for
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
points on
![\[0,1\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5B0%2C1%5D "[0,1]"),
we should choose nodes

![x_i = \cos \frac{(2i-1)\pi}{2n}, \quad i = 1,\ldots,n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_i%20%3D%20%5Ccos%20%5Cfrac%7B%282i-1%29%5Cpi%7D%7B2n%7D%2C%20%5Cquad%20i%20%3D%201%2C%5Cldots%2Cn "x_i = \cos \frac{(2i-1)\pi}{2n}, \quad i = 1,\ldots,n")

and the maximum magnitude of the relevant portion of the error term is

![\max\_{-1\leq x\leq 1}\left\|\prod\_{i=1}^n (x-x_i)\right\| = \frac{1}{2^{n-1}}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmax_%7B-1%5Cleq%20x%5Cleq%201%7D%5Cleft%7C%5Cprod_%7Bi%3D1%7D%5En%20%28x-x_i%29%5Cright%7C%20%3D%20%5Cfrac%7B1%7D%7B2%5E%7Bn-1%7D%7D. "\max_{-1\leq x\leq 1}\left|\prod_{i=1}^n (x-x_i)\right| = \frac{1}{2^{n-1}}.")

If we generalize these results to the interval
![\[a,b\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ba%2Cb%5D "[a,b]"),
we find

![x_i = \frac{b+a}{2} + \frac{b-a}{2}\cos \frac{(2i-1)\pi}{2n}, \quad i = 1,\ldots,n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_i%20%3D%20%5Cfrac%7Bb%2Ba%7D%7B2%7D%20%2B%20%5Cfrac%7Bb-a%7D%7B2%7D%5Ccos%20%5Cfrac%7B%282i-1%29%5Cpi%7D%7B2n%7D%2C%20%5Cquad%20i%20%3D%201%2C%5Cldots%2Cn "x_i = \frac{b+a}{2} + \frac{b-a}{2}\cos \frac{(2i-1)\pi}{2n}, \quad i = 1,\ldots,n")

with

![\max\_{a\leq x\leq b}\left\|\prod\_{i=1}^n (x-x_i)\right\| = \frac{\left(\frac{b-a}{2}\right)^n}{2^{n-1}}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmax_%7Ba%5Cleq%20x%5Cleq%20b%7D%5Cleft%7C%5Cprod_%7Bi%3D1%7D%5En%20%28x-x_i%29%5Cright%7C%20%3D%20%5Cfrac%7B%5Cleft%28%5Cfrac%7Bb-a%7D%7B2%7D%5Cright%29%5En%7D%7B2%5E%7Bn-1%7D%7D. "\max_{a\leq x\leq b}\left|\prod_{i=1}^n (x-x_i)\right| = \frac{\left(\frac{b-a}{2}\right)^n}{2^{n-1}}.")

The proof is not direct, which is why I have eliminated it here.

## Comparing interpolation methods

Let’s do an example comparing interpolation approaches. Consider
![f(x)=(1/\sqrt{2\pi})\mathrm{e}^{-x^2/2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28x%29%3D%281%2F%5Csqrt%7B2%5Cpi%7D%29%5Cmathrm%7Be%7D%5E%7B-x%5E2%2F2%7D "f(x)=(1/\sqrt{2\pi})\mathrm{e}^{-x^2/2}").
This is the standard normal distribution which plays a key role in
probability and statistics.

One thing to know is that the derivatives of this function grow with
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n").

``` r
n <- 0:10
maxderiv <- c(0.398942,0.241971,0.178032, 0.550588,1.19683,2.30711,4.24061,14.178,41.8889,115.091,302.425) # Computed in Mathematica
plot(n,log(maxderiv))
```

![](coursenotes_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

Given this, using equally-spaced nodes seems reckless, but we can try it
anyway, say, on
![\[-10,10\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5B-10%2C10%5D "[-10,10]")
with 30 data points to start with.

``` r
a <- -10
b <- 10
xexact <- seq(from=a,to=b,length=10000)
f <- function(x){exp(-x^2/2)/sqrt(2*pi)}
yexact <- f(xexact)
plot(xexact,yexact,type="l",lwd=3,xlim=c(a,b),ylim=c(-0.2,0.5))
n <- 30
xequal <- seq(from=a,to=b,length=n)
yequal <- lagrange(xequal,f(xequal),xexact)
lines(xexact,yequal,col="red",lwd=2)
```

![](coursenotes_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

Ok, that approach is not going to work! Let’s try the approach of a
lookup table with linear interpolation. Let’s suppose we wish to achieve
4 digit accuracy.

``` r
n <- 1
error <- Inf
while (error > 0.5e-4){
  n <- n+1
  xlookup <- seq(from=a,to=b,length=n)
  ylookup <- f(xlookup)
  ytable <- approx(xlookup,ylookup,xexact)$y
  error <- max(abs(yexact-ytable))
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
  odds <- seq(from=1,to=2*n-1,by=2)
  xcheb <- sort((b+a)/2 + (b-a)/2*cos(odds*pi/2/n))
  ycheb <- lagrange(xcheb, f(xcheb), xexact)
  error <- max(abs(yexact-ycheb))
}
ncheb <- n
print(ncheb)
```

    ## [1] 41

This is an improvement in compression by a factor of 632/41 =
15.4146341.

# Splines

## Big picture

I’ve tried to convince you that it can be problematic to construct
interpolating polynomials of high degree. When dealing with a lot of
data, an alternative approach can be to construct low degree
interpolating polynomials through successive sets of points. Typically
we use cubics, and these are called cubic splines.

## Goals

-   Explain advantages of interpolating data with cubic splines
-   Implement cubic spline interpolation

## Why cubic splines?

Before progressing to real data later on, let’s do an illustrative
example with a small amount of fake data. We make some data points and
connect them with linear splines. We can do this using the built-in
`approxfun` command which returns a function representing the linear
spline.

``` r
x <- c(-2,-1.5,-1,0.25,1,2,3.75,4,5)
y <- c(4,4.2,3,5,0,-2,2,1,1)
xplot <- seq(from=-2,to=5,length=200)
linearspline <- approxfun(x,y)
plot(x,y,ylim=c(-2.5,5.5))
lines(xplot,linearspline(xplot),col="red",lwd=2)
```

![](coursenotes_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

Your eyeball might be telling you that this is a very jagged graph. Most
things in nature and society are not this jagged, so it might feel
desirable to represent the data with smoother functions. Let me show you
what this looks like.

``` r
cubicspline <- splinefun(x,y)
plot(x,y,ylim=c(-2.5,5.5))
lines(xplot,linearspline(xplot),col="red",lwd=2)
lines(xplot,cubicspline(xplot),col="green",lwd=2)
```

![](coursenotes_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

This is smoother. What smoothness means here is continuity of
derivatives from one spline to the next. Using cubic rather than linear
splines gives us more coefficients, and using these coefficients we can
make the derivatives of successive splines match up. How exactly does
this work though?

## Mathematical conditions for cubic splines

Suppose we have data points
![(x_1,y_1),\ldots,(x_n,y_n)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28x_1%2Cy_1%29%2C%5Cldots%2C%28x_n%2Cy_n%29 "(x_1,y_1),\ldots,(x_n,y_n)").
We will connect successive points with cubic curves. With
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
points, there are
![n-1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n-1 "n-1")
curves. Each curve is cubic, having the form

![S_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S_i%28x%29%20%3D%20a_i%20%2B%20b_i%28x-x_i%29%20%2B%20c_i%28x-x_i%29%5E2%20%2B%20d_i%28x-x_i%29%5E3. "S_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3.")

With
![n-1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n-1 "n-1")
curves each having four coefficients, there are
![4(n-1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;4%28n-1%29 "4(n-1)")
coefficients.

First, we force each spline to pass through its two endpoints (which
also makes the overall spline curve continuous). By inspection, this
forces
![a_i=y_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_i%3Dy_i "a_i=y_i"),
and it also enforces a condition on the relationship between
![b_i, c_i, d_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b_i%2C%20c_i%2C%20d_i "b_i, c_i, d_i").
So doing the bookkeeping, that is
![2(n-1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2%28n-1%29 "2(n-1)")
conditions to enforce, which leaves
![2n-2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2n-2 "2n-2")
coefficients undetermined.

Next, we force continuity of the first derivatives of successive
splines. Note that at
![x_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_1 "x_1")
there’s no condition to enforce, since there’s no spline to the left of
it. And at
![x_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_n "x_n")
there’s no condition to enforce, since there’s no spline to the right of
it. Therefore, we enforce conditions at the
![n-2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n-2 "n-2")
points
![x_2,\ldots,x\_{n-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_2%2C%5Cldots%2Cx_%7Bn-1%7D "x_2,\ldots,x_{n-1}").
Specifically, the condition is
![S_i^'(x_i)=S\_{i-1}^'(x_i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S_i%5E%27%28x_i%29%3DS_%7Bi-1%7D%5E%27%28x_i%29 "S_i^'(x_i)=S_{i-1}^'(x_i)").
Subtracting these
![n-2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n-2 "n-2")
conditions from our previous count of
![2n-2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2n-2 "2n-2"),
there are
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
conditions left.

Finally, we enforce continuity of the second derivatives of successive
splines. This is very similar to enforcing the previous condition. The
condition is
![S_i^{''}(x_i)=S\_{i-1}^{''}(x_i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S_i%5E%7B%27%27%7D%28x_i%29%3DS_%7Bi-1%7D%5E%7B%27%27%7D%28x_i%29 "S_i^{''}(x_i)=S_{i-1}^{''}(x_i)")
at the points
![x_2,\ldots,x\_{n-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_2%2C%5Cldots%2Cx_%7Bn-1%7D "x_2,\ldots,x_{n-1}").
Substracting these
![n-2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n-2 "n-2")
conditions from our previous total of
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
conditions, there are two conditions left.

In short, thus far, the procedure for finding spline coefficients is
underetermined because there are
![4n-4](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;4n-4 "4n-4")
unknowns but only
![4n-2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;4n-2 "4n-2")
equations. To make the system solvable, we have to make a choice about
what conditions to enforce at the leftmost and rightmost points. Some of
the most common choices are:

-   **Natural spline**. The concavity at the left and right endpoints is
    zero, that is
    ![S_1^{''}(x_1)=S\_{n-1}^{''}(x_n)=0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S_1%5E%7B%27%27%7D%28x_1%29%3DS_%7Bn-1%7D%5E%7B%27%27%7D%28x_n%29%3D0 "S_1^{''}(x_1)=S_{n-1}^{''}(x_n)=0").
-   **Clamped sline**. The concavity at the left and right endpoints is
    set to a user specified value, that is,
    ![S_1^{''}(x_1)=m_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S_1%5E%7B%27%27%7D%28x_1%29%3Dm_1 "S_1^{''}(x_1)=m_1")
    and
    ![S\_{n-1}^{''}(x_n)=m_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S_%7Bn-1%7D%5E%7B%27%27%7D%28x_n%29%3Dm_n "S_{n-1}^{''}(x_n)=m_n").
-   **FMM (not-a-knot)**.
    ![S_1=S_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S_1%3DS_2 "S_1=S_2")
    is a single cubic equation that is run through the first 3 points,
    and
    ![S\_{n-2}=S\_{n-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S_%7Bn-2%7D%3DS_%7Bn-1%7D "S_{n-2}=S_{n-1}")
    is, similarly, a single cubic run through the last 3 points.

Crucially, since we find the coefficients by solving an square linear
system of dimension
![4n-4](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;4n-4 "4n-4"),
these boundary conditions don’t merely affect the first and last
splines, but influence ALL of the splines.

You don’t need to memorize the details of these different types of
splines. My main goals for you are to understand what they mean and to
be able to implement them in R.

## Splines and linear algebra

As we have been discussing, to find spline coefficients, we have to
solve a linear system. I won’t write down the whole system here because
we are going to use built-in tools to solve it. However, it’s good
mathematical literacy to know that by writing down the system of
equations, you can see that it is tridiagonal and strictly diagonally
dominant, which are nice numerical properties.

## Implementing cubic spline interpolation

Just to emphasize how splines avoid the problem of high-degree
polynomial interpolation, let’s do a cooked example.

``` r
set.seed(123)
n <- 30
x <- sort(runif(n))
y <- cumsum(abs(rnorm(n)))
plot(x,y,pch=19,ylim=c(-2,35),cex=1.5)
xx = seq(from=min(x),to=max(x),length=1000)
yy = barylag(x,y,xx)
lines(xx,yy,col="red",lwd=3)
cubicspline <- splinefun(x,y,method='natural')
lines(xx,cubicspline(xx),col="blue",lwd=3)
```

![](coursenotes_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

Now let’s work with some real data, let’s say, Tesla stock price for the
last 100 days (we won’t actually get 100 days because of days when the
markets were closed, including weekends).

``` r
# Set dates and stock symbol
first.date <- Sys.Date()-100
last.date <- Sys.Date()
tickers <- c('TSLA')
# Acquire data
l.out <- BatchGetSymbols(tickers = tickers, first.date = first.date, last.date = last.date)
```

    ## Warning: `BatchGetSymbols()` was deprecated in BatchGetSymbols 2.6.4.
    ## Please use `yfR::yf_get()` instead.
    ## 2022-05-01: Package BatchGetSymbols will soon be replaced by yfR. 
    ## More details about the change is available at github <<www.github.com/msperlin/yfR>
    ## You can install yfR by executing:
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
sampledday <- day[seq(from=1,to=length(price),by=5)]
sampledprice <- price[seq(from=1,to=length(price),by=5)]
# Fit interpolating polynomial
interpolatedprice <- barylag(sampledday,sampledprice,day[1:max(sampledday)])
plot(day[1:max(sampledday)],interpolatedprice,col="red",type="l",ylim=c(800,1200))
points(day,price,col="blue")
```

![](coursenotes_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

``` r
# Fit splines
TSLAspline <- splinefun(sampledday,sampledprice,method='natural')
plot(day,TSLAspline(day),col="blue",type="l",ylim=c(800,1200))
points(day,price,col="blue")
```

![](coursenotes_files/figure-gfm/unnamed-chunk-42-2.png)<!-- -->

# Least squares

## Big picture

All of our study of interpolation has been based on the idea that the
model (a polynomial or polynomnial spline) should precisely pass through
every data point. But what if we loosen this restriction? What if the
data is thought to have error or noise that we don’t want to represent
in our model? Or what if we just don’t care about precise interpolation,
and decide approximation even at the data points themselves is good
enough?

## Goals

-   Explain how least squares arises as a model fitting problem
-   Explain how the solution of least squares involves vector projection
-   Set up and solve least squares problems
-   Define, identify, and calculate crucial quantities such as basis,
    target, normal equations, least squares solution, residual,
    pseudoinverse, projection, projection operator
-   Recognize when a least squares approach is appropriate
-   Explain how least squares allows data compression

## Model fitting

By way of motivation, let’s examine a pedagogical data set. Suppose
![a](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a "a")
represents the amount of money (in $1,000’s) a company spent on
advertising during different quarters, and
![s](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s "s")
represents money the company earned on sales that quarter. We can plot
the data to explore it.

``` r
a <- c(3,4,5,6)
s <- c(105,117,141,152)
plot(a,s,xlab="advertising",ylab="sales")
```

![](coursenotes_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

The company would like to model this data so they can predict sales for
other levels of advertising. The data looks roughly linear, and we have
no reason to expect a complicated relationship, so let’s try modeling
the data with a line,
![s = x_0 + x_1 a](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s%20%3D%20x_0%20%2B%20x_1%20a "s = x_0 + x_1 a")
where
![x\_{0,1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_%7B0%2C1%7D "x_{0,1}")
are unknown coefficients. Plugging in to the model, we find

![\begin{align}
x_0 + 3x_1 &= 105\\\\
x_0 + 4x_1 &= 117\\\\
x_0 + 5x_1 &= 141\\\\
x_0 + 6x_1 &= 152
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0Ax_0%20%2B%203x_1%20%26%3D%20105%5C%5C%0Ax_0%20%2B%204x_1%20%26%3D%20117%5C%5C%0Ax_0%20%2B%205x_1%20%26%3D%20141%5C%5C%0Ax_0%20%2B%206x_1%20%26%3D%20152%0A%5Cend%7Balign%7D "\begin{align}
x_0 + 3x_1 &= 105\\
x_0 + 4x_1 &= 117\\
x_0 + 5x_1 &= 141\\
x_0 + 6x_1 &= 152
\end{align}")

which we can write in vector form as

![x_0 \begin{pmatrix} 1 \\\\ 1 \\\\ 1 \\\\ 1 \end{pmatrix} + x_1 \begin{pmatrix} 3 \\\\ 4 \\\\ 5 \\\\ 6 \end{pmatrix} = \begin{pmatrix} 105 \\\\ 117 \\\\ 141 \\\\ 152 \end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_0%20%5Cbegin%7Bpmatrix%7D%201%20%5C%5C%201%20%5C%5C%201%20%5C%5C%201%20%5Cend%7Bpmatrix%7D%20%2B%20x_1%20%5Cbegin%7Bpmatrix%7D%203%20%5C%5C%204%20%5C%5C%205%20%5C%5C%206%20%5Cend%7Bpmatrix%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20105%20%5C%5C%20117%20%5C%5C%20141%20%5C%5C%20152%20%5Cend%7Bpmatrix%7D. "x_0 \begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix} + x_1 \begin{pmatrix} 3 \\ 4 \\ 5 \\ 6 \end{pmatrix} = \begin{pmatrix} 105 \\ 117 \\ 141 \\ 152 \end{pmatrix}.")

By writing it this way, we can remember one interpretation of linear
systems. In this case we should imagine a four-dimensional space in
which we are trying to reach a particular target (the right hand side
vector) by taking linear combinations of the two vecors on the left.
Unfortunately, those two vectors only span a two dimensional subspace,
so the chance that we can make it to our target vector are pretty slim.
Stated differently: the sysem is overdetermined. But we would still like
to find a good model, so what should we do?

## Projection onto a vector

To examine the details, let’s start with an even more fundamental
example: a single vector in the plane. Suppose I hand you the vector
![\textbf{A} = (2,1)^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%3D%20%282%2C1%29%5ET "\textbf{A} = (2,1)^T")
and tell you to use it to reach the target vector
![\textbf{b}=(6,8)^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bb%7D%3D%286%2C8%29%5ET "\textbf{b}=(6,8)^T").
Well, you can’t do it exactly because there is no scalar
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
such that
![x\textbf{A} = \textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%5Ctextbf%7BA%7D%20%3D%20%5Ctextbf%7Bb%7D "x\textbf{A} = \textbf{b}").
So let us do the next best thing: let’s find the value of
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
such that
![x\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%5Ctextbf%7BA%7D "x\textbf{A}")
is as close as possible to
![\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bb%7D "\textbf{b}").
We can draw a picture to solve this problem.

![](coursenotes_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

Where should we stop on the dotted line? When we are perpendicular to
the end of
![\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bb%7D "\textbf{b}").
This results in the following picture.

![](coursenotes_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->

From this picture, two relationships arise:

![\begin{align}
x\textbf{A} + \mathbf{r}&=\textbf{b}\\\\
\textbf{A} \cdot \mathbf{r} &=0.
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0Ax%5Ctextbf%7BA%7D%20%2B%20%5Cmathbf%7Br%7D%26%3D%5Ctextbf%7Bb%7D%5C%5C%0A%5Ctextbf%7BA%7D%20%5Ccdot%20%5Cmathbf%7Br%7D%20%26%3D0.%0A%5Cend%7Balign%7D "\begin{align}
x\textbf{A} + \mathbf{r}&=\textbf{b}\\
\textbf{A} \cdot \mathbf{r} &=0.
\end{align}")

Let us dot the first equation with the vector
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}").
Then we have

![\textbf{A} \cdot \textbf{A}x + \textbf{A} \cdot \mathbf{r} = \textbf{A} \cdot \textbf{b}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%5Ccdot%20%5Ctextbf%7BA%7Dx%20%2B%20%5Ctextbf%7BA%7D%20%5Ccdot%20%5Cmathbf%7Br%7D%20%3D%20%5Ctextbf%7BA%7D%20%5Ccdot%20%5Ctextbf%7Bb%7D. "\textbf{A} \cdot \textbf{A}x + \textbf{A} \cdot \mathbf{r} = \textbf{A} \cdot \textbf{b}.")

From the second equation above, we can eliminate the second term on the
left hand side to write

![\textbf{A}^T \textbf{A} x = \textbf{A}^T \textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7BA%7D%20x%20%3D%20%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7Bb%7D "\textbf{A}^T \textbf{A} x = \textbf{A}^T \textbf{b}")

where we have used the fact that
![\mathbf{y}\cdot\mathbf{z}=\mathbf{y}^T\mathbf{z}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7By%7D%5Ccdot%5Cmathbf%7Bz%7D%3D%5Cmathbf%7By%7D%5ET%5Cmathbf%7Bz%7D "\mathbf{y}\cdot\mathbf{z}=\mathbf{y}^T\mathbf{z}").
We can solve for
![\textbf{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D "\textbf{x}")
by writing
![x = \left(\textbf{A}^T \textbf{A}\right)^{-1}\textbf{A}^T \textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%20%3D%20%5Cleft%28%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7BA%7D%5Cright%29%5E%7B-1%7D%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7Bb%7D "x = \left(\textbf{A}^T \textbf{A}\right)^{-1}\textbf{A}^T \textbf{b}").
We can also calculate the vector that was as close as possible to
![\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bb%7D "\textbf{b}").
We will call it
![\widehat{\textbf{b}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidehat%7B%5Ctextbf%7Bb%7D%7D "\widehat{\textbf{b}}")
and it is

![\begin{align}
\widehat{\textbf{b}} &= \textbf{A}x \\\\
&=\textbf{A}\left(\textbf{A}^T \textbf{A}\right)^{-1}\textbf{A}^T \textbf{b}\\\\
& = \mathbf{P} \textbf{b}.
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0A%5Cwidehat%7B%5Ctextbf%7Bb%7D%7D%20%26%3D%20%5Ctextbf%7BA%7Dx%20%5C%5C%0A%26%3D%5Ctextbf%7BA%7D%5Cleft%28%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7BA%7D%5Cright%29%5E%7B-1%7D%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7Bb%7D%5C%5C%0A%26%20%3D%20%5Cmathbf%7BP%7D%20%5Ctextbf%7Bb%7D.%0A%5Cend%7Balign%7D "\begin{align}
\widehat{\textbf{b}} &= \textbf{A}x \\
&=\textbf{A}\left(\textbf{A}^T \textbf{A}\right)^{-1}\textbf{A}^T \textbf{b}\\
& = \mathbf{P} \textbf{b}.
\end{align}")

where thre last equation defines a new quantity that we call
![\mathbf{P}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BP%7D "\mathbf{P}").

Let us know revisit what we have done and emphasize/introduce some
vocabulary. We started with a vector
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
that we used as a **basis** to try to reach the **target**
![\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bb%7D "\textbf{b}").
We couldn’t do it exactly, so we calcualted the closest we could come to
![\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bb%7D "\textbf{b}"),
which turned out to be
![\widehat{\textbf{b}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidehat%7B%5Ctextbf%7Bb%7D%7D "\widehat{\textbf{b}}").
This is called the **projection** of
![\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bb%7D "\textbf{b}")
into the subspace spanned by
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}").
We found
![\widehat{\textbf{b}} = \textbf{A}x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidehat%7B%5Ctextbf%7Bb%7D%7D%20%3D%20%5Ctextbf%7BA%7Dx "\widehat{\textbf{b}} = \textbf{A}x"),
where
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
is called the **least squares solution**, which solved the **normal
equations**
![\textbf{A}^T \textbf{A} x = \textbf{A}^T \textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7BA%7D%20x%20%3D%20%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7Bb%7D "\textbf{A}^T \textbf{A} x = \textbf{A}^T \textbf{b}").
We can summarize the calculation of
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
by remembering
![x = \left(\textbf{A}^T \textbf{A}\right)^{-1}\textbf{A}^T \textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%20%3D%20%5Cleft%28%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7BA%7D%5Cright%29%5E%7B-1%7D%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7Bb%7D "x = \left(\textbf{A}^T \textbf{A}\right)^{-1}\textbf{A}^T \textbf{b}")
where
![\left(\textbf{A}^T \textbf{A}\right)^{-1}\textbf{A}^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7BA%7D%5Cright%29%5E%7B-1%7D%5Ctextbf%7BA%7D%5ET "\left(\textbf{A}^T \textbf{A}\right)^{-1}\textbf{A}^T")
is called the **pseudoinverse** of
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}").
Also, we can summarize the calculation of
![\widehat{\textbf{b}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidehat%7B%5Ctextbf%7Bb%7D%7D "\widehat{\textbf{b}}")
as
![\widehat{\textbf{b}} = \mathbf{P} \textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidehat%7B%5Ctextbf%7Bb%7D%7D%20%3D%20%5Cmathbf%7BP%7D%20%5Ctextbf%7Bb%7D "\widehat{\textbf{b}} = \mathbf{P} \textbf{b}")
where
![\mathbf{P} = \textbf{A} (\textbf{A}^T\textbf{A})^{-1} \textbf{A}^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BP%7D%20%3D%20%5Ctextbf%7BA%7D%20%28%5Ctextbf%7BA%7D%5ET%5Ctextbf%7BA%7D%29%5E%7B-1%7D%20%5Ctextbf%7BA%7D%5ET "\mathbf{P} = \textbf{A} (\textbf{A}^T\textbf{A})^{-1} \textbf{A}^T")
is what we call a **projection operator** or a **projection matrix**.
Since we didn’t succeed in reaching
![\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bb%7D "\textbf{b}"),
there is some error, and we call this the **residual**,
![\mathbf{r} = \textbf{b}-\widehat{\textbf{b}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Br%7D%20%3D%20%5Ctextbf%7Bb%7D-%5Cwidehat%7B%5Ctextbf%7Bb%7D%7D "\mathbf{r} = \textbf{b}-\widehat{\textbf{b}}").

What are the words/ideas you should make sure you understand in the
narrative above?

-   basis
-   target
-   normal equations
-   least squares solution
-   pseudoinverse
-   projection
-   projection operator
-   residual

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
dot <- function(v1,v2){sum(v1*v2)}
dot(a,r)
```

    ## [1] 0

## Projection onto a plane

Let’s apply these same conceps to our original problem of predicting
sales from advertising. We had

![x_0 \begin{pmatrix} 1 \\\\ 1 \\\\ 1 \\\\ 1 \end{pmatrix} + x_1 \begin{pmatrix} 3 \\\\ 4 \\\\ 5 \\\\ 6 \end{pmatrix} = \begin{pmatrix} 105 \\\\ 117 \\\\ 141 \\\\ 152 \end{pmatrix}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_0%20%5Cbegin%7Bpmatrix%7D%201%20%5C%5C%201%20%5C%5C%201%20%5C%5C%201%20%5Cend%7Bpmatrix%7D%20%2B%20x_1%20%5Cbegin%7Bpmatrix%7D%203%20%5C%5C%204%20%5C%5C%205%20%5C%5C%206%20%5Cend%7Bpmatrix%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20105%20%5C%5C%20117%20%5C%5C%20141%20%5C%5C%20152%20%5Cend%7Bpmatrix%7D "x_0 \begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix} + x_1 \begin{pmatrix} 3 \\ 4 \\ 5 \\ 6 \end{pmatrix} = \begin{pmatrix} 105 \\ 117 \\ 141 \\ 152 \end{pmatrix}")

which we’ll write symbolically as

![x_0 \textbf{A}\_0 + x_1 \textbf{A}\_1 = \textbf{b}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_0%20%5Ctextbf%7BA%7D_0%20%2B%20x_1%20%5Ctextbf%7BA%7D_1%20%3D%20%5Ctextbf%7Bb%7D. "x_0 \textbf{A}_0 + x_1 \textbf{A}_1 = \textbf{b}.")

We’d like to reach
![\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bb%7D "\textbf{b}")
using the basis vectors
![\textbf{A}\_{0,1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D_%7B0%2C1%7D "\textbf{A}_{0,1}"),
but we can’t, so let’s consider getting as close as possible. The
picture looks something like this.

![](PlaneProjection.png)

This picture gives rise to the equation

![\begin{align}
x_0 \textbf{A}\_0 + x_1 \mathbf{a_1} + \mathbf{r} &= \textbf{b} \\\\
a_0 \cdot \mathbf{r} &= 0 \\\\
a_1 \cdot \mathbf{r} &= 0.
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0Ax_0%20%5Ctextbf%7BA%7D_0%20%2B%20x_1%20%5Cmathbf%7Ba_1%7D%20%2B%20%5Cmathbf%7Br%7D%20%26%3D%20%5Ctextbf%7Bb%7D%20%5C%5C%0Aa_0%20%5Ccdot%20%5Cmathbf%7Br%7D%20%26%3D%200%20%5C%5C%0Aa_1%20%5Ccdot%20%5Cmathbf%7Br%7D%20%26%3D%200.%0A%5Cend%7Balign%7D "\begin{align}
x_0 \textbf{A}_0 + x_1 \mathbf{a_1} + \mathbf{r} &= \textbf{b} \\
a_0 \cdot \mathbf{r} &= 0 \\
a_1 \cdot \mathbf{r} &= 0.
\end{align}")

Let’s take the dot product of the first equaton with each basis vector.
We find

![\begin{align}
x_0 \textbf{A}\_0\cdot\textbf{A}\_0 + x_1 \textbf{A}\_0\cdot\textbf{A}\_1 &= \textbf{A}\_0 \cdot \textbf{b}\\\\
x_0 \textbf{A}\_1\cdot\textbf{A}\_0 + x_1 \textbf{A}\_1\cdot\textbf{A}\_1 &= \textbf{A}\_1 \cdot \textbf{b}
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0Ax_0%20%5Ctextbf%7BA%7D_0%5Ccdot%5Ctextbf%7BA%7D_0%20%2B%20x_1%20%5Ctextbf%7BA%7D_0%5Ccdot%5Ctextbf%7BA%7D_1%20%26%3D%20%5Ctextbf%7BA%7D_0%20%5Ccdot%20%5Ctextbf%7Bb%7D%5C%5C%0Ax_0%20%5Ctextbf%7BA%7D_1%5Ccdot%5Ctextbf%7BA%7D_0%20%2B%20x_1%20%5Ctextbf%7BA%7D_1%5Ccdot%5Ctextbf%7BA%7D_1%20%26%3D%20%5Ctextbf%7BA%7D_1%20%5Ccdot%20%5Ctextbf%7Bb%7D%0A%5Cend%7Balign%7D "\begin{align}
x_0 \textbf{A}_0\cdot\textbf{A}_0 + x_1 \textbf{A}_0\cdot\textbf{A}_1 &= \textbf{A}_0 \cdot \textbf{b}\\
x_0 \textbf{A}_1\cdot\textbf{A}_0 + x_1 \textbf{A}_1\cdot\textbf{A}_1 &= \textbf{A}_1 \cdot \textbf{b}
\end{align}")

where we have elminated terms that turn out to be zero thanks to the
second and third equations previously. Note that there is a matrix way
to write this. we can write

![\begin{pmatrix} \textbf{A}\_0^T \\\\ \textbf{A}\_1^T \end{pmatrix} \begin{pmatrix} \textbf{A}\_0 & \textbf{A}\_1 \end{pmatrix} \begin{pmatrix} x_0 \\\\ x_1 \end{pmatrix} = \begin{pmatrix} \textbf{A}\_0^T \\\\ \textbf{A}\_1^T \end{pmatrix} \textbf{b}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Bpmatrix%7D%20%5Ctextbf%7BA%7D_0%5ET%20%5C%5C%20%5Ctextbf%7BA%7D_1%5ET%20%5Cend%7Bpmatrix%7D%20%5Cbegin%7Bpmatrix%7D%20%5Ctextbf%7BA%7D_0%20%26%20%5Ctextbf%7BA%7D_1%20%5Cend%7Bpmatrix%7D%20%5Cbegin%7Bpmatrix%7D%20x_0%20%5C%5C%20x_1%20%5Cend%7Bpmatrix%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20%5Ctextbf%7BA%7D_0%5ET%20%5C%5C%20%5Ctextbf%7BA%7D_1%5ET%20%5Cend%7Bpmatrix%7D%20%5Ctextbf%7Bb%7D. "\begin{pmatrix} \textbf{A}_0^T \\ \textbf{A}_1^T \end{pmatrix} \begin{pmatrix} \textbf{A}_0 & \textbf{A}_1 \end{pmatrix} \begin{pmatrix} x_0 \\ x_1 \end{pmatrix} = \begin{pmatrix} \textbf{A}_0^T \\ \textbf{A}_1^T \end{pmatrix} \textbf{b}.")

If we let
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
represent the matrix with columns
![\textbf{A}\_{0,1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D_%7B0%2C1%7D "\textbf{A}_{0,1}")
and if we let
![\textbf{x}=(x_0,x_1)^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D%3D%28x_0%2Cx_1%29%5ET "\textbf{x}=(x_0,x_1)^T")
then we can write the last equation as

![\textbf{A}^T \textbf{A} \textbf{x} = \textbf{A}^T \textbf{b}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7BA%7D%20%5Ctextbf%7Bx%7D%20%3D%20%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7Bb%7D. "\textbf{A}^T \textbf{A} \textbf{x} = \textbf{A}^T \textbf{b}.")

These are the normal equations. For us, concretely, it looks like

![\begin{align}
\begin{pmatrix} 1 & 1 & 1 & 1\\\\ 3 & 4 & 5 & 6\end{pmatrix} \begin{pmatrix} 1 & 3\\\1 & 4\\\\ 1 & 5\\\\ 1& 6\end{pmatrix}\begin{pmatrix}x_1 \\\\ x_2 \end{pmatrix} &= \begin{pmatrix} 1 & 1 & 1 & 1\\\\ 3 & 4 & 5 & 6\end{pmatrix} \begin{pmatrix} 105 \\\\ 117 \\\\ 141 \\\\ 152\end{pmatrix}\\\\
\begin{pmatrix} 4 & 18 \\\\ 18 & 86 \end{pmatrix} \begin{pmatrix} x_1 \\\\ x_2 \end{pmatrix} &= \begin{pmatrix} 515 \\\\ 2400 \end{pmatrix}
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0A%5Cbegin%7Bpmatrix%7D%201%20%26%201%20%26%201%20%26%201%5C%5C%203%20%26%204%20%26%205%20%26%206%5Cend%7Bpmatrix%7D%20%5Cbegin%7Bpmatrix%7D%201%20%26%203%5C%5C1%20%26%204%5C%5C%201%20%26%205%5C%5C%201%26%206%5Cend%7Bpmatrix%7D%5Cbegin%7Bpmatrix%7Dx_1%20%5C%5C%20x_2%20%5Cend%7Bpmatrix%7D%20%26%3D%20%5Cbegin%7Bpmatrix%7D%201%20%26%201%20%26%201%20%26%201%5C%5C%203%20%26%204%20%26%205%20%26%206%5Cend%7Bpmatrix%7D%20%5Cbegin%7Bpmatrix%7D%20105%20%5C%5C%20117%20%5C%5C%20141%20%5C%5C%20152%5Cend%7Bpmatrix%7D%5C%5C%0A%5Cbegin%7Bpmatrix%7D%204%20%26%2018%20%5C%5C%2018%20%26%2086%20%5Cend%7Bpmatrix%7D%20%5Cbegin%7Bpmatrix%7D%20x_1%20%5C%5C%20x_2%20%5Cend%7Bpmatrix%7D%20%26%3D%20%5Cbegin%7Bpmatrix%7D%20515%20%5C%5C%202400%20%5Cend%7Bpmatrix%7D%0A%5Cend%7Balign%7D "\begin{align}
\begin{pmatrix} 1 & 1 & 1 & 1\\ 3 & 4 & 5 & 6\end{pmatrix} \begin{pmatrix} 1 & 3\\1 & 4\\ 1 & 5\\ 1& 6\end{pmatrix}\begin{pmatrix}x_1 \\ x_2 \end{pmatrix} &= \begin{pmatrix} 1 & 1 & 1 & 1\\ 3 & 4 & 5 & 6\end{pmatrix} \begin{pmatrix} 105 \\ 117 \\ 141 \\ 152\end{pmatrix}\\
\begin{pmatrix} 4 & 18 \\ 18 & 86 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} &= \begin{pmatrix} 515 \\ 2400 \end{pmatrix}
\end{align}")

which has solution

![\begin{pmatrix} x_1 \\\\ x_2  \end{pmatrix} = \begin{pmatrix}  54.5 \\\\ 16.5 \end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Bpmatrix%7D%20x_1%20%5C%5C%20x_2%20%20%5Cend%7Bpmatrix%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20%2054.5%20%5C%5C%2016.5%20%5Cend%7Bpmatrix%7D. "\begin{pmatrix} x_1 \\ x_2  \end{pmatrix} = \begin{pmatrix}  54.5 \\ 16.5 \end{pmatrix}.")

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
dot(r,a0)
```

    ## [1] 1.84741111298e-13

``` r
dot(r,a1)
```

    ## [1] 8.81072992343e-13

``` r
plot(a1,b,xlab="advertising",ylab="sales")
xx <- seq(from=0,to=6,length=200)
lines(xx,horner(as.numeric(x),xx)$y)
```

![](coursenotes_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->

Symbolically, we calculated the least squares solution
![\textbf{x} = \left(\textbf{A}^T \textbf{A}\right)^{-1}\textbf{A}^T \textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D%20%3D%20%5Cleft%28%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7BA%7D%5Cright%29%5E%7B-1%7D%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7Bb%7D "\textbf{x} = \left(\textbf{A}^T \textbf{A}\right)^{-1}\textbf{A}^T \textbf{b}")
where
![\left(\textbf{A}^T \textbf{A}\right)^{-1}\textbf{A}^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7BA%7D%5Cright%29%5E%7B-1%7D%5Ctextbf%7BA%7D%5ET "\left(\textbf{A}^T \textbf{A}\right)^{-1}\textbf{A}^T")
is the pseudoinverse. The projection is
![\widehat{\textbf{b}} = \textbf{A} \textbf{x} = \textbf{A} \left(\textbf{A}^T \textbf{A}\right)^{-1}\textbf{A}^T \textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidehat%7B%5Ctextbf%7Bb%7D%7D%20%3D%20%5Ctextbf%7BA%7D%20%5Ctextbf%7Bx%7D%20%3D%20%5Ctextbf%7BA%7D%20%5Cleft%28%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7BA%7D%5Cright%29%5E%7B-1%7D%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7Bb%7D "\widehat{\textbf{b}} = \textbf{A} \textbf{x} = \textbf{A} \left(\textbf{A}^T \textbf{A}\right)^{-1}\textbf{A}^T \textbf{b}")
where
![\textbf{A} \left(\textbf{A}^T \textbf{A}\right)^{-1}\textbf{A}^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%5Cleft%28%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7BA%7D%5Cright%29%5E%7B-1%7D%5Ctextbf%7BA%7D%5ET "\textbf{A} \left(\textbf{A}^T \textbf{A}\right)^{-1}\textbf{A}^T")
is the projection matrix.

Though our example here has used merely two basis vectors, the ideas
extend to any number.

## Model fitting

The technique we have developed here of solving a least squares problem
by vector projection works in a model fitting context whenever the
unknowns
![\textbf{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D "\textbf{x}")
appear linearly. It doesn’t require that other quantities appear
linearly.

For instance, suppose that we wanted to fit a parabola to the data
![(0,6)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%280%2C6%29 "(0,6)"),
![(1,5)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%281%2C5%29 "(1,5)"),
![(2,2)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%282%2C2%29 "(2,2)"),
and
![(3,2)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%283%2C2%29 "(3,2)").
We choose the model
![y = c_0 + c_1 x + c_2 x^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y%20%3D%20c_0%20%2B%20c_1%20x%20%2B%20c_2%20x%5E2 "y = c_0 + c_1 x + c_2 x^2")
where our unknown is the vector
![\mathbf{c}=(c_0,c_1,c_2)^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bc%7D%3D%28c_0%2Cc_1%2Cc_2%29%5ET "\mathbf{c}=(c_0,c_1,c_2)^T").
Though
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
appears nonlinearly, the coefficients
![c\_{0,1,2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c_%7B0%2C1%2C2%7D "c_{0,1,2}")
don’t, so we are fine! Let’s solve this problem using linear algebra. We
want to solve

![\begin{align}
c_0 + c_1 \cdot 0 + c_2 \cdot 0^2 &= 6 \\\\
c_0 + c_1 \cdot 1 + c_2 \cdot 1^2 &= 5 \\\\
c_0 + c_1 \cdot 2 + c_2 \cdot 2^2 &= 2 \\\\
c_0 + c_1 \cdot 3 + c_2 \cdot 3^2 &= 2
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0Ac_0%20%2B%20c_1%20%5Ccdot%200%20%2B%20c_2%20%5Ccdot%200%5E2%20%26%3D%206%20%5C%5C%0Ac_0%20%2B%20c_1%20%5Ccdot%201%20%2B%20c_2%20%5Ccdot%201%5E2%20%26%3D%205%20%5C%5C%0Ac_0%20%2B%20c_1%20%5Ccdot%202%20%2B%20c_2%20%5Ccdot%202%5E2%20%26%3D%202%20%5C%5C%0Ac_0%20%2B%20c_1%20%5Ccdot%203%20%2B%20c_2%20%5Ccdot%203%5E2%20%26%3D%202%0A%5Cend%7Balign%7D "\begin{align}
c_0 + c_1 \cdot 0 + c_2 \cdot 0^2 &= 6 \\
c_0 + c_1 \cdot 1 + c_2 \cdot 1^2 &= 5 \\
c_0 + c_1 \cdot 2 + c_2 \cdot 2^2 &= 2 \\
c_0 + c_1 \cdot 3 + c_2 \cdot 3^2 &= 2
\end{align}")

So, we calculuate in `R`:

``` r
x <- c(0,1,2,3)
A <- matrix(cbind(x^0,x^1,x^2),nrow=4)
b <- c(6,5,2,12)
pseudoinv <- solve(t(A) %*% A) %*% t(A)
c <- pseudoinv %*% b
bhat <- A %*% c
r <- b - bhat
print(c)
```

    ##       [,1]
    ## [1,]  6.75
    ## [2,] -6.75
    ## [3,]  2.75

``` r
plot(x,b)
xx <- seq(from=0,to=3,length=200)
lines(xx,horner(c,xx)$y)
```

![](coursenotes_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->

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
fit data to
![y = C \mathrm{e}^{kx}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y%20%3D%20C%20%5Cmathrm%7Be%7D%5E%7Bkx%7D "y = C \mathrm{e}^{kx}"),
you could take the log of both sides to obtain
![\ln y = \ln C + kx](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cln%20y%20%3D%20%5Cln%20C%20%2B%20kx "\ln y = \ln C + kx").
By considering the data
![(x,\ln y)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28x%2C%5Cln%20y%29 "(x,\ln y)")
you could take a least squares approach to find
![\ln C](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cln%20C "\ln C")
and
![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k").

Overall, here’s the process one would follow for model fitting in cases
where the model parameters appear in a lienar fashion.

1.  Look at data.
2.  Propose a model.
3.  Force the model, resulting in
    ![\textbf{A}\textbf{x}=\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%3D%5Ctextbf%7Bb%7D "\textbf{A}\textbf{x}=\textbf{b}").
4.  Solve the normal equations
    ![\textbf{A}^T\textbf{A}\textbf{x}=\textbf{A}^T \textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5ET%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%3D%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7Bb%7D "\textbf{A}^T\textbf{A}\textbf{x}=\textbf{A}^T \textbf{b}").
5.  Assess the fit of the model visually and/or using the residual
    vector.

A statistics class would provide much more sophisticated ways of
analyzing the error.

## Least squares and data compression

Adopting a least squares approach allows, potentially, massive
compression of data. Suppose we had 100 points that looked like this

``` r
x <- seq(from=0,to=1,length=100)
y <- 3*x + 0.4*(2*runif(100)-1)
plot(x,y)
```

![](coursenotes_files/figure-gfm/unnamed-chunk-49-1.png)<!-- -->

If we decided to represent this data with a line, we’d go down from
having 100 pieces of information (the original data points) to merely 2
(a slope and an intercept).

## An optimization viewpoint on least squares

So far, we’ve taken a geometric approach to solving least-squares
problems, but there’s another way to get the same result: optimization.
We can start directly with the square of the norm of the residual vector
(a factor of 1/2 is included for algebraic convenience, but it doesn’t
change the result)

![\frac{1}{2}\|\|\mathbf{r}\|\|^2 = \frac{1}{2}\|\|\textbf{A}\textbf{x}-\textbf{b}\|\|^2.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B1%7D%7B2%7D%7C%7C%5Cmathbf%7Br%7D%7C%7C%5E2%20%3D%20%5Cfrac%7B1%7D%7B2%7D%7C%7C%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D-%5Ctextbf%7Bb%7D%7C%7C%5E2. "\frac{1}{2}||\mathbf{r}||^2 = \frac{1}{2}||\textbf{A}\textbf{x}-\textbf{b}||^2.")

More formally, we can define the scalar **objective function**

![f(\textbf{x}) = \frac{1}{2} \|\|\textbf{A}\textbf{x}-\textbf{b}\|\|^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28%5Ctextbf%7Bx%7D%29%20%3D%20%5Cfrac%7B1%7D%7B2%7D%20%7C%7C%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D-%5Ctextbf%7Bb%7D%7C%7C%5E2 "f(\textbf{x}) = \frac{1}{2} ||\textbf{A}\textbf{x}-\textbf{b}||^2")

and define our least squares solution
![\textbf{x}\_{LS}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D_%7BLS%7D "\textbf{x}_{LS}")as
the value of
![\textbf{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D "\textbf{x}")
that minimizes this objective function, that is,

![\textbf{x}\_{LS} = \mathop{\mathrm{arg\\,min}}\_{\textbf{x}} f(\textbf{x}).](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D_%7BLS%7D%20%3D%20%5Cmathop%7B%5Cmathrm%7Barg%5C%2Cmin%7D%7D_%7B%5Ctextbf%7Bx%7D%7D%20f%28%5Ctextbf%7Bx%7D%29. "\textbf{x}_{LS} = \mathop{\mathrm{arg\,min}}_{\textbf{x}} f(\textbf{x}).")

Now you know where the terminology least squares comes from! We can
apply calculus to minimize this expression and the normal equations will
result. This is a problem on your in-class exercise and your homework,
and it is critical that you work through the details.

# QR Factorization

## Big picture

We’ve talked about how to fit a model to data using least squares, and
we’ve examined the theoretical aspect of this process, but thus far
we’ve ignored computational issues. It turns out that solving the normal
equations can be difficult on a computer, and a technique called QR
factorization provides a potentially better way to solve a least squares
problem.

## Goals

-   Define orthonormal matrices
-   Perform Graham-Schmidt Orthogonalization
-   Solve least squares problems with **QR** decomposition

## Orthogonality is nice

The eventual goal of this lesson is to show you how to solve the least
squares problem

![\textbf{x}\_{LS} = \mathop{\mathrm{arg\\,min}}\_x \frac{1}{2}\|\|\textbf{A}\textbf{x}-\textbf{b}\|\|^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D_%7BLS%7D%20%3D%20%5Cmathop%7B%5Cmathrm%7Barg%5C%2Cmin%7D%7D_x%20%5Cfrac%7B1%7D%7B2%7D%7C%7C%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D-%5Ctextbf%7Bb%7D%7C%7C%5E2 "\textbf{x}_{LS} = \mathop{\mathrm{arg\,min}}_x \frac{1}{2}||\textbf{A}\textbf{x}-\textbf{b}||^2")

by writing the matrix
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
in a convenient way. But to build up to that, we need to introduce a
number of ideas. First off: orthogonality.

You might remember that two vectors
![\mathbf{v}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D "\mathbf{v}")
and
![\mathbf{w}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bw%7D "\mathbf{w}")
are **orthogonal** if
![\mathbf{v} \cdot \mathbf{w} = \mathbf{v}^T \mathbf{w} = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D%20%5Ccdot%20%5Cmathbf%7Bw%7D%20%3D%20%5Cmathbf%7Bv%7D%5ET%20%5Cmathbf%7Bw%7D%20%3D%200 "\mathbf{v} \cdot \mathbf{w} = \mathbf{v}^T \mathbf{w} = 0").
We say that two vectors are **orthonormal** if they are orthogonal and
each have norm of one, that is
![\mathbf{v} \cdot \mathbf{v} = \mathbf{w} \cdot \mathbf{w} = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D%20%5Ccdot%20%5Cmathbf%7Bv%7D%20%3D%20%5Cmathbf%7Bw%7D%20%5Ccdot%20%5Cmathbf%7Bw%7D%20%3D%201 "\mathbf{v} \cdot \mathbf{v} = \mathbf{w} \cdot \mathbf{w} = 1").

A matrix is called **orthogonal** or **orthonormal** (they are used
interchangeably for matrices, sometimes) if its columns are orthonormal
vectors. A cool property arises from this. Suppose these orthonormal
columns are
![\mathbf{v}\_1,\ldots,\mathbf{v}\_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_1%2C%5Cldots%2C%5Cmathbf%7Bv%7D_n "\mathbf{v}_1,\ldots,\mathbf{v}_n").
Then we can consider the quantity
![\textbf{A}^T \textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7BA%7D "\textbf{A}^T \textbf{A}"):

![\begin{align}
\textbf{A}^T \textbf{A} &= \begin{pmatrix} \mathbf{v}\_1^T\\\\ \mathbf{v}\_2^T \\\\ \vdots \\\\ \mathbf{v}\_n^T \end{pmatrix} \begin{pmatrix} \mathbf{v}\_1 & \mathbf{v}\_2 & \cdots & \mathbf{v}\_n \end{pmatrix} \\\\
&= \begin{pmatrix} \mathbf{v}\_1^T \mathbf{v}\_1 & \mathbf{v}\_1^T \mathbf{v}\_2 & \cdots & \mathbf{v}\_1^T \mathbf{v}\_n \\\\ \mathbf{v}\_2^T \mathbf{v}\_1 & \mathbf{v}\_2^T \mathbf{v}\_2 & \cdots &  \mathbf{v}\_2^T \mathbf{v}\_n \\\\
\vdots & \vdots & \cdots & \vdots \\\\
\mathbf{v}\_n^T \mathbf{v}\_1 & \mathbf{v}\_n^T \mathbf{v}\_2 & \cdots & \mathbf{v}\_n^T \mathbf{v}\_n
\end{pmatrix}\\\\
&= \begin{pmatrix} 1 & 0 & \cdots & 0 \\\\ 0 & 1 & \cdots & 0\\\\ \vdots & \vdots & \cdots & \vdots \\\\ 0 & 0 & \cdots & 1\end{pmatrix}\\\\
&= \mathbf{I}\_n.
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0A%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7BA%7D%20%26%3D%20%5Cbegin%7Bpmatrix%7D%20%5Cmathbf%7Bv%7D_1%5ET%5C%5C%20%5Cmathbf%7Bv%7D_2%5ET%20%5C%5C%20%5Cvdots%20%5C%5C%20%5Cmathbf%7Bv%7D_n%5ET%20%5Cend%7Bpmatrix%7D%20%5Cbegin%7Bpmatrix%7D%20%5Cmathbf%7Bv%7D_1%20%26%20%5Cmathbf%7Bv%7D_2%20%26%20%5Ccdots%20%26%20%5Cmathbf%7Bv%7D_n%20%5Cend%7Bpmatrix%7D%20%5C%5C%0A%26%3D%20%5Cbegin%7Bpmatrix%7D%20%5Cmathbf%7Bv%7D_1%5ET%20%5Cmathbf%7Bv%7D_1%20%26%20%5Cmathbf%7Bv%7D_1%5ET%20%5Cmathbf%7Bv%7D_2%20%26%20%5Ccdots%20%26%20%5Cmathbf%7Bv%7D_1%5ET%20%5Cmathbf%7Bv%7D_n%20%5C%5C%20%5Cmathbf%7Bv%7D_2%5ET%20%5Cmathbf%7Bv%7D_1%20%26%20%5Cmathbf%7Bv%7D_2%5ET%20%5Cmathbf%7Bv%7D_2%20%26%20%5Ccdots%20%26%20%20%5Cmathbf%7Bv%7D_2%5ET%20%5Cmathbf%7Bv%7D_n%20%5C%5C%0A%5Cvdots%20%26%20%5Cvdots%20%26%20%5Ccdots%20%26%20%5Cvdots%20%5C%5C%0A%5Cmathbf%7Bv%7D_n%5ET%20%5Cmathbf%7Bv%7D_1%20%26%20%5Cmathbf%7Bv%7D_n%5ET%20%5Cmathbf%7Bv%7D_2%20%26%20%5Ccdots%20%26%20%5Cmathbf%7Bv%7D_n%5ET%20%5Cmathbf%7Bv%7D_n%0A%5Cend%7Bpmatrix%7D%5C%5C%0A%26%3D%20%5Cbegin%7Bpmatrix%7D%201%20%26%200%20%26%20%5Ccdots%20%26%200%20%5C%5C%200%20%26%201%20%26%20%5Ccdots%20%26%200%5C%5C%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Ccdots%20%26%20%5Cvdots%20%5C%5C%200%20%26%200%20%26%20%5Ccdots%20%26%201%5Cend%7Bpmatrix%7D%5C%5C%0A%26%3D%20%5Cmathbf%7BI%7D_n.%0A%5Cend%7Balign%7D "\begin{align}
\textbf{A}^T \textbf{A} &= \begin{pmatrix} \mathbf{v}_1^T\\ \mathbf{v}_2^T \\ \vdots \\ \mathbf{v}_n^T \end{pmatrix} \begin{pmatrix} \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n \end{pmatrix} \\
&= \begin{pmatrix} \mathbf{v}_1^T \mathbf{v}_1 & \mathbf{v}_1^T \mathbf{v}_2 & \cdots & \mathbf{v}_1^T \mathbf{v}_n \\ \mathbf{v}_2^T \mathbf{v}_1 & \mathbf{v}_2^T \mathbf{v}_2 & \cdots &  \mathbf{v}_2^T \mathbf{v}_n \\
\vdots & \vdots & \cdots & \vdots \\
\mathbf{v}_n^T \mathbf{v}_1 & \mathbf{v}_n^T \mathbf{v}_2 & \cdots & \mathbf{v}_n^T \mathbf{v}_n
\end{pmatrix}\\
&= \begin{pmatrix} 1 & 0 & \cdots & 0 \\ 0 & 1 & \cdots & 0\\ \vdots & \vdots & \cdots & \vdots \\ 0 & 0 & \cdots & 1\end{pmatrix}\\
&= \mathbf{I}_n.
\end{align}")

If
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
is square and orthonormal, since
![\textbf{A}^T\textbf{A} = \mathbf{I}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5ET%5Ctextbf%7BA%7D%20%3D%20%5Cmathbf%7BI%7D "\textbf{A}^T\textbf{A} = \mathbf{I}"),
then by definition of inverse,
![\textbf{A}^T = \textbf{A}^{-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5ET%20%3D%20%5Ctextbf%7BA%7D%5E%7B-1%7D "\textbf{A}^T = \textbf{A}^{-1}").
This is a pretty great way to calculate an inverse! Remember, though,
that in our least squares context
,![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
is generally *not* a square matrix.

Another helpful property of orthogonal matrices is that they preserve
length, meaning that have norm of one. That is, suppose
![\mathbf{Q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BQ%7D "\mathbf{Q}")
is orthogonal. Then

![\begin{align}
\|\|\mathbf{Q} \mathbf{v}\|\|^2 &= \mathbf{Q} \mathbf{v} \cdot \mathbf{Q} \mathbf{v} \\\\
& = (\mathbf{Q}\mathbf{v})^T (\mathbf{Q}\mathbf{v})\\\\
& = \mathbf{v}^T \mathbf{Q}^T \mathbf{Q} \mathbf{v}\\\\
& = \mathbf{v}^T \mathbf{Q}^{-1} \mathbf{Q} \mathbf{v}\\\\
& = \mathbf{v}^T \mathbf{v}\\\\
& = \|\|\mathbf{v}\|\|^2.
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0A%7C%7C%5Cmathbf%7BQ%7D%20%5Cmathbf%7Bv%7D%7C%7C%5E2%20%26%3D%20%5Cmathbf%7BQ%7D%20%5Cmathbf%7Bv%7D%20%5Ccdot%20%5Cmathbf%7BQ%7D%20%5Cmathbf%7Bv%7D%20%5C%5C%0A%26%20%3D%20%28%5Cmathbf%7BQ%7D%5Cmathbf%7Bv%7D%29%5ET%20%28%5Cmathbf%7BQ%7D%5Cmathbf%7Bv%7D%29%5C%5C%0A%26%20%3D%20%5Cmathbf%7Bv%7D%5ET%20%5Cmathbf%7BQ%7D%5ET%20%5Cmathbf%7BQ%7D%20%5Cmathbf%7Bv%7D%5C%5C%0A%26%20%3D%20%5Cmathbf%7Bv%7D%5ET%20%5Cmathbf%7BQ%7D%5E%7B-1%7D%20%5Cmathbf%7BQ%7D%20%5Cmathbf%7Bv%7D%5C%5C%0A%26%20%3D%20%5Cmathbf%7Bv%7D%5ET%20%5Cmathbf%7Bv%7D%5C%5C%0A%26%20%3D%20%7C%7C%5Cmathbf%7Bv%7D%7C%7C%5E2.%0A%5Cend%7Balign%7D "\begin{align}
||\mathbf{Q} \mathbf{v}||^2 &= \mathbf{Q} \mathbf{v} \cdot \mathbf{Q} \mathbf{v} \\
& = (\mathbf{Q}\mathbf{v})^T (\mathbf{Q}\mathbf{v})\\
& = \mathbf{v}^T \mathbf{Q}^T \mathbf{Q} \mathbf{v}\\
& = \mathbf{v}^T \mathbf{Q}^{-1} \mathbf{Q} \mathbf{v}\\
& = \mathbf{v}^T \mathbf{v}\\
& = ||\mathbf{v}||^2.
\end{align}")

Since
![\|\|\mathbf{Q}\mathbf{v}\|\|^2 = \|\|\mathbf{v}\|\|^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%7C%5Cmathbf%7BQ%7D%5Cmathbf%7Bv%7D%7C%7C%5E2%20%3D%20%7C%7C%5Cmathbf%7Bv%7D%7C%7C%5E2 "||\mathbf{Q}\mathbf{v}||^2 = ||\mathbf{v}||^2"),
we know
![\|\|\mathbf{Q}\mathbf{v}\|\| = \|\|\mathbf{v}\|\|](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%7C%5Cmathbf%7BQ%7D%5Cmathbf%7Bv%7D%7C%7C%20%3D%20%7C%7C%5Cmathbf%7Bv%7D%7C%7C "||\mathbf{Q}\mathbf{v}|| = ||\mathbf{v}||")
and therefore
![\|\|\mathbf{Q}\|\|=1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%7C%5Cmathbf%7BQ%7D%7C%7C%3D1 "||\mathbf{Q}||=1").
A similar calculuation shows that
![\|\|\mathbf{Q}^{-1}\|\|=1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%7C%5Cmathbf%7BQ%7D%5E%7B-1%7D%7C%7C%3D1 "||\mathbf{Q}^{-1}||=1").
And then, by definition of condition number,
![\kappa(\mathbf{Q}) = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ckappa%28%5Cmathbf%7BQ%7D%29%20%3D%201 "\kappa(\mathbf{Q}) = 1").
This means that orthogonal matrices are incredibly well-conditioned.

Finally, let’s consider projecting a vector into a subspace spanned by
orthonormal vectors
![\mathbf{q}\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bq%7D_i "\mathbf{q}_i"),
![i=1,\ldots,n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%3D1%2C%5Cldots%2Cn "i=1,\ldots,n").
Define
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
as

![\textbf{A} = \begin{pmatrix} \| & \| &  & \| \\\\ \mathbf{q}\_1 & \mathbf{q}\_2 & \cdots & \\  \mathbf{q}\_n \\\\\| & \| &  & \| \\\\ \end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20%7C%20%26%20%7C%20%26%20%20%26%20%7C%20%5C%5C%20%5Cmathbf%7Bq%7D_1%20%26%20%5Cmathbf%7Bq%7D_2%20%26%20%5Ccdots%20%26%20%5C%20%20%5Cmathbf%7Bq%7D_n%20%5C%5C%7C%20%26%20%7C%20%26%20%20%26%20%7C%20%5C%5C%20%5Cend%7Bpmatrix%7D. "\textbf{A} = \begin{pmatrix} | & | &  & | \\ \mathbf{q}_1 & \mathbf{q}_2 & \cdots & \  \mathbf{q}_n \\| & | &  & | \\ \end{pmatrix}.")

Let’s project a vector
![\mathbf{w}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bw%7D "\mathbf{w}")
onto the subspace spanned by the columns of
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}").
By definition of the projection operator, and using the fact that
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
is orthogonal,

![\begin{align}
\mathbf{P}\mathbf{w} &= \textbf{A} (\textbf{A}^T \textbf{A})^{-1} \textbf{A}^T \mathbf{w}\\\\
&= \textbf{A} \textbf{A}^{-1} \textbf{A} \textbf{A}^T \mathbf{w} \\\\
&= \textbf{A} \textbf{A}^T \mathbf{w} \\\\
&= \mathbf{q}\_1 \mathbf{q}\_1^T \mathbf{w} + \mathbf{q}\_2 \mathbf{q}\_2^T \mathbf{w} + \cdots + \mathbf{q}\_n \mathbf{q}\_n^T \mathbf{w}.
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0A%5Cmathbf%7BP%7D%5Cmathbf%7Bw%7D%20%26%3D%20%5Ctextbf%7BA%7D%20%28%5Ctextbf%7BA%7D%5ET%20%5Ctextbf%7BA%7D%29%5E%7B-1%7D%20%5Ctextbf%7BA%7D%5ET%20%5Cmathbf%7Bw%7D%5C%5C%0A%26%3D%20%5Ctextbf%7BA%7D%20%5Ctextbf%7BA%7D%5E%7B-1%7D%20%5Ctextbf%7BA%7D%20%5Ctextbf%7BA%7D%5ET%20%5Cmathbf%7Bw%7D%20%5C%5C%0A%26%3D%20%5Ctextbf%7BA%7D%20%5Ctextbf%7BA%7D%5ET%20%5Cmathbf%7Bw%7D%20%5C%5C%0A%26%3D%20%5Cmathbf%7Bq%7D_1%20%5Cmathbf%7Bq%7D_1%5ET%20%5Cmathbf%7Bw%7D%20%2B%20%5Cmathbf%7Bq%7D_2%20%5Cmathbf%7Bq%7D_2%5ET%20%5Cmathbf%7Bw%7D%20%2B%20%5Ccdots%20%2B%20%5Cmathbf%7Bq%7D_n%20%5Cmathbf%7Bq%7D_n%5ET%20%5Cmathbf%7Bw%7D.%0A%5Cend%7Balign%7D "\begin{align}
\mathbf{P}\mathbf{w} &= \textbf{A} (\textbf{A}^T \textbf{A})^{-1} \textbf{A}^T \mathbf{w}\\
&= \textbf{A} \textbf{A}^{-1} \textbf{A} \textbf{A}^T \mathbf{w} \\
&= \textbf{A} \textbf{A}^T \mathbf{w} \\
&= \mathbf{q}_1 \mathbf{q}_1^T \mathbf{w} + \mathbf{q}_2 \mathbf{q}_2^T \mathbf{w} + \cdots + \mathbf{q}_n \mathbf{q}_n^T \mathbf{w}.
\end{align}")

So, when projecting onto the space spanned by orthonormal vectors, you
can just project onto each vector separately and add up the results. In
the calculation above, remember that
![\mathbf{q}\_i \mathbf{q}\_i^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bq%7D_i%20%5Cmathbf%7Bq%7D_i%5ET "\mathbf{q}_i \mathbf{q}_i^T")
is a matrix.

As an example, take

![\mathbf{q}\_1 = \left(\frac{1}{3},\frac{2}{3},\frac{2}{3}\right)^T, \quad \mathbf{q}\_2 = \left(\frac{2}{15},\frac{2}{3},-\frac{11}{15}\right)^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bq%7D_1%20%3D%20%5Cleft%28%5Cfrac%7B1%7D%7B3%7D%2C%5Cfrac%7B2%7D%7B3%7D%2C%5Cfrac%7B2%7D%7B3%7D%5Cright%29%5ET%2C%20%5Cquad%20%5Cmathbf%7Bq%7D_2%20%3D%20%5Cleft%28%5Cfrac%7B2%7D%7B15%7D%2C%5Cfrac%7B2%7D%7B3%7D%2C-%5Cfrac%7B11%7D%7B15%7D%5Cright%29%5ET "\mathbf{q}_1 = \left(\frac{1}{3},\frac{2}{3},\frac{2}{3}\right)^T, \quad \mathbf{q}_2 = \left(\frac{2}{15},\frac{2}{3},-\frac{11}{15}\right)^T")

so that

![\textbf{A} = \begin{pmatrix} \frac{1}{3} & \frac{2}{15} \\\\ \frac{2}{3} & \frac{2}{3} \\\\ \frac{2}{3} & -\frac{11}{15} \end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20%5Cfrac%7B1%7D%7B3%7D%20%26%20%5Cfrac%7B2%7D%7B15%7D%20%5C%5C%20%5Cfrac%7B2%7D%7B3%7D%20%26%20%5Cfrac%7B2%7D%7B3%7D%20%5C%5C%20%5Cfrac%7B2%7D%7B3%7D%20%26%20-%5Cfrac%7B11%7D%7B15%7D%20%5Cend%7Bpmatrix%7D. "\textbf{A} = \begin{pmatrix} \frac{1}{3} & \frac{2}{15} \\ \frac{2}{3} & \frac{2}{3} \\ \frac{2}{3} & -\frac{11}{15} \end{pmatrix}.")

Then

![\begin{align}
\textbf{A} \textbf{A}^T &= \begin{pmatrix} \frac{1}{3} & \frac{2}{15} \\\\ \frac{2}{3} & \frac{2}{3} \\\\ \frac{2}{3} & -\frac{11}{15} \end{pmatrix}
\begin{pmatrix} \frac{1}{3} &  \frac{2}{3} & \frac{2}{3} \\\\ \frac{2}{15} & \frac{2}{3} & -\frac{11}{15} \end{pmatrix} \\\\
&= \frac{1}{225}
\begin{pmatrix}
29 & 70 & 28 \\\\
70 & 200 & -10 \\\\
28 & -10 & 221
\end{pmatrix}\\\\
& = \underbrace{
\begin{pmatrix}
\frac{1}{9} & \frac{2}{9} & \frac{2}{9} \\\\
\frac{2}{9} & \frac{4}{9} & \frac{4}{9} \\\\
\frac{2}{9} & \frac{4}{9} & \frac{4}{9}
\end{pmatrix}}\_{\mathbf{q}\_1 \mathbf{q}\_1^T}  +  \underbrace{
\begin{pmatrix}
\frac{4}{225} & \frac{4}{45} & -\frac{22}{225} \\\\
\frac{4}{45} & \frac{4}{9} & -\frac{22}{45} \\\\
-\frac{22}{225} & -\frac{22}{45} & \frac{121}{225}
\end{pmatrix}}\_{\mathbf{q}\_2 \mathbf{q}\_2^T}.
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0A%5Ctextbf%7BA%7D%20%5Ctextbf%7BA%7D%5ET%20%26%3D%20%5Cbegin%7Bpmatrix%7D%20%5Cfrac%7B1%7D%7B3%7D%20%26%20%5Cfrac%7B2%7D%7B15%7D%20%5C%5C%20%5Cfrac%7B2%7D%7B3%7D%20%26%20%5Cfrac%7B2%7D%7B3%7D%20%5C%5C%20%5Cfrac%7B2%7D%7B3%7D%20%26%20-%5Cfrac%7B11%7D%7B15%7D%20%5Cend%7Bpmatrix%7D%0A%5Cbegin%7Bpmatrix%7D%20%5Cfrac%7B1%7D%7B3%7D%20%26%20%20%5Cfrac%7B2%7D%7B3%7D%20%26%20%5Cfrac%7B2%7D%7B3%7D%20%5C%5C%20%5Cfrac%7B2%7D%7B15%7D%20%26%20%5Cfrac%7B2%7D%7B3%7D%20%26%20-%5Cfrac%7B11%7D%7B15%7D%20%5Cend%7Bpmatrix%7D%20%5C%5C%0A%26%3D%20%5Cfrac%7B1%7D%7B225%7D%0A%5Cbegin%7Bpmatrix%7D%0A29%20%26%2070%20%26%2028%20%5C%5C%0A70%20%26%20200%20%26%20-10%20%5C%5C%0A28%20%26%20-10%20%26%20221%0A%5Cend%7Bpmatrix%7D%5C%5C%0A%26%20%3D%20%5Cunderbrace%7B%0A%5Cbegin%7Bpmatrix%7D%0A%5Cfrac%7B1%7D%7B9%7D%20%26%20%5Cfrac%7B2%7D%7B9%7D%20%26%20%5Cfrac%7B2%7D%7B9%7D%20%5C%5C%0A%5Cfrac%7B2%7D%7B9%7D%20%26%20%5Cfrac%7B4%7D%7B9%7D%20%26%20%5Cfrac%7B4%7D%7B9%7D%20%5C%5C%0A%5Cfrac%7B2%7D%7B9%7D%20%26%20%5Cfrac%7B4%7D%7B9%7D%20%26%20%5Cfrac%7B4%7D%7B9%7D%0A%5Cend%7Bpmatrix%7D%7D_%7B%5Cmathbf%7Bq%7D_1%20%5Cmathbf%7Bq%7D_1%5ET%7D%20%20%2B%20%20%5Cunderbrace%7B%0A%5Cbegin%7Bpmatrix%7D%0A%5Cfrac%7B4%7D%7B225%7D%20%26%20%5Cfrac%7B4%7D%7B45%7D%20%26%20-%5Cfrac%7B22%7D%7B225%7D%20%5C%5C%0A%5Cfrac%7B4%7D%7B45%7D%20%26%20%5Cfrac%7B4%7D%7B9%7D%20%26%20-%5Cfrac%7B22%7D%7B45%7D%20%5C%5C%0A-%5Cfrac%7B22%7D%7B225%7D%20%26%20-%5Cfrac%7B22%7D%7B45%7D%20%26%20%5Cfrac%7B121%7D%7B225%7D%0A%5Cend%7Bpmatrix%7D%7D_%7B%5Cmathbf%7Bq%7D_2%20%5Cmathbf%7Bq%7D_2%5ET%7D.%0A%5Cend%7Balign%7D "\begin{align}
\textbf{A} \textbf{A}^T &= \begin{pmatrix} \frac{1}{3} & \frac{2}{15} \\ \frac{2}{3} & \frac{2}{3} \\ \frac{2}{3} & -\frac{11}{15} \end{pmatrix}
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
\end{align}")

## Gram-Schmidt orthogonalization

Now that I’ve convinced you that it’s nice to have an orthonormal basis,
let’s talk about how to get one. Gram-Schmidt orthogonalization is an
iterative way of creating an orthonormal basis from a set of vectors.

Let’s see how it works via an example. Suppose we have
![\mathbf{v}\_1 = (3,4,0)^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_1%20%3D%20%283%2C4%2C0%29%5ET "\mathbf{v}_1 = (3,4,0)^T")
and
![\mathbf{v}\_2 = (1,1,0)^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_2%20%3D%20%281%2C1%2C0%29%5ET "\mathbf{v}_2 = (1,1,0)^T").

``` r
v1 <- c(1,2,2)
v2 <- c(2,1,-2)
v3 <- c(1,1,0)
```

1.  Take the first vector and turn it into a unit vector.

``` r
y1 <- v1
r11 <- Norm(y1,2)
q1 <- y1/r11
```

2.  Think of
    ![\mathbf{v}\_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_2 "\mathbf{v}_2")
    as made up of stuff in the subspace spanned by
    ![\mathbf{q}\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bq%7D_1 "\mathbf{q}_1")
    and stuff orthogonal to it. Through away stuff in the span of
    ![\mathbf{q}\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bq%7D_1 "\mathbf{q}_1")
    since we have it covered already. Take what’s left of
    ![\mathbf{v}\_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_2 "\mathbf{v}_2")
    and turn it unto a unit vector.

``` r
y2 <- v2 - q1%*%t(q1)%*%v2
r22 <- Norm(y2,2)
q2 <- y2/r22
```

3.  Think of
    ![\mathbf{v}\_3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_3 "\mathbf{v}_3")
    as made up of stuff in the subspace spanned by
    ![\mathbf{q}\_{1,2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bq%7D_%7B1%2C2%7D "\mathbf{q}_{1,2}")
    and stuff orthogonal to it. Through away stuff in the span of
    ![\mathbf{q}\_{1,2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bq%7D_%7B1%2C2%7D "\mathbf{q}_{1,2}")
    since we have it covered already. Take what’s left of
    ![\mathbf{v}\_3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_3 "\mathbf{v}_3")
    and turn it unto a unit vector.

``` r
y3 <- v3 - q1%*%t(q1)%*%v3 - q2%*%t(q2)%*%v3
r33 <- Norm(y3,2)
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
Norm(q1,2)
```

    ## [1] 1

``` r
Norm(q2,2)
```

    ## [1] 1

``` r
dot(q1,q2)
```

    ## [1] 0

## QR decomposition

Just like LU decomposition is a way of using matrices to encode the
process of Gaussian elimination, QR decomposition is a way of using
matrices to encode the process of Gram-Schmidt orthogonalization. Let’s
take the equations we implemented above:

![\begin{align}
\mathbf{y}\_1 &= \mathbf{v}\_1\\\\
\mathbf{q}\_1 &= \frac{1}{\|\|\mathbf{y}\_1\|\|} \mathbf{y}\_1 \\\\
\mathbf{y}\_2 &= \mathbf{v}\_2 - \mathbf{q}\_1 \mathbf{q}\_1^T \mathbf{v}\_2 \\\\
\mathbf{q}\_2 &= \frac{1}{\|\|\mathbf{y}\_2\|\|} \mathbf{y}\_2\\\\
\mathbf{y}\_3 &= \mathbf{v}\_3 - ( \mathbf{q}\_1 \mathbf{q}\_1^T \mathbf{v}\_3 + \mathbf{q}\_2 \mathbf{q}\_2^T \mathbf{v}\_3)\\\\
\mathbf{q}\_3 &= \frac{1}{\|\|\mathbf{y}\_3\|\|} \mathbf{y}\_3.
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0A%5Cmathbf%7By%7D_1%20%26%3D%20%5Cmathbf%7Bv%7D_1%5C%5C%0A%5Cmathbf%7Bq%7D_1%20%26%3D%20%5Cfrac%7B1%7D%7B%7C%7C%5Cmathbf%7By%7D_1%7C%7C%7D%20%5Cmathbf%7By%7D_1%20%5C%5C%0A%5Cmathbf%7By%7D_2%20%26%3D%20%5Cmathbf%7Bv%7D_2%20-%20%5Cmathbf%7Bq%7D_1%20%5Cmathbf%7Bq%7D_1%5ET%20%5Cmathbf%7Bv%7D_2%20%5C%5C%0A%5Cmathbf%7Bq%7D_2%20%26%3D%20%5Cfrac%7B1%7D%7B%7C%7C%5Cmathbf%7By%7D_2%7C%7C%7D%20%5Cmathbf%7By%7D_2%5C%5C%0A%5Cmathbf%7By%7D_3%20%26%3D%20%5Cmathbf%7Bv%7D_3%20-%20%28%20%5Cmathbf%7Bq%7D_1%20%5Cmathbf%7Bq%7D_1%5ET%20%5Cmathbf%7Bv%7D_3%20%2B%20%5Cmathbf%7Bq%7D_2%20%5Cmathbf%7Bq%7D_2%5ET%20%5Cmathbf%7Bv%7D_3%29%5C%5C%0A%5Cmathbf%7Bq%7D_3%20%26%3D%20%5Cfrac%7B1%7D%7B%7C%7C%5Cmathbf%7By%7D_3%7C%7C%7D%20%5Cmathbf%7By%7D_3.%0A%5Cend%7Balign%7D "\begin{align}
\mathbf{y}_1 &= \mathbf{v}_1\\
\mathbf{q}_1 &= \frac{1}{||\mathbf{y}_1||} \mathbf{y}_1 \\
\mathbf{y}_2 &= \mathbf{v}_2 - \mathbf{q}_1 \mathbf{q}_1^T \mathbf{v}_2 \\
\mathbf{q}_2 &= \frac{1}{||\mathbf{y}_2||} \mathbf{y}_2\\
\mathbf{y}_3 &= \mathbf{v}_3 - ( \mathbf{q}_1 \mathbf{q}_1^T \mathbf{v}_3 + \mathbf{q}_2 \mathbf{q}_2^T \mathbf{v}_3)\\
\mathbf{q}_3 &= \frac{1}{||\mathbf{y}_3||} \mathbf{y}_3.
\end{align}")

Let’s rewrite these equations, solving for the
![\mathbf{v}\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_i "\mathbf{v}_i"):

![\begin{align}
\mathbf{v}\_1 &= \|\|\mathbf{y}\_1\|\| \mathbf{q}\_1 \\\\
\mathbf{v}\_2 &= (\mathbf{q}\_1\cdot \mathbf{v}\_2) \mathbf{q}\_1 + \|\|\mathbf{y}\_2\|\| \mathbf{q}\_2 \\\\
\mathbf{v}\_3 &= (\mathbf{q}\_1\cdot \mathbf{v}\_3) \mathbf{q}\_1 + (\mathbf{q}\_2\cdot \mathbf{v}\_3) \mathbf{q}\_2  + \|\|\mathbf{y}\_3\|\| \mathbf{q}\_3 
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0A%5Cmathbf%7Bv%7D_1%20%26%3D%20%7C%7C%5Cmathbf%7By%7D_1%7C%7C%20%5Cmathbf%7Bq%7D_1%20%5C%5C%0A%5Cmathbf%7Bv%7D_2%20%26%3D%20%28%5Cmathbf%7Bq%7D_1%5Ccdot%20%5Cmathbf%7Bv%7D_2%29%20%5Cmathbf%7Bq%7D_1%20%2B%20%7C%7C%5Cmathbf%7By%7D_2%7C%7C%20%5Cmathbf%7Bq%7D_2%20%5C%5C%0A%5Cmathbf%7Bv%7D_3%20%26%3D%20%28%5Cmathbf%7Bq%7D_1%5Ccdot%20%5Cmathbf%7Bv%7D_3%29%20%5Cmathbf%7Bq%7D_1%20%2B%20%28%5Cmathbf%7Bq%7D_2%5Ccdot%20%5Cmathbf%7Bv%7D_3%29%20%5Cmathbf%7Bq%7D_2%20%20%2B%20%7C%7C%5Cmathbf%7By%7D_3%7C%7C%20%5Cmathbf%7Bq%7D_3%20%0A%5Cend%7Balign%7D "\begin{align}
\mathbf{v}_1 &= ||\mathbf{y}_1|| \mathbf{q}_1 \\
\mathbf{v}_2 &= (\mathbf{q}_1\cdot \mathbf{v}_2) \mathbf{q}_1 + ||\mathbf{y}_2|| \mathbf{q}_2 \\
\mathbf{v}_3 &= (\mathbf{q}_1\cdot \mathbf{v}_3) \mathbf{q}_1 + (\mathbf{q}_2\cdot \mathbf{v}_3) \mathbf{q}_2  + ||\mathbf{y}_3|| \mathbf{q}_3 
\end{align}")

or, as matrices,

![\begin{pmatrix} \| & \| &   \| \\\\ \mathbf{v}\_1 & \mathbf{v}\_2  & \mathbf{v}\_3 \\\\\| & \| &   \| \\\\ \end{pmatrix} = 
\underbrace{ \begin{pmatrix} \| & \| &   \| \\\\ \mathbf{q}\_1 & \mathbf{q}\_2  & \mathbf{q}\_3 \\\\\| & \| &   \| \\\\ \end{pmatrix}}\_Q   \underbrace{\begin{pmatrix} r\_{1,1} & r\_{1,2} &   r\_{1,3} \\\\ 0 &r\_{2,2}  & r\_{2,3} \\\0 & 0 &   r\_{3,3} \\\\ \end{pmatrix} }\_R](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Bpmatrix%7D%20%7C%20%26%20%7C%20%26%20%20%20%7C%20%5C%5C%20%5Cmathbf%7Bv%7D_1%20%26%20%5Cmathbf%7Bv%7D_2%20%20%26%20%5Cmathbf%7Bv%7D_3%20%5C%5C%7C%20%26%20%7C%20%26%20%20%20%7C%20%5C%5C%20%5Cend%7Bpmatrix%7D%20%3D%20%0A%5Cunderbrace%7B%20%5Cbegin%7Bpmatrix%7D%20%7C%20%26%20%7C%20%26%20%20%20%7C%20%5C%5C%20%5Cmathbf%7Bq%7D_1%20%26%20%5Cmathbf%7Bq%7D_2%20%20%26%20%5Cmathbf%7Bq%7D_3%20%5C%5C%7C%20%26%20%7C%20%26%20%20%20%7C%20%5C%5C%20%5Cend%7Bpmatrix%7D%7D_Q%20%20%20%5Cunderbrace%7B%5Cbegin%7Bpmatrix%7D%20r_%7B1%2C1%7D%20%26%20r_%7B1%2C2%7D%20%26%20%20%20r_%7B1%2C3%7D%20%5C%5C%200%20%26r_%7B2%2C2%7D%20%20%26%20r_%7B2%2C3%7D%20%5C%5C0%20%26%200%20%26%20%20%20r_%7B3%2C3%7D%20%5C%5C%20%5Cend%7Bpmatrix%7D%20%7D_R "\begin{pmatrix} | & | &   | \\ \mathbf{v}_1 & \mathbf{v}_2  & \mathbf{v}_3 \\| & | &   | \\ \end{pmatrix} = 
\underbrace{ \begin{pmatrix} | & | &   | \\ \mathbf{q}_1 & \mathbf{q}_2  & \mathbf{q}_3 \\| & | &   | \\ \end{pmatrix}}_Q   \underbrace{\begin{pmatrix} r_{1,1} & r_{1,2} &   r_{1,3} \\ 0 &r_{2,2}  & r_{2,3} \\0 & 0 &   r_{3,3} \\ \end{pmatrix} }_R")

with

![r\_{i,i} = \|\|\mathbf{y}\_i\|\|, \qquad r\_{i,j} = \mathbf{q}\_i \cdot \mathbf{v}\_j.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r_%7Bi%2Ci%7D%20%3D%20%7C%7C%5Cmathbf%7By%7D_i%7C%7C%2C%20%5Cqquad%20r_%7Bi%2Cj%7D%20%3D%20%5Cmathbf%7Bq%7D_i%20%5Ccdot%20%5Cmathbf%7Bv%7D_j. "r_{i,i} = ||\mathbf{y}_i||, \qquad r_{i,j} = \mathbf{q}_i \cdot \mathbf{v}_j.")

The QR decomposition is a matrix decomposition that writes a
![m \times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m%20%5Ctimes%20n "m \times n")
matrix
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
as a product
![\textbf{A} = \mathbf{Q} \mathbf{R}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%3D%20%5Cmathbf%7BQ%7D%20%5Cmathbf%7BR%7D "\textbf{A} = \mathbf{Q} \mathbf{R}")
where:

-   ![\mathbf{Q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BQ%7D "\mathbf{Q}")
    is an
    ![m \times r](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m%20%5Ctimes%20r "m \times r")
    matrix with orthonormal columns, where
    ![r](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r "r")
    is the number of linearly independent columns of
    ![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}").
-   ![\mathbf{R}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BR%7D "\mathbf{R}")
    is an
    ![r \times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r%20%5Ctimes%20n "r \times n")
    matrix which is upper triangular if
    ![r=n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r%3Dn "r=n"),
    or the top portion of an upper triangular matrix if
    ![r\<n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r%3Cn "r<n").
-   The columns of
    ![\mathbf{Q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BQ%7D "\mathbf{Q}")
    span the same space as the columns of
    ![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}").
-   The matrix
    ![\mathbf{R}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BR%7D "\mathbf{R}")
    gives the change of basis between the vectors in
    ![\mathbf{Q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BQ%7D "\mathbf{Q}")
    and the vectors in
    ![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}").
-   The decompisition is unique up to some sign changes, so if we
    require
    ![R\_{ii}\geq 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R_%7Bii%7D%5Cgeq%200 "R_{ii}\geq 0"),
    it is unique.
-   If the columns of
    ![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
    are independent, then
    ![R\_{ii}\neq 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R_%7Bii%7D%5Cneq%200 "R_{ii}\neq 0").
-   On the other hand, if column
    ![j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;j "j")
    of
    ![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
    can be written as a linear combination of columns to the left, then
    ![R\_{jj}=0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R_%7Bjj%7D%3D0 "R_{jj}=0").

The QR decomposition we’ve done so far is actually called the partial QR
decomposition We distinguish this from the **full** or **complete** QR
decomposition. In the latter, we include vectors than span parts of the
space not spanned by
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
itself. Of course, these contribute nothing to the matrix
![\mathbf{Q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BQ%7D "\mathbf{Q}"),
so it results in a bunch of 0’s in
![\mathbf{R}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BR%7D "\mathbf{R}").
Below is the key picture to understand. Here,
![\mathbf{Q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BQ%7D "\mathbf{Q}")
and
![\mathbf{R}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BR%7D "\mathbf{R}")
are the (partial) QR decomposition and
![\overline{Q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Coverline%7BQ%7D "\overline{Q}")
and
![\overline{R}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Coverline%7BR%7D "\overline{R}")
are the full version.

![](fullQR.png)

Let’s make this clear with numerous examples!

First example. Let’s start with a
![3 \times 3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;3%20%5Ctimes%203 "3 \times 3")
matrix and perform Gram-Schmidt orthogonalization to obtain the QR
decomposition.

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
r12 <- dot(q1,v2)

# Step 3
y3 <- v3 - (q1%*%t(q1))%*%v3 - (q2%*%t(q2))%*%v3
r33 <- Norm(y3,2)
q3 <- y3/r33
r13 <- dot(q1,v3)
r23 <- dot(q2,v3)

# Assemble into Qbar and Rbar
Qbar <- cbind(q1,q2,q3)
Rbar <- rbind(c(r11,r12,r13),c(0,r22,r23),c(0,0,r33))
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

In this example, notice that
![\overline{\mathbf{Q}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Coverline%7B%5Cmathbf%7BQ%7D%7D "\overline{\mathbf{Q}}")
has three columns. This reflects the fact that the columns of
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
are linearly independent. Hence, it spans all of
![\mathbb{R}^3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BR%7D%5E3 "\mathbb{R}^3")
and the QR and full QR decompositions are the same.

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
r11 <- Norm(y1,2)
q1 <- y1/r11

# Step2
y2 <- v2 - (q1%*%t(q1))%*%v2
r22 <- Norm(y2,2)
q2 <- y2/r22
r12 <- dot(q1,v2)

# Step 3
y3 <- v3 - (q1%*%t(q1))%*%v3 - (q2%*%t(q2))%*%v3
y3
```

    ##                   [,1]
    ## [1,] 1.11022302463e-16
    ## [2,] 5.55111512313e-17
    ## [3,] 0.00000000000e+00

Oh! it turns out there’s nothing left. If we want to compute the full QR
decomposition, we need to find something orthogonal to the span of
![q\_{1,2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q_%7B1%2C2%7D "q_{1,2}").
An easy way to do that is to pick a vector outside of their span and
continue the orthogonalization procedure. You can algorithmically use
your linear algebra skills to do this, but for little cases like ours,
you can eyeball it.

``` r
V3 <- c(1,1,1)
y3 <- V3 - (q1%*%t(q1))%*%V3 - (q2%*%t(q2))%*%V3
q3 <- y3/Norm(y3,2)
r13 <- dot(q1,v3)
r23 <- dot(q2,v3)
r33 <- 0

# Assemble into Qbar and Rbar
Qbar <- cbind(q1,q2,q3)
Rbar <- rbind(c(r11,r12,r13),c(0,r22,r23),c(0,0,r33))
Qbar
```

    ##                  q1                                
    ## [1,] 0.333333333333  0.666666666667  0.666666666667
    ## [2,] 0.666666666667  0.333333333333 -0.666666666667
    ## [3,] 0.666666666667 -0.666666666667  0.333333333333

``` r
Rbar
```

    ##      [,1] [,2] [,3]
    ## [1,]    3    0    1
    ## [2,]    0    3    1
    ## [3,]    0    0    0

To reiterate, we ended up with a row of zeros at the bottom of
![\overline{\mathbf{R}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Coverline%7B%5Cmathbf%7BR%7D%7D "\overline{\mathbf{R}}").
That’s because the columns of
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
are linearly dependent and don’t span
![\mathbb{R}^3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BR%7D%5E3 "\mathbb{R}^3").
At any rate, let’s go ahead and check our result.

``` r
Qbar%*%Rbar - A
```

    ##      v1 v2 v3
    ## [1,]  0  0  0
    ## [2,]  0  0  0
    ## [3,]  0  0  0

In this case, the full and reduced QR are not the same. Let’s see how
this works.

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

    ##      [,1] [,2] [,3]
    ## [1,]    3    0    1
    ## [2,]    0    3    1

``` r
Q%*%R - A
```

    ##      v1 v2 v3
    ## [1,]  0  0  0
    ## [2,]  0  0  0
    ## [3,]  0  0  0

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
r11 <- Norm(y1,2)
q1 <- y1/r11

# Step 2
y2 <- v2 - (q1%*%t(q1))%*%v2
r22 <- Norm(y2,2)
q2 <- y2/r22
r12 <- dot(q1,v2)

# Step 3
y3 <- v3 - (q1%*%t(q1))%*%v3 - (q2%*%t(q2))%*%v3
r33 <- Norm(y3,2)
q3 <- y3/(r33)
r13 <- dot(q1,v3)
r23 <- dot(q2,v3)
```

But wait! We are living in
![\mathbb{R}^4](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BR%7D%5E4 "\mathbb{R}^4")
and we only have three vectors to far,
![\mathbf{q}\_{1,2,3}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bq%7D_%7B1%2C2%2C3%7D "\mathbf{q}_{1,2,3}").
If we want the full decomposition, we have to find a basis for the
orthogonal complement of
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}").

``` r
# Choose a vector not in the span of q1, q2, q3
v4 <- c(1,2,3,4)
y4 <- v4 - (q1%*%t(q1))%*%v4 - (q2%*%t(q2))%*%v4  - (q3%*%t(q3))%*%v4
r44 <- Norm(y4,2) 
q4 <- y4/(r44)
r14 <- dot(q1,v4)
r24 <- dot(q2,v4)
r34 <- dot(q3,v4)

# Assemble into Qbar and Rbar, check answer
Qbar <- cbind(q1,q2,q3,q4)
Rbar <- rbind(c(r11,r12,r13),c(0,r22,r23),c(0,0,r33),c(0,0,0))
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
Qbar%*%Rbar - A
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
Q%*%R - A
```

    ##      v1 v2                 v3
    ## [1,]  0  0 -1.23259516441e-32
    ## [2,]  0  0  0.00000000000e+00
    ## [3,]  0  0  0.00000000000e+00
    ## [4,]  0  0  0.00000000000e+00

## Computational considerations

So far, we have calculated the QR decomposition using a methodology that
is based on Gram-Schmidt orthogonalization. This is because Gram-Schmidt
is the most conceptually straightforward method. For an
![m \times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m%20%5Ctimes%20n "m \times n")
matrix, it turns out that the computational cost of factorization is
![\mathcal{O})(2mn^2)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathcal%7BO%7D%29%282mn%5E2%29 "\mathcal{O})(2mn^2)")
multiplications and additions. An alternative choice for factorization
uses objects called Householder reflections. I won’t go in to these
here, but the method requires fewer operations and is the one actually
implemented in many software packages.

## What is the point of QR factorization

Finally, we ask: why have we bothered to do all of this? It is really,
really convenient for least squares, and turns out to have very nice
numerical properties because of small condition numbers. You’ll work
with the numerical issue on your activities and/or homework, but for
now, here’s how least squares works when you use QR decomposition on
![\textbf{A}\textbf{x}=\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%3D%5Ctextbf%7Bb%7D "\textbf{A}\textbf{x}=\textbf{b}").
Note

![\begin{align}
\|\|\textbf{A}\textbf{x}-\textbf{b}\|\|^2 &= \|\|\bar{\mathbf{Q}}\bar{\mathbf{R}} \textbf{x} - \textbf{b}\|\|^2 \\\\
&= \|\|\bar{\mathbf{Q}}^T(\bar{\mathbf{Q}}\bar{\mathbf{R}} \textbf{x} - \textbf{b})\|\|^2 \\\\
&= \|\|\bar{\mathbf{Q}}^T \bar{\mathbf{Q}}\bar{\mathbf{R}} \textbf{x} - \bar{\mathbf{Q}}^T \textbf{b}\|\|^2 \\\\
&= \|\|\bar{\mathbf{R}} \textbf{x} - \bar{\mathbf{Q}}^T \textbf{b}\|\|^2 \\\\
& = \|\|\mathbf{R} \textbf{x}-\mathbf{Q}^T \textbf{b}\|\|^2+ \|\|\widehat{\mathbf{Q}}^T \textbf{b}\|\|^2. 
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0A%7C%7C%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D-%5Ctextbf%7Bb%7D%7C%7C%5E2%20%26%3D%20%7C%7C%5Cbar%7B%5Cmathbf%7BQ%7D%7D%5Cbar%7B%5Cmathbf%7BR%7D%7D%20%5Ctextbf%7Bx%7D%20-%20%5Ctextbf%7Bb%7D%7C%7C%5E2%20%5C%5C%0A%26%3D%20%7C%7C%5Cbar%7B%5Cmathbf%7BQ%7D%7D%5ET%28%5Cbar%7B%5Cmathbf%7BQ%7D%7D%5Cbar%7B%5Cmathbf%7BR%7D%7D%20%5Ctextbf%7Bx%7D%20-%20%5Ctextbf%7Bb%7D%29%7C%7C%5E2%20%5C%5C%0A%26%3D%20%7C%7C%5Cbar%7B%5Cmathbf%7BQ%7D%7D%5ET%20%5Cbar%7B%5Cmathbf%7BQ%7D%7D%5Cbar%7B%5Cmathbf%7BR%7D%7D%20%5Ctextbf%7Bx%7D%20-%20%5Cbar%7B%5Cmathbf%7BQ%7D%7D%5ET%20%5Ctextbf%7Bb%7D%7C%7C%5E2%20%5C%5C%0A%26%3D%20%7C%7C%5Cbar%7B%5Cmathbf%7BR%7D%7D%20%5Ctextbf%7Bx%7D%20-%20%5Cbar%7B%5Cmathbf%7BQ%7D%7D%5ET%20%5Ctextbf%7Bb%7D%7C%7C%5E2%20%5C%5C%0A%26%20%3D%20%7C%7C%5Cmathbf%7BR%7D%20%5Ctextbf%7Bx%7D-%5Cmathbf%7BQ%7D%5ET%20%5Ctextbf%7Bb%7D%7C%7C%5E2%2B%20%7C%7C%5Cwidehat%7B%5Cmathbf%7BQ%7D%7D%5ET%20%5Ctextbf%7Bb%7D%7C%7C%5E2.%20%0A%5Cend%7Balign%7D "\begin{align}
||\textbf{A}\textbf{x}-\textbf{b}||^2 &= ||\bar{\mathbf{Q}}\bar{\mathbf{R}} \textbf{x} - \textbf{b}||^2 \\
&= ||\bar{\mathbf{Q}}^T(\bar{\mathbf{Q}}\bar{\mathbf{R}} \textbf{x} - \textbf{b})||^2 \\
&= ||\bar{\mathbf{Q}}^T \bar{\mathbf{Q}}\bar{\mathbf{R}} \textbf{x} - \bar{\mathbf{Q}}^T \textbf{b}||^2 \\
&= ||\bar{\mathbf{R}} \textbf{x} - \bar{\mathbf{Q}}^T \textbf{b}||^2 \\
& = ||\mathbf{R} \textbf{x}-\mathbf{Q}^T \textbf{b}||^2+ ||\widehat{\mathbf{Q}}^T \textbf{b}||^2. 
\end{align}")

Proceeding from the second to last line to the very last line is not
obvious, and takes you a little bit of writing out expressions more
explicitly. This is a good exercise for you to try.

In any case, to minimize this, note that the second term doesn’t even
include
![\textbf{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bx%7D "\textbf{x}")
so there’s nothing we can do about it (in fact, it is the norm of the
residual). To minimize the quantity, we can then just force the first
term to be zero. This is fine because it is a square system! There’s a
command that does this all automatically for you called `qr.solve`.

# Eigenvalues

## Big picture

Along with solving the linear system
![\textbf{A}\textbf{x}=\textbf{b}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Ctextbf%7Bx%7D%3D%5Ctextbf%7Bb%7D "\textbf{A}\textbf{x}=\textbf{b}"),
finding the eigenvalues and eigenvectors of a matrix
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
is one of the most important problems in linear algebra. Knowing the
eigenpairs can help simplify a problem and reveal important information
about systems modeled with linear algebra.

## Goals

-   Define eigenvalues and eigenvectors, and calculate them by hand for
    small matrices
-   Diagonalize matrices and state the algebraic and geometric
    multiplicity of eigenvalues
-   Apply the power iteration technique to find an eigenvalue

## Eigenvalue fundamentals

For an
![n \times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%20n "n \times n")
matrix
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}"),
a scalar
![\lambda \in \mathbb{C}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda%20%5Cin%20%5Cmathbb%7BC%7D "\lambda \in \mathbb{C}"),
and vector
![\mathbf{v} \in \mathbb{R}^n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D%20%5Cin%20%5Cmathbb%7BR%7D%5En "\mathbf{v} \in \mathbb{R}^n"),
![\mathbf{v} \neq \mathbf{0}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D%20%5Cneq%20%5Cmathbf%7B0%7D "\mathbf{v} \neq \mathbf{0}"),
then we say
![\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda")
is an **eigenvalue** of
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
and
![\mathbf{v}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D "\mathbf{v}")
is an **eigenvector** of
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
if
![\textbf{A}\mathbf{v}=\lambda \mathbf{v}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Cmathbf%7Bv%7D%3D%5Clambda%20%5Cmathbf%7Bv%7D "\textbf{A}\mathbf{v}=\lambda \mathbf{v}").

Stated in words: an eigenvector and eigenvalue are the magical vector
![\mathbf{v}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D "\mathbf{v}")
and scalar
![\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda")
such that if you hit
![\mathbf{v}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D "\mathbf{v}")
with
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}"),
you get back the same vector
![\mathbf{v}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D "\mathbf{v}")
but multiplied by a constant
![\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda").

How do we calculate them? Let’s take the definition
![\textbf{A} \mathbf{v} = \lambda \mathbf{v}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%5Cmathbf%7Bv%7D%20%3D%20%5Clambda%20%5Cmathbf%7Bv%7D "\textbf{A} \mathbf{v} = \lambda \mathbf{v}")
and rearrange it to write
![(\textbf{A} - \mathbf{I} \lambda) \mathbf{v} = \mathbf{0}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28%5Ctextbf%7BA%7D%20-%20%5Cmathbf%7BI%7D%20%5Clambda%29%20%5Cmathbf%7Bv%7D%20%3D%20%5Cmathbf%7B0%7D "(\textbf{A} - \mathbf{I} \lambda) \mathbf{v} = \mathbf{0}").
There are only two ways this can happen. One choice is
![\mathbf{v}=0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D%3D0 "\mathbf{v}=0"),
but that’s trivial because it works for any
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}").
The other choice, by the Invertible Matrix Theorem, is that the matrix
on the left is singular. Also from the Invertible Matrix Theorem, if it
is singular, then it has determinant zero, that is
![\det (\textbf{A}- \mathbf{I} \lambda) = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cdet%20%28%5Ctextbf%7BA%7D-%20%5Cmathbf%7BI%7D%20%5Clambda%29%20%3D%200 "\det (\textbf{A}- \mathbf{I} \lambda) = 0").
This equation is a polynomial in
![\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda")
and is called the **characteristic polynomial**. When calculating by
hand, we find the characteristic polynomial first and then solve it to
find the eigenvalues. To find eigenvectors, we remember that
![\textbf{A}\mathbf{v} = \lambda \mathbf{v} \rightarrow (\textbf{A}-\mathbf{I}\lambda)\mathbf{v}=\mathbf{0}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Cmathbf%7Bv%7D%20%3D%20%5Clambda%20%5Cmathbf%7Bv%7D%20%5Crightarrow%20%28%5Ctextbf%7BA%7D-%5Cmathbf%7BI%7D%5Clambda%29%5Cmathbf%7Bv%7D%3D%5Cmathbf%7B0%7D "\textbf{A}\mathbf{v} = \lambda \mathbf{v} \rightarrow (\textbf{A}-\mathbf{I}\lambda)\mathbf{v}=\mathbf{0}")
and solve for
![\mathbf{v}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D "\mathbf{v}").

To recap:

1.  Solve the characteristic equation
    ![\det (\textbf{A}- \mathbf{I} \lambda) = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cdet%20%28%5Ctextbf%7BA%7D-%20%5Cmathbf%7BI%7D%20%5Clambda%29%20%3D%200 "\det (\textbf{A}- \mathbf{I} \lambda) = 0")
    to find the
    ![\lambda_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_i "\lambda_i").
2.  Solve
    ![(\textbf{A}-\mathbf{I}\lambda_i)\mathbf{v_i}=\mathbf{0}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28%5Ctextbf%7BA%7D-%5Cmathbf%7BI%7D%5Clambda_i%29%5Cmathbf%7Bv_i%7D%3D%5Cmathbf%7B0%7D "(\textbf{A}-\mathbf{I}\lambda_i)\mathbf{v_i}=\mathbf{0}")
    to find the
    ![\mathbf{v}\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_i "\mathbf{v}_i").

For example, let’s find the eigenvalues and eigenvectors of

![\textbf{A} = \begin{pmatrix} -3 & 2 \\\\ 2 & -3 \end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20-3%20%26%202%20%5C%5C%202%20%26%20-3%20%5Cend%7Bpmatrix%7D. "\textbf{A} = \begin{pmatrix} -3 & 2 \\ 2 & -3 \end{pmatrix}.")

Using the result above, we can write the characteristic polynomial:

![\begin{align}
\det (\textbf{A} - \mathbf{I} \lambda) &= 0\\\\
\det \begin{pmatrix} -3 - \lambda & 2 \\\\ 2 & -3 - \lambda \end{pmatrix} &= 0\\\\
(\lambda+3)^2-4 &= 0\\\\
\lambda^2 + 6\lambda + 5 &= 0\\\\
(\lambda+5)(\lambda+1) &= 0
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0A%5Cdet%20%28%5Ctextbf%7BA%7D%20-%20%5Cmathbf%7BI%7D%20%5Clambda%29%20%26%3D%200%5C%5C%0A%5Cdet%20%5Cbegin%7Bpmatrix%7D%20-3%20-%20%5Clambda%20%26%202%20%5C%5C%202%20%26%20-3%20-%20%5Clambda%20%5Cend%7Bpmatrix%7D%20%26%3D%200%5C%5C%0A%28%5Clambda%2B3%29%5E2-4%20%26%3D%200%5C%5C%0A%5Clambda%5E2%20%2B%206%5Clambda%20%2B%205%20%26%3D%200%5C%5C%0A%28%5Clambda%2B5%29%28%5Clambda%2B1%29%20%26%3D%200%0A%5Cend%7Balign%7D "\begin{align}
\det (\textbf{A} - \mathbf{I} \lambda) &= 0\\
\det \begin{pmatrix} -3 - \lambda & 2 \\ 2 & -3 - \lambda \end{pmatrix} &= 0\\
(\lambda+3)^2-4 &= 0\\
\lambda^2 + 6\lambda + 5 &= 0\\
(\lambda+5)(\lambda+1) &= 0
\end{align}")

and therefore
![\lambda\_{1,2} = -1, -5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_%7B1%2C2%7D%20%3D%20-1%2C%20-5 "\lambda_{1,2} = -1, -5").
To find
![\mathbf{v}\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_1 "\mathbf{v}_1"),
we solve:

![\begin{align}
(\textbf{A}-\mathbf{I}\lambda_1)\mathbf{v}\_1&=\mathbf{0}\\\\
\begin{pmatrix} -2 & 2 \\\\ 2 & -2 \end{pmatrix}\mathbf{v}\_1 &= \mathbf{0} \\\\
\mathbf{v}\_1 &= \begin{pmatrix} 1 \\\\ 1 \end{pmatrix}
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0A%28%5Ctextbf%7BA%7D-%5Cmathbf%7BI%7D%5Clambda_1%29%5Cmathbf%7Bv%7D_1%26%3D%5Cmathbf%7B0%7D%5C%5C%0A%5Cbegin%7Bpmatrix%7D%20-2%20%26%202%20%5C%5C%202%20%26%20-2%20%5Cend%7Bpmatrix%7D%5Cmathbf%7Bv%7D_1%20%26%3D%20%5Cmathbf%7B0%7D%20%5C%5C%0A%5Cmathbf%7Bv%7D_1%20%26%3D%20%5Cbegin%7Bpmatrix%7D%201%20%5C%5C%201%20%5Cend%7Bpmatrix%7D%0A%5Cend%7Balign%7D "\begin{align}
(\textbf{A}-\mathbf{I}\lambda_1)\mathbf{v}_1&=\mathbf{0}\\
\begin{pmatrix} -2 & 2 \\ 2 & -2 \end{pmatrix}\mathbf{v}_1 &= \mathbf{0} \\
\mathbf{v}_1 &= \begin{pmatrix} 1 \\ 1 \end{pmatrix}
\end{align}")

or any scalar multiple of this vector. Similarly, we find

![\mathbf{v}\_2 = \begin{pmatrix} 1 \\\\ -1 \end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_2%20%3D%20%5Cbegin%7Bpmatrix%7D%201%20%5C%5C%20-1%20%5Cend%7Bpmatrix%7D. "\mathbf{v}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}.")

## Algebraic multiplicity, geometric multiplicity, and diagonalization

Many applications of eigenvalues are intimiately tied up with the idea
of **diagonalization** of matrices. Suppose
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
has eigenpairs
![\lambda_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_i "\lambda_i"),
![\mathbf{v}\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_i "\mathbf{v}_i"),
![i = 1,\ldots,n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%20%3D%201%2C%5Cldots%2Cn "i = 1,\ldots,n").
Then we can write down the definition of eigenpair for all pairs
simultaneously:
![\textbf{A} \mathbf{v}\_i = \lambda_i \mathbf{v}\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%5Cmathbf%7Bv%7D_i%20%3D%20%5Clambda_i%20%5Cmathbf%7Bv%7D_i "\textbf{A} \mathbf{v}_i = \lambda_i \mathbf{v}_i")
implies

![\begin{align}
\textbf{A} \underbrace{\begin{pmatrix}
\vert & \vert & \cdots & \vert \\\\
\mathbf{v}\_1 & \mathbf{v}\_2 & \cdots & \mathbf{v}\_n \\\\
\vert & \vert & \cdots & \vert 
\end{pmatrix}}\_\mathbf{S} &= 
{\begin{pmatrix} \vert & \vert & \cdots & \vert \\\\
\lambda_1 \mathbf{v}\_1 & \lambda_2 \mathbf{v}\_2 & \cdots & \lambda_n \mathbf{v}\_n \\\\
\vert & \vert & \cdots & \vert 
\end{pmatrix}} \\\\
&=\underbrace{\begin{pmatrix} \vert & \vert & \cdots & \vert \\\\
\mathbf{v}\_1 & \mathbf{v}\_2 & \cdots & \mathbf{v}\_n \\\\
\vert & \vert & \cdots & \vert 
\end{pmatrix}}\_\mathbf{S}
\underbrace{\begin{pmatrix} \lambda_1 &&& \\\\ & \lambda_2 & & \\\\ && \ddots & \\\\ &&& \lambda_n   \end{pmatrix}}\_\mathbf{\Lambda}.
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0A%5Ctextbf%7BA%7D%20%5Cunderbrace%7B%5Cbegin%7Bpmatrix%7D%0A%5Cvert%20%26%20%5Cvert%20%26%20%5Ccdots%20%26%20%5Cvert%20%5C%5C%0A%5Cmathbf%7Bv%7D_1%20%26%20%5Cmathbf%7Bv%7D_2%20%26%20%5Ccdots%20%26%20%5Cmathbf%7Bv%7D_n%20%5C%5C%0A%5Cvert%20%26%20%5Cvert%20%26%20%5Ccdots%20%26%20%5Cvert%20%0A%5Cend%7Bpmatrix%7D%7D_%5Cmathbf%7BS%7D%20%26%3D%20%0A%7B%5Cbegin%7Bpmatrix%7D%20%5Cvert%20%26%20%5Cvert%20%26%20%5Ccdots%20%26%20%5Cvert%20%5C%5C%0A%5Clambda_1%20%5Cmathbf%7Bv%7D_1%20%26%20%5Clambda_2%20%5Cmathbf%7Bv%7D_2%20%26%20%5Ccdots%20%26%20%5Clambda_n%20%5Cmathbf%7Bv%7D_n%20%5C%5C%0A%5Cvert%20%26%20%5Cvert%20%26%20%5Ccdots%20%26%20%5Cvert%20%0A%5Cend%7Bpmatrix%7D%7D%20%5C%5C%0A%26%3D%5Cunderbrace%7B%5Cbegin%7Bpmatrix%7D%20%5Cvert%20%26%20%5Cvert%20%26%20%5Ccdots%20%26%20%5Cvert%20%5C%5C%0A%5Cmathbf%7Bv%7D_1%20%26%20%5Cmathbf%7Bv%7D_2%20%26%20%5Ccdots%20%26%20%5Cmathbf%7Bv%7D_n%20%5C%5C%0A%5Cvert%20%26%20%5Cvert%20%26%20%5Ccdots%20%26%20%5Cvert%20%0A%5Cend%7Bpmatrix%7D%7D_%5Cmathbf%7BS%7D%0A%5Cunderbrace%7B%5Cbegin%7Bpmatrix%7D%20%5Clambda_1%20%26%26%26%20%5C%5C%20%26%20%5Clambda_2%20%26%20%26%20%5C%5C%20%26%26%20%5Cddots%20%26%20%5C%5C%20%26%26%26%20%5Clambda_n%20%20%20%5Cend%7Bpmatrix%7D%7D_%5Cmathbf%7B%5CLambda%7D.%0A%5Cend%7Balign%7D "\begin{align}
\textbf{A} \underbrace{\begin{pmatrix}
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
\end{align}")

Since
![\textbf{A}\mathbf{S} = \mathbf{S} \mathbf{\Lambda}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%5Cmathbf%7BS%7D%20%3D%20%5Cmathbf%7BS%7D%20%5Cmathbf%7B%5CLambda%7D "\textbf{A}\mathbf{S} = \mathbf{S} \mathbf{\Lambda}"),
we can write
![\textbf{A} = \mathbf{S} \mathbf{\Lambda} \mathbf{S}^{-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%3D%20%5Cmathbf%7BS%7D%20%5Cmathbf%7B%5CLambda%7D%20%5Cmathbf%7BS%7D%5E%7B-1%7D "\textbf{A} = \mathbf{S} \mathbf{\Lambda} \mathbf{S}^{-1}").
If we think of
![\mathbf{S}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BS%7D "\mathbf{S}")
as describing a change of basis, this equation says that the action of
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
is like going into another basis, multiplying by a diagonal matrix, and
then changing back to the original basis.

Let’s show that we can do this with our example matrix
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}")
from before,

![\textbf{A} = \begin{pmatrix} -3 & 2 \\\\ 2 & -3 \end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20-3%20%26%202%20%5C%5C%202%20%26%20-3%20%5Cend%7Bpmatrix%7D. "\textbf{A} = \begin{pmatrix} -3 & 2 \\ 2 & -3 \end{pmatrix}.")

``` r
A <- matrix(c(-3,2,2,-3),byrow=TRUE,nrow=2)
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
S%*%Lambda%*%solve(S)
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

![A_0 = 2,\quad A_n = 5A\_{n-1}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A_0%20%3D%202%2C%5Cquad%20A_n%20%3D%205A_%7Bn-1%7D. "A_0 = 2,\quad A_n = 5A_{n-1}.")

What is a formula for
![A_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A_n "A_n")?
Well, notice that

![\begin{align}
A_1 &= 2 \cdot 5\\\\
A_2 &= 2 \cdot 5 \cdot 5\\\\
&\vdots\\\\
A_n &= 2 \cdot 5^n.
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0AA_1%20%26%3D%202%20%5Ccdot%205%5C%5C%0AA_2%20%26%3D%202%20%5Ccdot%205%20%5Ccdot%205%5C%5C%0A%26%5Cvdots%5C%5C%0AA_n%20%26%3D%202%20%5Ccdot%205%5En.%0A%5Cend%7Balign%7D "\begin{align}
A_1 &= 2 \cdot 5\\
A_2 &= 2 \cdot 5 \cdot 5\\
&\vdots\\
A_n &= 2 \cdot 5^n.
\end{align}")

We’ll need this result! Now recall the Fibonacci sequence,

![F_0 = 1,\quad F_1 = 1, \quad F_n = F\_{n-1} + F\_{n-2}, \quad n \geq 2.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_0%20%3D%201%2C%5Cquad%20F_1%20%3D%201%2C%20%5Cquad%20F_n%20%3D%20F_%7Bn-1%7D%20%2B%20F_%7Bn-2%7D%2C%20%5Cquad%20n%20%5Cgeq%202. "F_0 = 1,\quad F_1 = 1, \quad F_n = F_{n-1} + F_{n-2}, \quad n \geq 2.")

If I asked you for
![F\_{1000}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_%7B1000%7D "F_{1000}"),
then using the recursive definition above, you’d need to compute the
![999](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;999 "999")
terms that come before it. But is there a way to avoid doing this? That
is, is there a way to come up with an explicit (non-recursive) formula
for
![F_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_n "F_n")?

So that we can use tools of linear algebra, let’s transform this to a
problem involving a matrix. Define
![G_n = F\_{n-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;G_n%20%3D%20F_%7Bn-1%7D "G_n = F_{n-1}").
Then substituting into the definition of
![F_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_n "F_n"),
we can write
![F_n = F\_{n-1} + G\_{n-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_n%20%3D%20F_%7Bn-1%7D%20%2B%20G_%7Bn-1%7D "F_n = F_{n-1} + G_{n-1}").
Putting it all together,

![\begin{align}
F_n &= F\_{n-1} + G\_{n-1} \\\\
G_n &= F\_{n-1} \\\\
\begin{pmatrix} F_n \\\\ G_n \end{pmatrix} &= \begin{pmatrix} 1 & 1 \\\\ 1 & 0 \end{pmatrix} \begin{pmatrix} F\_{n-1} \\\\ G\_{n-1} \end{pmatrix} \\\\
\mathbf{F}\_n &= \textbf{A} \mathbf{F}\_{n-1}.
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0AF_n%20%26%3D%20F_%7Bn-1%7D%20%2B%20G_%7Bn-1%7D%20%5C%5C%0AG_n%20%26%3D%20F_%7Bn-1%7D%20%5C%5C%0A%5Cbegin%7Bpmatrix%7D%20F_n%20%5C%5C%20G_n%20%5Cend%7Bpmatrix%7D%20%26%3D%20%5Cbegin%7Bpmatrix%7D%201%20%26%201%20%5C%5C%201%20%26%200%20%5Cend%7Bpmatrix%7D%20%5Cbegin%7Bpmatrix%7D%20F_%7Bn-1%7D%20%5C%5C%20G_%7Bn-1%7D%20%5Cend%7Bpmatrix%7D%20%5C%5C%0A%5Cmathbf%7BF%7D_n%20%26%3D%20%5Ctextbf%7BA%7D%20%5Cmathbf%7BF%7D_%7Bn-1%7D.%0A%5Cend%7Balign%7D "\begin{align}
F_n &= F_{n-1} + G_{n-1} \\
G_n &= F_{n-1} \\
\begin{pmatrix} F_n \\ G_n \end{pmatrix} &= \begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix} \begin{pmatrix} F_{n-1} \\ G_{n-1} \end{pmatrix} \\
\mathbf{F}_n &= \textbf{A} \mathbf{F}_{n-1}.
\end{align}")

Let’s diagonalize
![\textbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BA%7D "\textbf{A}"),
and momentarily, you’ll see why.

``` r
A <- matrix(c(1,1,1,0),byrow=TRUE,nrow=2)
lambdap <- (1+sqrt(5))/2
lambdap
```

    ## [1] 1.61803398875

``` r
lambdam <- (1-sqrt(5))/2
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
v1 <- S[,1]/S[2,1]
v2 <- S[,2]/S[2,2]
S <- cbind(v1,v2)
S
```

    ##                 v1             v2
    ## [1,] 1.61803398875 -0.61803398875
    ## [2,] 1.00000000000  1.00000000000

``` r
S%*%Lambda%*%solve(S)
```

    ##      [,1]               [,2]
    ## [1,]    1  1.00000000000e+00
    ## [2,]    1 -1.11022302463e-16

So to recap, we have

![\mathbf{\Lambda} = \begin{pmatrix} \lambda\_+ & 0 \\\\ 0 & \lambda\_- \end{pmatrix}, \quad \mathbf{S} = \begin{pmatrix} \lambda\_+ & \lambda\_- \\\\ 1 & 1 \end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7B%5CLambda%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20%5Clambda_%2B%20%26%200%20%5C%5C%200%20%26%20%5Clambda_-%20%5Cend%7Bpmatrix%7D%2C%20%5Cquad%20%5Cmathbf%7BS%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20%5Clambda_%2B%20%26%20%5Clambda_-%20%5C%5C%201%20%26%201%20%5Cend%7Bpmatrix%7D. "\mathbf{\Lambda} = \begin{pmatrix} \lambda_+ & 0 \\ 0 & \lambda_- \end{pmatrix}, \quad \mathbf{S} = \begin{pmatrix} \lambda_+ & \lambda_- \\ 1 & 1 \end{pmatrix}.")

Why would we do this? Let’s use the diagonalization to re-write our
problem:

![\begin{align}
\mathbf{F}\_n &= \textbf{A} \mathbf{F}\_{n-1}\\\\
\mathbf{F}\_n &= \mathbf{S} \mathbf{\Lambda} \mathbf{S}^{-1} \mathbf{F}\_{n-1}\\\\
\mathbf{S}^{-1} \mathbf{F} &= \mathbf{\Lambda} \mathbf{S}^{-1} \mathbf{F}\_{n-1}.
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0A%5Cmathbf%7BF%7D_n%20%26%3D%20%5Ctextbf%7BA%7D%20%5Cmathbf%7BF%7D_%7Bn-1%7D%5C%5C%0A%5Cmathbf%7BF%7D_n%20%26%3D%20%5Cmathbf%7BS%7D%20%5Cmathbf%7B%5CLambda%7D%20%5Cmathbf%7BS%7D%5E%7B-1%7D%20%5Cmathbf%7BF%7D_%7Bn-1%7D%5C%5C%0A%5Cmathbf%7BS%7D%5E%7B-1%7D%20%5Cmathbf%7BF%7D%20%26%3D%20%5Cmathbf%7B%5CLambda%7D%20%5Cmathbf%7BS%7D%5E%7B-1%7D%20%5Cmathbf%7BF%7D_%7Bn-1%7D.%0A%5Cend%7Balign%7D "\begin{align}
\mathbf{F}_n &= \textbf{A} \mathbf{F}_{n-1}\\
\mathbf{F}_n &= \mathbf{S} \mathbf{\Lambda} \mathbf{S}^{-1} \mathbf{F}_{n-1}\\
\mathbf{S}^{-1} \mathbf{F} &= \mathbf{\Lambda} \mathbf{S}^{-1} \mathbf{F}_{n-1}.
\end{align}")

Now define

![\mathbf{S}^{-1} \mathbf{F} = \mathbf{S}^{-1} \begin{pmatrix} F_n \\\\ G_n \end{pmatrix} \equiv  \begin{pmatrix} A_n \\\\ B_n \end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BS%7D%5E%7B-1%7D%20%5Cmathbf%7BF%7D%20%3D%20%5Cmathbf%7BS%7D%5E%7B-1%7D%20%5Cbegin%7Bpmatrix%7D%20F_n%20%5C%5C%20G_n%20%5Cend%7Bpmatrix%7D%20%5Cequiv%20%20%5Cbegin%7Bpmatrix%7D%20A_n%20%5C%5C%20B_n%20%5Cend%7Bpmatrix%7D. "\mathbf{S}^{-1} \mathbf{F} = \mathbf{S}^{-1} \begin{pmatrix} F_n \\ G_n \end{pmatrix} \equiv  \begin{pmatrix} A_n \\ B_n \end{pmatrix}.")

Then our problem becomes

![\begin{align}
\begin{pmatrix} A_n \\\\ B_n \end{pmatrix} &= \mathbf{\Lambda} \begin{pmatrix} A\_{n-1} \\\\ B\_{n-1} \end{pmatrix} \\\\
\begin{pmatrix} A_n \\\\ B_n \end{pmatrix} &= \begin{pmatrix} \lambda\_+ &  0 \\\\ 0 & \lambda\_- \end{pmatrix} \begin{pmatrix} A\_{n-1} \\\\ B\_{n-1}. \end{pmatrix}
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0A%5Cbegin%7Bpmatrix%7D%20A_n%20%5C%5C%20B_n%20%5Cend%7Bpmatrix%7D%20%26%3D%20%5Cmathbf%7B%5CLambda%7D%20%5Cbegin%7Bpmatrix%7D%20A_%7Bn-1%7D%20%5C%5C%20B_%7Bn-1%7D%20%5Cend%7Bpmatrix%7D%20%5C%5C%0A%5Cbegin%7Bpmatrix%7D%20A_n%20%5C%5C%20B_n%20%5Cend%7Bpmatrix%7D%20%26%3D%20%5Cbegin%7Bpmatrix%7D%20%5Clambda_%2B%20%26%20%200%20%5C%5C%200%20%26%20%5Clambda_-%20%5Cend%7Bpmatrix%7D%20%5Cbegin%7Bpmatrix%7D%20A_%7Bn-1%7D%20%5C%5C%20B_%7Bn-1%7D.%20%5Cend%7Bpmatrix%7D%0A%5Cend%7Balign%7D "\begin{align}
\begin{pmatrix} A_n \\ B_n \end{pmatrix} &= \mathbf{\Lambda} \begin{pmatrix} A_{n-1} \\ B_{n-1} \end{pmatrix} \\
\begin{pmatrix} A_n \\ B_n \end{pmatrix} &= \begin{pmatrix} \lambda_+ &  0 \\ 0 & \lambda_- \end{pmatrix} \begin{pmatrix} A_{n-1} \\ B_{n-1}. \end{pmatrix}
\end{align}")

Because this problem is diagonal, the top and bottom equations are
uncoupled, which makes them much easier to solve. In fact, they are
geometrtic, just like the warm-up problem we did. We can write

![\begin{pmatrix} A_n \\\\ B_n \end{pmatrix} = \begin{pmatrix}  A_0 \lambda\_+^n \\\\ B_0 \lambda\_-^n \end{pmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Bpmatrix%7D%20A_n%20%5C%5C%20B_n%20%5Cend%7Bpmatrix%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20%20A_0%20%5Clambda_%2B%5En%20%5C%5C%20B_0%20%5Clambda_-%5En%20%5Cend%7Bpmatrix%7D. "\begin{pmatrix} A_n \\ B_n \end{pmatrix} = \begin{pmatrix}  A_0 \lambda_+^n \\ B_0 \lambda_-^n \end{pmatrix}.")

But we didn’t want to know
![A_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A_n "A_n")
and
![B_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;B_n "B_n").
We wanted to know
![F_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_n "F_n")
and
![G_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;G_n "G_n").
So we have to undo the transformation we did before. We have

![\begin{align}
\mathbf{S}^{-1} \mathbf{F} &= \begin{pmatrix} A_n \\\\ B_n \end{pmatrix} \\\\
\mathbf{S}^{-1} \mathbf{F} &= \begin{pmatrix}  A_0 \lambda_1^n \\\\ B_0 \lambda_2^n \end{pmatrix} \\\\
\mathbf{F} &= \mathbf{S} \begin{pmatrix}  A_0 \lambda_1^n \\\\ B_0 \lambda_2^n \end{pmatrix} \\\\
\mathbf{F} &= \begin{pmatrix} \lambda_1 & \lambda_2 \\\\ 1 & 1 \end{pmatrix} \begin{pmatrix}  A_0 \lambda_1^n \\\\ B_0 \lambda_2^n \end{pmatrix} \\\\
\mathbf{F} &= \begin{pmatrix} \lambda_1 A_0 \lambda_1^n + \lambda_2 B_0 \lambda^2n \\\\ \mathrm{doesn't\\ matter}\end{pmatrix}.
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0A%5Cmathbf%7BS%7D%5E%7B-1%7D%20%5Cmathbf%7BF%7D%20%26%3D%20%5Cbegin%7Bpmatrix%7D%20A_n%20%5C%5C%20B_n%20%5Cend%7Bpmatrix%7D%20%5C%5C%0A%5Cmathbf%7BS%7D%5E%7B-1%7D%20%5Cmathbf%7BF%7D%20%26%3D%20%5Cbegin%7Bpmatrix%7D%20%20A_0%20%5Clambda_1%5En%20%5C%5C%20B_0%20%5Clambda_2%5En%20%5Cend%7Bpmatrix%7D%20%5C%5C%0A%5Cmathbf%7BF%7D%20%26%3D%20%5Cmathbf%7BS%7D%20%5Cbegin%7Bpmatrix%7D%20%20A_0%20%5Clambda_1%5En%20%5C%5C%20B_0%20%5Clambda_2%5En%20%5Cend%7Bpmatrix%7D%20%5C%5C%0A%5Cmathbf%7BF%7D%20%26%3D%20%5Cbegin%7Bpmatrix%7D%20%5Clambda_1%20%26%20%5Clambda_2%20%5C%5C%201%20%26%201%20%5Cend%7Bpmatrix%7D%20%5Cbegin%7Bpmatrix%7D%20%20A_0%20%5Clambda_1%5En%20%5C%5C%20B_0%20%5Clambda_2%5En%20%5Cend%7Bpmatrix%7D%20%5C%5C%0A%5Cmathbf%7BF%7D%20%26%3D%20%5Cbegin%7Bpmatrix%7D%20%5Clambda_1%20A_0%20%5Clambda_1%5En%20%2B%20%5Clambda_2%20B_0%20%5Clambda%5E2n%20%5C%5C%20%5Cmathrm%7Bdoesn%27t%5C%20matter%7D%5Cend%7Bpmatrix%7D.%0A%5Cend%7Balign%7D "\begin{align}
\mathbf{S}^{-1} \mathbf{F} &= \begin{pmatrix} A_n \\ B_n \end{pmatrix} \\
\mathbf{S}^{-1} \mathbf{F} &= \begin{pmatrix}  A_0 \lambda_1^n \\ B_0 \lambda_2^n \end{pmatrix} \\
\mathbf{F} &= \mathbf{S} \begin{pmatrix}  A_0 \lambda_1^n \\ B_0 \lambda_2^n \end{pmatrix} \\
\mathbf{F} &= \begin{pmatrix} \lambda_1 & \lambda_2 \\ 1 & 1 \end{pmatrix} \begin{pmatrix}  A_0 \lambda_1^n \\ B_0 \lambda_2^n \end{pmatrix} \\
\mathbf{F} &= \begin{pmatrix} \lambda_1 A_0 \lambda_1^n + \lambda_2 B_0 \lambda^2n \\ \mathrm{doesn't\ matter}\end{pmatrix}.
\end{align}")

Therefore,

![F_n = A_0 \lambda\_+^{n+1} + B_0 \lambda\_-^{n+1}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_n%20%3D%20A_0%20%5Clambda_%2B%5E%7Bn%2B1%7D%20%2B%20B_0%20%5Clambda_-%5E%7Bn%2B1%7D. "F_n = A_0 \lambda_+^{n+1} + B_0 \lambda_-^{n+1}.")

To find the unknown constants, we plug in the initial conditions.

![\begin{align}
F_0 &= A_0 \lambda\_+ + B_0 \lambda\_- = 1\\\\
F_1 &= A_0 \lambda\_+^2 + B_0 \lambda\_-^2 = 1.
\end{align}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Balign%7D%0AF_0%20%26%3D%20A_0%20%5Clambda_%2B%20%2B%20B_0%20%5Clambda_-%20%3D%201%5C%5C%0AF_1%20%26%3D%20A_0%20%5Clambda_%2B%5E2%20%2B%20B_0%20%5Clambda_-%5E2%20%3D%201.%0A%5Cend%7Balign%7D "\begin{align}
F_0 &= A_0 \lambda_+ + B_0 \lambda_- = 1\\
F_1 &= A_0 \lambda_+^2 + B_0 \lambda_-^2 = 1.
\end{align}")

Solving, we find
![A_0 = 1/\sqrt{5}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A_0%20%3D%201%2F%5Csqrt%7B5%7D "A_0 = 1/\sqrt{5}"),
![B_0 = -1/\sqrt{5}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;B_0%20%3D%20-1%2F%5Csqrt%7B5%7D "B_0 = -1/\sqrt{5}").
The final answer, then, is

![F_n = \frac{1}{\sqrt{5}}\lambda\_+^{n+1} - \frac{1}{\sqrt{5}}\lambda\_-^{n+1},\quad \lambda\_\pm = \frac{1 \pm \sqrt{5}}{2}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_n%20%3D%20%5Cfrac%7B1%7D%7B%5Csqrt%7B5%7D%7D%5Clambda_%2B%5E%7Bn%2B1%7D%20-%20%5Cfrac%7B1%7D%7B%5Csqrt%7B5%7D%7D%5Clambda_-%5E%7Bn%2B1%7D%2C%5Cquad%20%5Clambda_%5Cpm%20%3D%20%5Cfrac%7B1%20%5Cpm%20%5Csqrt%7B5%7D%7D%7B2%7D. "F_n = \frac{1}{\sqrt{5}}\lambda_+^{n+1} - \frac{1}{\sqrt{5}}\lambda_-^{n+1},\quad \lambda_\pm = \frac{1 \pm \sqrt{5}}{2}.")

Let’s test this out.

``` r
fib <- function(n){
  lambdap <- (1+sqrt(5))/2
  lambdam <- (1-sqrt(5))/2
  Fn <- 1/sqrt(5)*lambdap^(n+1) - 1/sqrt(5)*lambdam^(n+1)
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
![F_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_n "F_n")
for large
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
in a simpler way. Since
![\|\lambda\_-\| \< 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%5Clambda_-%7C%20%3C%201 "|\lambda_-| < 1"),
after many repeated iterations, the term involving
![\lambda\_-](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_- "\lambda_-")
will die out and we can write, for large
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n"),

![F_n \approx \frac{1}{\sqrt{5}}\lambda\_+^{n+1}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_n%20%5Capprox%20%5Cfrac%7B1%7D%7B%5Csqrt%7B5%7D%7D%5Clambda_%2B%5E%7Bn%2B1%7D. "F_n \approx \frac{1}{\sqrt{5}}\lambda_+^{n+1}.")

Let’s test this out!

``` r
fibapprox <- function(n){
  lambdap <- (1+sqrt(5))/2
  Fn <- 1/sqrt(5)*lambdap^(n+1)
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
error <- abs(fib(n)-fibapprox(n))
plot(n,log10(error))
```

![](coursenotes_files/figure-gfm/unnamed-chunk-66-1.png)<!-- -->

It’s important to remember from linear algebra that not every matrix can
be diagonalized. For a matrix to be diagonalizable, you need for the
multiplicity of each eigenvalue to be the same as the dimension of the
eigenspace. Here are a couple of examples.

``` r
A <- matrix(c(5,-4,4,12,-11,12,4,-4,5),byrow=TRUE,nrow=3)
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

Here, the **algebraic multiplicity** of
![-3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-3 "-3")
is 1 because it is only an eigenvalue once, and the **geometric
multiplicity** is 1 because it only has one eigenvector. The algebraic
multiplicity of
![1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1 "1")
is
![2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2 "2")
because it is an eigenvalue twice. Since it has two independent
eigenvectors, the geometric multiplicity is also 2. Therefore, this
matrix is diagonalizable. Stated differently: we need enough independent
eigenvectors to form the
![\mathbf{S}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BS%7D "\mathbf{S}")
matrix.

By way of counterexample, consider this problem.

``` r
A <- matrix(c(1,-1,0,1),byrow=TRUE,nrow=2)
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

The eigenvalue
![1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1 "1")
has algebraic multiplicity 2, but there is only one eigenvector (the
trivial eigenvector doesn’t count) so it has geometric multiplicity 1.
We don’t have enough eigenvectors to make
![\mathbf{S}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BS%7D "\mathbf{S}"),
so the matrix is not diagonalizable.

## Power iteration

Please find the eigenvalues of this matrix:

``` r
set.seed(123)
A <- matrix(sample(-100:100,64),nrow=8)
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
A <- matrix(c(-13,170,240,19,-224,-320,-14,166,237),byrow=TRUE,nrow=3)
A
```

    ##      [,1] [,2] [,3]
    ## [1,]  -13  170  240
    ## [2,]   19 -224 -320
    ## [3,]  -14  166  237

``` r
v <- c(1,1,1)
v <- A%*%v
v
```

    ##      [,1]
    ## [1,]  397
    ## [2,] -525
    ## [3,]  389

``` r
v <- v/Norm(v,2)
v
```

    ##                 [,1]
    ## [1,]  0.519251568355
    ## [2,] -0.686667691150
    ## [3,]  0.508788060681

``` r
v <- A%*%v
v
```

    ##                 [,1]
    ## [1,] -1.374643320759
    ## [2,]  0.867163198538
    ## [3,] -0.673588306557

``` r
v <- v/Norm(v,2)
for (i in 1:100){
  v <- A%*%v
  v <- v/Norm(v,2)
}
(A%*%v)/v
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

Interesting. Using this iteration, we found the eigenvalue
![-3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-3 "-3").
Any guesses why?

Let’s try another example.

``` r
A <- matrix(c(2,-520,8,-90,468,-360,40,-698,160),byrow=TRUE,nrow=3)
v <- c(1,1,1)
for (i in 1:100){
  v <- A%*%v
  v <- v/Norm(v,2)
}
(A%*%v)/v
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

This method if called **power iteration**. How and why, exactly, does it
work though? This is the subject of your activity for today.

You might not realize it, but you use power iteration nearly every day.
Google searches are based on an algorithm called PageRank, which is
applies power iteration to a matrix that encodes links between web
pages. We may investigate this in class at some point. But as you can
imagine, the matrix in question is gigantic, which is why numerical
methods are necessary.
