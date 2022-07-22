Remember to show your work / give your reasoning, as applicable.

    library(pracma)

### Problem 1

Consider the points (-2,0), (-1,14), (2,-4), (3,10).

a\. Find the interpolating polynomial in R, using the Vandermonde
method. Also, evaluate your polynomial at the point x = 1.5 by using the
build-in `horner` command. Make sure you pay attention to what order
your coefficients are in after solving the Vandermonde problem, and what
order the `horner` function expects them it.

b\. Find the interpolating polynomial by hand, using the Lagrange
method. Simplify your answer to show that it is the same polynomial you
obtained in part a. You can use WolframAlpha or a similar tool to help
you simplify it.

c\. Write code to perform interpolation using the Lagrange method. Your
code should take the following as inputs:

-   a vector `xsample` containing coordinates of the interpolation
    points
-   a vector `ysample` containing the y coordinates of the interpolation
    points
-   a vector `xgrid` containing the x coordinates of the points at which
    you want to evaluate the interpolating polynomial

The code should output a vector the values of the Lagrange interpolating
polynomial at the points in `xgrid`.

Plot the interpolating polynomial using 30 equally spaced points from x
= -2 to x = 3. Make the polynomial appear as a curve and add the
original four sampled values as points.

### Problem 1 solution

a\. Your solution goes here.

b\. Your solution goes here.

c\. Your solution goes here.

### Problem 2

a\. Construct the interpolating polynomial *P*(*x*) for the function
*f*(*x*) ≡ 1/(1+*x*<sup>2</sup>) on \[−2,2\] using 100 equally spaced
points. Make a plot of the original function in blue (using 1000 equally
spaced points on the interval) and the interpolating polynomial in red
on the same axes, setting the vertical range of the plot to be \[0,1\].
State the (approximate) infinity norm of the error, that is,
||*f*(*x*) − *P*(*x*)||<sub>∞</sub>.

b\. Repeat part (a) but using 100 Chebyshev nodes.

### Problem 2 Solution

a\. Your solution goes here.

b\. Your solution goes here.
