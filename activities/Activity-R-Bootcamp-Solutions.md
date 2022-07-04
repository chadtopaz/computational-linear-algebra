Activity - R Bootcamp
================
Solutions

# Objectives

-   Create your first in-class R Markdown
-   Create matrices in R
-   Evaluate polynomials in R

# Problem 1

Generate the following matrices with as little code as possible.

1.  

![\\begin{bmatrix}
1 & 1 & 1 & 1 & 1 & 1 & 1 \\\\
1 & 1 & 1 & 1 & 1 & 1 & 1 \\\\
1 & 1 & 1 & 1 & 1 & 1 & 1 \\\\
1 & 1 & 1 & 1 & 1 & 1 & 1 \\\\
1 & 1 & 1 & 1 & 1 & 1 & 1 \\\\
1 & 1 & 1 & 5 & 5 & 5 & 1 \\\\
1 & 1 & 1 & 5 & 5 & 5 & 1 \\\\
1 & 1 & 1 & 1 & 1 & 1 & 1 
\\end{bmatrix}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Bbmatrix%7D%0A1%20%26%201%20%26%201%20%26%201%20%26%201%20%26%201%20%26%201%20%5C%5C%0A1%20%26%201%20%26%201%20%26%201%20%26%201%20%26%201%20%26%201%20%5C%5C%0A1%20%26%201%20%26%201%20%26%201%20%26%201%20%26%201%20%26%201%20%5C%5C%0A1%20%26%201%20%26%201%20%26%201%20%26%201%20%26%201%20%26%201%20%5C%5C%0A1%20%26%201%20%26%201%20%26%201%20%26%201%20%26%201%20%26%201%20%5C%5C%0A1%20%26%201%20%26%201%20%26%205%20%26%205%20%26%205%20%26%201%20%5C%5C%0A1%20%26%201%20%26%201%20%26%205%20%26%205%20%26%205%20%26%201%20%5C%5C%0A1%20%26%201%20%26%201%20%26%201%20%26%201%20%26%201%20%26%201%20%0A%5Cend%7Bbmatrix%7D "\begin{bmatrix}
1 & 1 & 1 & 1 & 1 & 1 & 1 \\
1 & 1 & 1 & 1 & 1 & 1 & 1 \\
1 & 1 & 1 & 1 & 1 & 1 & 1 \\
1 & 1 & 1 & 1 & 1 & 1 & 1 \\
1 & 1 & 1 & 1 & 1 & 1 & 1 \\
1 & 1 & 1 & 5 & 5 & 5 & 1 \\
1 & 1 & 1 & 5 & 5 & 5 & 1 \\
1 & 1 & 1 & 1 & 1 & 1 & 1 
\end{bmatrix}")

2.  

![\\begin{bmatrix}
1 & 2 & 3 & \\cdots & 10 \\\\
2 & 4 & 6 & \\cdots & 20 \\\\
3 & 6 & 9 & \\cdots & 30 \\\\
\\vdots & \\vdots & \\vdots & \\ddots & \\vdots \\\\
8 & 16 & 24 & \\cdots & 80
\\end{bmatrix}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Bbmatrix%7D%0A1%20%26%202%20%26%203%20%26%20%5Ccdots%20%26%2010%20%5C%5C%0A2%20%26%204%20%26%206%20%26%20%5Ccdots%20%26%2020%20%5C%5C%0A3%20%26%206%20%26%209%20%26%20%5Ccdots%20%26%2030%20%5C%5C%0A%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cddots%20%26%20%5Cvdots%20%5C%5C%0A8%20%26%2016%20%26%2024%20%26%20%5Ccdots%20%26%2080%0A%5Cend%7Bbmatrix%7D%0A "\begin{bmatrix}
1 & 2 & 3 & \cdots & 10 \\
2 & 4 & 6 & \cdots & 20 \\
3 & 6 & 9 & \cdots & 30 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
8 & 16 & 24 & \cdots & 80
\end{bmatrix}
")

3.  

![\\begin{bmatrix}
1 & 1 & 1 & 1 & 1 \\\\
4 & 4 & 4 & 4 & 4 \\\\
9 & 9 & 9 & 9 & 9 \\\\
\\vdots & \\vdots & \\vdots & \\vdots & \\vdots \\\\
121 & 121 & 121 & 121 & 121
\\end{bmatrix}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Bbmatrix%7D%0A1%20%26%201%20%26%201%20%26%201%20%26%201%20%5C%5C%0A4%20%26%204%20%26%204%20%26%204%20%26%204%20%5C%5C%0A9%20%26%209%20%26%209%20%26%209%20%26%209%20%5C%5C%0A%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%5C%5C%0A121%20%26%20121%20%26%20121%20%26%20121%20%26%20121%0A%5Cend%7Bbmatrix%7D%0A "\begin{bmatrix}
1 & 1 & 1 & 1 & 1 \\
4 & 4 & 4 & 4 & 4 \\
9 & 9 & 9 & 9 & 9 \\
\vdots & \vdots & \vdots & \vdots & \vdots \\
121 & 121 & 121 & 121 & 121
\end{bmatrix}
")

# Problem 1 Solution

1.  

``` r
A <- matrix(1, nrow = 8, ncol = 7)
A[6:7,4:6] <- 5
print(A)
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
    ## [1,]    1    1    1    1    1    1    1
    ## [2,]    1    1    1    1    1    1    1
    ## [3,]    1    1    1    1    1    1    1
    ## [4,]    1    1    1    1    1    1    1
    ## [5,]    1    1    1    1    1    1    1
    ## [6,]    1    1    1    5    5    5    1
    ## [7,]    1    1    1    5    5    5    1
    ## [8,]    1    1    1    1    1    1    1

2.  

``` r
v1 <- 1:8
v2 <- 1:10
B <- v1 %*% t(v2) # t() is the transpose operator
print(B)
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
    ## [1,]    1    2    3    4    5    6    7    8    9    10
    ## [2,]    2    4    6    8   10   12   14   16   18    20
    ## [3,]    3    6    9   12   15   18   21   24   27    30
    ## [4,]    4    8   12   16   20   24   28   32   36    40
    ## [5,]    5   10   15   20   25   30   35   40   45    50
    ## [6,]    6   12   18   24   30   36   42   48   54    60
    ## [7,]    7   14   21   28   35   42   49   56   63    70
    ## [8,]    8   16   24   32   40   48   56   64   72    80

A shortcut is to use the outer product notation %o%:

``` r
B <- v1 %o% v2
```

3.  If we think about matrix multiplication as taking different linear
    combinations of the columns, we can see that each of the columns of
    this matrix is a linear combination (in fact, the same linear
    combination) of a single column vector. Note also how we can take an
    exponent of each of the elements of a vector in R and how we can
    create a vector with the repeated entries.

``` r
w1 <- 1:11
w2 <- rep(1,5)
C <- w1^2 %o% w2
```

# Problem 2

Here is an example of how to create a function in R:

``` r
P <- function(x){
  y <- 2*x^5 - 3*x^4 + x^3 - 10*x^2 + 1
  return(y)
}
```

I can evaluate the function at a specific value of
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x"):

``` r
P(2)
```

    ## [1] -15

or at multiple values of
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x"):

``` r
xx <- 3:7
P(xx)
```

    ## [1]   181  1185  4251 11521 26265

1.  Write a function to evaluate the same polynomial,
    ![P(x)=2x^5-3x^4+x^3-10x^2+1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;P%28x%29%3D2x%5E5-3x%5E4%2Bx%5E3-10x%5E2%2B1 "P(x)=2x^5-3x^4+x^3-10x^2+1"),
    but this time, let the inputs to the function be a vector containing
    the coefficients of
    ![P](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;P "P")
    and then a vector of values of
    ![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
    at which you want to evaluate the polynomial. It should look like

``` r
P <- function(a,x) {
  
  # Your code goes here
  
}
```

where I have used
![a](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a "a")
to represent the vector of coefficients in order of highest degree to
lowest degree (making sure to account for terms that have coefficients
of zero). The point of doing this is that now don’t need to change our
function each time we want to evaluate a different polynomial. We can
just pass in a different set of coefficients. **Before you start
coding**, make a plan (in words) for how to do the calculation. Report
your plan, write the code, use it to evaluate the polynomial at
![x = 3, 4, 5, 6, 7](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%20%3D%203%2C%204%2C%205%2C%206%2C%207 "x = 3, 4, 5, 6, 7"),
and verify that the answers match what we computed above

Hints:

-   Actually create a plan first. The real intellectual challenge of the
    problem is figuring out how to do the calculation in a nice, tight,
    efficient way. Turning the idea into code is only the second most
    important part.
-   Avoid using a for-loop at all costs.
-   The command `outer` with `FUN = "^"` will be helpful.
-   My function to evaluate
    ![P(x)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;P%28x%29 "P(x)")
    was under 10 lines of brief code (not including comments).

# Problem 2 Solution

``` r
P <- function(a,x) {
  
  # See how many terms there are in our polynomial
  numterms <- length(a)
  
  # See how many values of x have been input for evaluation
  numx <- length(x)
  
  # Since we will carry out the calculation using matrices,
  # create matrix with one copy of the coefficients list in each
  # row, with one (identical) row for each value of x
  A <- matrix(rep(a,numx), nrow = numx, byrow = TRUE)
  
  # Create vector of powers associated to each term of polynomial
  p <- (numterms-1):0
  
  # Evaluate each value of x raised to each power p
  xp <- outer(x, p, `^`)

  # Multiply xp by coefficients and then
  # take row sums to obtain final answer
  y <- rowSums(A*xp)
  return(y)
  
}

mycoeffs <- c(2, -3, 1, -10, 0, 1)
x <- 3:7
P(mycoeffs,x)
```

    ## [1]   181  1185  4251 11521 26265

# Key Commands

-   `matrix` to create a matrix
-   `:` (colon operator) to create a vector of integers
-   `function` to define a function
-   `length` to get length of a vector
-   `rep` to copy a number/vector multiple times
-   `t` to take transpose of a matrix or vector
-   `%o%` to take outer product of two vectors
-   `outer` to take outer product of two vectors (more flexible command)
-   `rowSums` to get the sum of each row of a matrix
