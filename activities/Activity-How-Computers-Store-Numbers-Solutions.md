Activity - How Computers Store Numbers
================
Solutions

### Objectives

-   Practice with binary and machine representation of numbers
-   Recognize situations when the limitations of floating point
    arithmetic are realized
-   Reformulate those situations to get improved numerical accuracy

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

Convert the binary number 110111.001 to base 10.

### Problem 1 Solution

![1/8 + 1 + 2 + 4 + 16 + 32 = 55.125](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%2F8%20%2B%201%20%2B%202%20%2B%204%20%2B%2016%20%2B%2032%20%3D%2055.125 "1/8 + 1 + 2 + 4 + 16 + 32 = 55.125")

### Problem 2

Convert the base 10 number 1/3 to binary.

### Problem 2 Solution

![1/3 = 1/4 + 1/16 + 1/64 + \\ldots](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%2F3%20%3D%201%2F4%20%2B%201%2F16%20%2B%201%2F64%20%2B%20%5Cldots "1/3 = 1/4 + 1/16 + 1/64 + \ldots")

which in binary is
![0.0101010101\\dots](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.0101010101%5Cdots "0.0101010101\dots").

### Problem 3

By hand, calculate the IEEE double precision representation fl(1/3) and
find the relative error in the computer storage of this number.

### Problem 3 Solution

The sign is
![s=0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s%3D0 "s=0"),
the exponent is
![1021 - 1023 = -2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1021%20-%201023%20%3D%20-2 "1021 - 1023 = -2")
and the mantissa is
![1.0101010101](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1.0101010101 "1.0101010101")
where the
![01](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;01 "01")
pattern after the decimal point is repeated 26 times. This mantissa is

![1 + 1/4 + (1/4)^2 + \\ldots + (1/4)^26 = \\frac{6004799503160661}{4503599627370496}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%20%2B%201%2F4%20%2B%20%281%2F4%29%5E2%20%2B%20%5Cldots%20%2B%20%281%2F4%29%5E26%20%3D%20%5Cfrac%7B6004799503160661%7D%7B4503599627370496%7D "1 + 1/4 + (1/4)^2 + \ldots + (1/4)^26 = \frac{6004799503160661}{4503599627370496}")

So, our final machine number is this numnber times
![1/4](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%2F4 "1/4"),
which is

![fl(1/3) = \\frac{6004799503160661}{18014398509481984}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;fl%281%2F3%29%20%3D%20%5Cfrac%7B6004799503160661%7D%7B18014398509481984%7D "fl(1/3) = \frac{6004799503160661}{18014398509481984}")

The relative error is

![\\frac{1/3 -fl(1/3)}{1/3}= \\frac{1}{54043195528445952} \\approx 1.85\\times10^{-17}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B1%2F3%20-fl%281%2F3%29%7D%7B1%2F3%7D%3D%20%5Cfrac%7B1%7D%7B54043195528445952%7D%20%5Capprox%201.85%5Ctimes10%5E%7B-17%7D "\frac{1/3 -fl(1/3)}{1/3}= \frac{1}{54043195528445952} \approx 1.85\times10^{-17}")

### Problem 4

Consider the function

![f(x)=\\frac{1-(1-x)^3}{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28x%29%3D%5Cfrac%7B1-%281-x%29%5E3%7D%7Bx%7D "f(x)=\frac{1-(1-x)^3}{x}")

for x =
![10^{-1}, 10^{-2}, \\ldots, 10^{-14}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;10%5E%7B-1%7D%2C%2010%5E%7B-2%7D%2C%20%5Cldots%2C%2010%5E%7B-14%7D "10^{-1}, 10^{-2}, \ldots, 10^{-14}").
Then use an alternative form of the expression that doesn’t suffer from
subtracting nearly equal numbers. Consider the results from the second
calculation to be exact, and make a log-log plot of the relative error
in the first calculation as a function of x. What is the general trend?

### Problem 4 Solution

``` r
f1 <- function(x) {(1-(1-x)^3)/x}
f2 <- function(x) {x^2 - 3*x + 3}
x <- 10^seq(-1,-14)
y1 <- f1(x)
y2 <- f2(x)
relerr <- (abs(y2-y1)/y2)
plot(x,relerr,log="xy")
```

![](Activity-How-Computers-Store-Numbers-Solutions_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

The relative error increases drastically as x decreases.

### Problem 5

Consider a right triangle whose legs are of length 3344556600 and
1.2222222. How much longer is the hypotenuse than the longer leg?

### Problem 5 Solution

We are tempted to compute this, which does not end well:

``` r
a <- 3344556600
b <- 1.2222222
sqrt(a^2+b^2)-a
```

    ## [1] 0

The problem is due to the subtraction of nearly equal numbers. To avoid
this, we can rationalize the expression,

![\\left(\\sqrt{a^2+b^2}-a\\right) \\times \\frac{\\sqrt{a^2+b^2}+a}{\\sqrt{a^2+b^2}+a}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28%5Csqrt%7Ba%5E2%2Bb%5E2%7D-a%5Cright%29%20%5Ctimes%20%5Cfrac%7B%5Csqrt%7Ba%5E2%2Bb%5E2%7D%2Ba%7D%7B%5Csqrt%7Ba%5E2%2Bb%5E2%7D%2Ba%7D "\left(\sqrt{a^2+b^2}-a\right) \times \frac{\sqrt{a^2+b^2}+a}{\sqrt{a^2+b^2}+a}")

which simplifies to

![\\frac{b^2}{\\sqrt{a^2+b^2}+a}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7Bb%5E2%7D%7B%5Csqrt%7Ba%5E2%2Bb%5E2%7D%2Ba%7D. "\frac{b^2}{\sqrt{a^2+b^2}+a}.")

We compute as

``` r
b^2/(sqrt(a^2+b^2)+a)
```

    ## [1] 2.233221e-10
