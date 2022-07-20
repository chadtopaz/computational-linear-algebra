Pset - How Computers Store Numbers Solutions
================

### Problem 1

Convert the following base-10 numbers to binary.

1.  11.25

2.  2/3

3.  99.9

### Problem 1 Solution

1.  11.25 = 11 + 0.25. Integer part:
    ![(11)\_{10} = (8 + 2 + 1)\_{10} = (1011)\_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%2811%29_%7B10%7D%20%3D%20%288%20%2B%202%20%2B%201%29_%7B10%7D%20%3D%20%281011%29_2 "(11)_{10} = (8 + 2 + 1)_{10} = (1011)_2").
    Fractional part:
    ![(0.25)\_{10} = (0.01)\_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%280.25%29_%7B10%7D%20%3D%20%280.01%29_2 "(0.25)_{10} = (0.01)_2"),
    so
    ![(11.25)\_{10} = (1011.01)\_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%2811.25%29_%7B10%7D%20%3D%20%281011.01%29_2 "(11.25)_{10} = (1011.01)_2").

2.  We adopt the strategy of successively multiplying by two and pulling
    off a bit.

    ![\frac{2}{3}\*2 = \frac{1}{3} + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B2%7D%7B3%7D%2A2%20%3D%20%5Cfrac%7B1%7D%7B3%7D%20%2B%201 "\frac{2}{3}*2 = \frac{1}{3} + 1")

    ![\frac{1}{3}\*2 = \frac{2}{3} + 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B1%7D%7B3%7D%2A2%20%3D%20%5Cfrac%7B2%7D%7B3%7D%20%2B%200 "\frac{1}{3}*2 = \frac{2}{3} + 0")

    ![\frac{2}{3}\*2 = \frac{1}{3} + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B2%7D%7B3%7D%2A2%20%3D%20%5Cfrac%7B1%7D%7B3%7D%20%2B%201 "\frac{2}{3}*2 = \frac{1}{3} + 1")

    ![\vdots](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cvdots "\vdots")

    So
    ![(\frac{2}{3})\_{10}=(0.\overline{10})\_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28%5Cfrac%7B2%7D%7B3%7D%29_%7B10%7D%3D%280.%5Coverline%7B10%7D%29_2 "(\frac{2}{3})_{10}=(0.\overline{10})_2")

3.  ![(99.9)\_{10} = (99)\_{10} + (0.9)\_{10}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%2899.9%29_%7B10%7D%20%3D%20%2899%29_%7B10%7D%20%2B%20%280.9%29_%7B10%7D "(99.9)_{10} = (99)_{10} + (0.9)_{10}").
    For the first part,
    ![(99)\_{10} = (64 + 32 + 2 + 1)\_{10} = (1100011)\_{2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%2899%29_%7B10%7D%20%3D%20%2864%20%2B%2032%20%2B%202%20%2B%201%29_%7B10%7D%20%3D%20%281100011%29_%7B2%7D "(99)_{10} = (64 + 32 + 2 + 1)_{10} = (1100011)_{2}").
    For the fractional part, we adopt the strategy of successively
    multiplying by two and pulling off a bit.

![0.9\*2 = 0.8 + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.9%2A2%20%3D%200.8%20%2B%201 "0.9*2 = 0.8 + 1")

![0.8\*2 = 0.6 + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.8%2A2%20%3D%200.6%20%2B%201 "0.8*2 = 0.6 + 1")

![0.6\*2 = 0.2 + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.6%2A2%20%3D%200.2%20%2B%201 "0.6*2 = 0.2 + 1")

![0.2\*2 = 0.4 + 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.2%2A2%20%3D%200.4%20%2B%200 "0.2*2 = 0.4 + 0")

![0.4\*2 = 0.8 + 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.4%2A2%20%3D%200.8%20%2B%200 "0.4*2 = 0.8 + 0")

![0.8\*2 = 0.6 + 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.8%2A2%20%3D%200.6%20%2B%201 "0.8*2 = 0.6 + 1")

![\vdots](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cvdots "\vdots")

So
![(0.9)\_{10}=(0.1\overline{1100})\_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%280.9%29_%7B10%7D%3D%280.1%5Coverline%7B1100%7D%29_2 "(0.9)_{10}=(0.1\overline{1100})_2"),
and therefore
![(99.9)\_{10}=(1100011.1\overline{1100})\_{2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%2899.9%29_%7B10%7D%3D%281100011.1%5Coverline%7B1100%7D%29_%7B2%7D "(99.9)_{10}=(1100011.1\overline{1100})_{2}")

### Problem 2

1.  Explain why you can determine machine epsilon on a computer using
    IEEE double precision and the IEEE Rounding to Nearest Rule by
    calculating (7/3 − 4/3) − 1. I’m not asking for (just) a conceptual
    explanation; you should actually calculate some things on paper to
    make your argument.

2.  Does (4/3 − 1/3) − 1 also give
    ![\epsilon\_{mach}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7Bmach%7D "\epsilon_{mach}")?
    Explain by converting to floating point numbers and carrying out the
    machine arithmetic.

### Problem 2 Solution

1.  Let’s convert the first two terms:
    ![\left(\frac{7}{3}\right)\_{10} = 1.00\overline{10}\*2^1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28%5Cfrac%7B7%7D%7B3%7D%5Cright%29_%7B10%7D%20%3D%201.00%5Coverline%7B10%7D%2A2%5E1 "\left(\frac{7}{3}\right)_{10} = 1.00\overline{10}*2^1"),
    and after rounding,
    fl![\left(\frac{7}{3}\right) = 1.0010 . . . 1010 1011 \* 2^1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28%5Cfrac%7B7%7D%7B3%7D%5Cright%29%20%3D%201.0010%20.%20.%20.%201010%201011%20%2A%202%5E1 "\left(\frac{7}{3}\right) = 1.0010 . . . 1010 1011 * 2^1");
    and
    ![\left(\frac{4}{3}\right)\_{10} = 1.\overline{01} \* 2^0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28%5Cfrac%7B4%7D%7B3%7D%5Cright%29_%7B10%7D%20%3D%201.%5Coverline%7B01%7D%20%2A%202%5E0 "\left(\frac{4}{3}\right)_{10} = 1.\overline{01} * 2^0"),
    and after rounding,
    fl![\left(\frac{4}{3}\right) = 1.01 . . . 0101 0101 \* 2^0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28%5Cfrac%7B4%7D%7B3%7D%5Cright%29%20%3D%201.01%20.%20.%20.%200101%200101%20%2A%202%5E0 "\left(\frac{4}{3}\right) = 1.01 . . . 0101 0101 * 2^0").
    Subtracting the second from the first yields:

    ![\~\~\~\~1.0010\~1010\~1010\~1010\~1010\~1010\~1010\~1010\~1010\~1010\~1010\~1010\~1011\~0\*2^1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;~~~~1.0010~1010~1010~1010~1010~1010~1010~1010~1010~1010~1010~1010~1011~0%2A2%5E1 "~~~~1.0010~1010~1010~1010~1010~1010~1010~1010~1010~1010~1010~1010~1011~0*2^1")

    ![-0.1010\~1010\~1010\~1010\~1010\~1010\~1010\~1010\~1010\~1010\~1010\~1010\~1010\~1\*2^1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-0.1010~1010~1010~1010~1010~1010~1010~1010~1010~1010~1010~1010~1010~1%2A2%5E1 "-0.1010~1010~1010~1010~1010~1010~1010~1010~1010~1010~1010~1010~1010~1*2^1")

    ![-----------------------------------](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;----------------------------------- "-----------------------------------")

    ![=0.1000\~0000\~0000\~0000\~0000\~0000\~0000\~0000\~0000\~0000\~0000\~0000\~0000\~1\*2^1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%3D0.1000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~1%2A2%5E1 "=0.1000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~1*2^1")

    ![=1.0000\~0000\~0000\~0000\~0000\~0000\~0000\~0000\~0000\~0000\~0000\~0000\~0001\*2^0.\~\~](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%3D1.0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0001%2A2%5E0.~~ "=1.0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0001*2^0.~~")

    When we subtract 1 from this amount, we are left with
    ![\epsilon\_{mach}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7Bmach%7D "\epsilon_{mach}").

2.  Now
    ![\left(\frac{1}{3}\right)\_{10} = 1.\overline{01}\*2^{-2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28%5Cfrac%7B1%7D%7B3%7D%5Cright%29_%7B10%7D%20%3D%201.%5Coverline%7B01%7D%2A2%5E%7B-2%7D "\left(\frac{1}{3}\right)_{10} = 1.\overline{01}*2^{-2}"),
    and after rounding,
    fl![\left(\frac{1}{3}\right) = 1.01 . . . 0101 0101 \* 2^{−2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%28%5Cfrac%7B1%7D%7B3%7D%5Cright%29%20%3D%201.01%20.%20.%20.%200101%200101%20%2A%202%5E%7B%E2%88%922%7D "\left(\frac{1}{3}\right) = 1.01 . . . 0101 0101 * 2^{−2}").
    Subtracting this from
    ![\frac{4}{3}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cfrac%7B4%7D%7B3%7D "\frac{4}{3}"),
    we have:

    ![\~\~1.0101\~0101\~0101\~0101\~0101\~0101\~0101\~0101\~0101\~0101\~0101\~0101\~0101\~00\*2^0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;~~1.0101~0101~0101~0101~0101~0101~0101~0101~0101~0101~0101~0101~0101~00%2A2%5E0 "~~1.0101~0101~0101~0101~0101~0101~0101~0101~0101~0101~0101~0101~0101~00*2^0")

    ![-0.0101\~0101\~0101\~0101\~0101\~0101\~0101\~0101\~0101\~0101\~0101\~0101\~0101\~01\*2^0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;-0.0101~0101~0101~0101~0101~0101~0101~0101~0101~0101~0101~0101~0101~01%2A2%5E0 "-0.0101~0101~0101~0101~0101~0101~0101~0101~0101~0101~0101~0101~0101~01*2^0")

    ![-----------------------------------](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;----------------------------------- "-----------------------------------")

    ![=0.1111\~1111\~1111\~1111\~1111\~1111\~1111\~1111\~1111\~1111\~1111\~1111\~1111\~11\*2^0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%3D0.1111~1111~1111~1111~1111~1111~1111~1111~1111~1111~1111~1111~1111~11%2A2%5E0 "=0.1111~1111~1111~1111~1111~1111~1111~1111~1111~1111~1111~1111~1111~11*2^0")

    ![=1.1111\~1111\~1111\~1111\~1111\~1111\~1111\~1111\~1111\~1111\~1111\~1111\~1111\~1\*2^{-1},\~\~](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%3D1.1111~1111~1111~1111~1111~1111~1111~1111~1111~1111~1111~1111~1111~1%2A2%5E%7B-1%7D%2C~~ "=1.1111~1111~1111~1111~1111~1111~1111~1111~1111~1111~1111~1111~1111~1*2^{-1},~~")

    which rounds to
    ![10\*2^{-1}=1.0\*2^0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;10%2A2%5E%7B-1%7D%3D1.0%2A2%5E0 "10*2^{-1}=1.0*2^0").
    Subtracting 1 yields 0, not
    ![\epsilon\_{mach}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7Bmach%7D. "\epsilon_{mach}.")

Let’s test it out:

``` r
options(digits=20)
.Machine$double.eps
```

    ## [1] 2.2204460492503130808e-16

``` r
(7/3-4/3)-1
```

    ## [1] 2.2204460492503130808e-16

``` r
(4/3-1/3)-1
```

    ## [1] 0

### Problem 3

Using R, compute: 1000000000000000000 + 100 - 1000000000000000000. What
happens? Show some hand calculations to explain the result.

### Problem 3 Solution

``` r
1000000000000000000 + 100 - 1000000000000000000
```

    ## [1] 128

In binary, 1000000000000000000 is equal to (use an online converter)

![1101\~1110\~0000\~1011\~0110\~1011\~0011\~1010\~0111\~0110\~0100\~0000\~0000\~0000\~0000=1.10111100000101101101011001110100111011001000000000000000000\*2^{59}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1101~1110~0000~1011~0110~1011~0011~1010~0111~0110~0100~0000~0000~0000~0000%3D1.10111100000101101101011001110100111011001000000000000000000%2A2%5E%7B59%7D "1101~1110~0000~1011~0110~1011~0011~1010~0111~0110~0100~0000~0000~0000~0000=1.10111100000101101101011001110100111011001000000000000000000*2^{59}")

When we convert to double precision, we only keep 52 bits after the
decimal and the last 7 get rounded, but since they are all zeros, we can
exactly represent 1000000000000000000 as

![1.1011110000010110110101100111010011101100100000000000\*2^{59}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1.1011110000010110110101100111010011101100100000000000%2A2%5E%7B59%7D "1.1011110000010110110101100111010011101100100000000000*2^{59}")

Meanwhile, in binary, 1000000000000000100 is equal to

![1101\~1110\~0000\~1011\~0110\~1011\~0011\~1010\~0111\~0110\~0100\~0000\~0000\~0110\~0100=1.10111100000101101101011001110100111011001000000000001100100 \* 2^{59}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1101~1110~0000~1011~0110~1011~0011~1010~0111~0110~0100~0000~0000~0110~0100%3D1.10111100000101101101011001110100111011001000000000001100100%20%2A%202%5E%7B59%7D "1101~1110~0000~1011~0110~1011~0011~1010~0111~0110~0100~0000~0000~0110~0100=1.10111100000101101101011001110100111011001000000000001100100 * 2^{59}")

When we convert this to double precision, the last seven digits
(1100100) get chopped off, and the 0 in the 8th to last position gets
rounded up to 1, so we have

![1.1011110000010110110101100111010011101100100000000001 \* 2^{59}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1.1011110000010110110101100111010011101100100000000001%20%2A%202%5E%7B59%7D "1.1011110000010110110101100111010011101100100000000001 * 2^{59}")

Now the difference between 1000000000000000100 and 1000000000000000000
is

![0.000000000000000000000000000000000000000000000000001 \* 2^{59} = (10000000)\_2 = (128)\_{10}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.000000000000000000000000000000000000000000000000001%20%2A%202%5E%7B59%7D%20%3D%20%2810000000%29_2%20%3D%20%28128%29_%7B10%7D "0.000000000000000000000000000000000000000000000000001 * 2^{59} = (10000000)_2 = (128)_{10}")

This is another example of cancellation error, which occurs because the
function evaluation of
![x_1-x_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_1-x_2 "x_1-x_2")
is bad when
![x_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_1 "x_1")
and
![x_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_2 "x_2")
are relatively close.

### Problem 4

Find the smallest positive integer
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i")
such that
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i")
is not exactly representable using the IEEE standard in double
precision; i.e.,
![\hbox{fl}(i)\neq i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chbox%7Bfl%7D%28i%29%5Cneq%20i "\hbox{fl}(i)\neq i").

### Problem 4 Solution

Everything that has 52 bits in the mantissa can be represented exactly,
so we need 1.0\*2^53 at least, and then we need to put a 1 in the 53 bit
of the mantissa so that some rounding will occur:

![1.0000\~0000\~0000\~0000\~0000\~0000\~0000\~0000\~0000\~0000\~0000\~0000\~0000\~1 \* 2^{53},](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1.0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~1%20%2A%202%5E%7B53%7D%2C "1.0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~0000~1 * 2^{53},")

which is equal to 9007199254740993. This gets rounded down to
9007199254740992 (recall that the rounding rule when the next digit is 1
followed by all zeros is to round so that the new final digit is 0,
which in this case is rounding down). You can try this out!

``` r
print(9007199254740993,digits=16)
```

    ## [1] 9007199254740992

### Problem 5

Choose one section from 2.1 - 2.4 of [Bits and
Bugs](https://epubs.siam.org/doi/book/10.1137/1.9781611975567). Explain
the context of the problem you read about and give a mathematical
explanation of the error/challenge that occurred.

### Problem 5 Solution

Your solution goes here.
