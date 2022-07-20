### Problem 1

Convert the following base-10 numbers to binary.

a\. 11.25

b\. 2/3

c\. 99.9

### Problem 1 Solution

a\. 11.25 = 11 + 0.25. Integer part:
(11)<sub>10</sub> = (8+2+1)<sub>10</sub> = (1011)<sub>2</sub>.
Fractional part: (0.25)<sub>10</sub> = (0.01)<sub>2</sub>, so
(11.25)<sub>10</sub> = (1011.01)<sub>2</sub>.

b\. We adopt the strategy of successively multiplying by two and pulling
off a bit.

$$
\frac{2}{3}\*2 = \frac{1}{3} + 1
$$

$$
\frac{1}{3}\*2 = \frac{2}{3} + 0
$$

$$
\frac{2}{3}\*2 = \frac{1}{3} + 1
$$

⋮

So $(\frac{2}{3})\_{10}=(0.\overline{10})\_2$

c\. (99.9)<sub>10</sub> = (99)<sub>10</sub> + (0.9)<sub>10</sub>. For
the first part,
(99)<sub>10</sub> = (64+32+2+1)<sub>10</sub> = (1100011)<sub>2</sub>.
For the fractional part, we adopt the strategy of successively
multiplying by two and pulling off a bit.

0.9 \* 2 = 0.8 + 1

0.8 \* 2 = 0.6 + 1

0.6 \* 2 = 0.2 + 1

0.2 \* 2 = 0.4 + 0

0.4 \* 2 = 0.8 + 0

0.8 \* 2 = 0.6 + 1

⋮

So $(0.9)\_{10}=(0.1\overline{1100})\_2$, and therefore
$(99.9)\_{10}=(1100011.1\overline{1100})\_{2}$

### Problem 2

a\. Explain why you can determine machine epsilon on a computer using
IEEE double precision and the IEEE Rounding to Nearest Rule by
calculating (7/3 − 4/3) − 1. I’m not asking for (just) a conceptual
explanation; you should actually calculate some things on paper to make
your argument.

b\. Does (4/3 − 1/3) − 1 also give *ϵ*<sub>*m**a**c**h*</sub>? Explain
by converting to floating point numbers and carrying out the machine
arithmetic.

### Problem 2 Solution

a\. Let’s convert the first two terms:
$\left(\frac{7}{3}\right)\_{10} = 1.00\overline{10}\*2^1$, and after
rounding, fl$\left(\frac{7}{3}\right) = 1.0010 . . . 1010 1011 \* 2^1$;
and $\left(\frac{4}{3}\right)\_{10} = 1.\overline{01} \* 2^0$, and after
rounding, fl$\left(\frac{4}{3}\right) = 1.01 . . . 0101 0101 \* 2^0$.
Subtracting the second from the first yields:

    1.0010 1010 1010 1010 1010 1010 1010 1010 1010 1010 1010 1010 1011 0 \* 2<sup>1</sup>

 − 0.1010 1010 1010 1010 1010 1010 1010 1010 1010 1010 1010 1010 1010 1 \* 2<sup>1</sup>

<table>
<thead>
<tr class="header">
<th style="text-align: left;">$$</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">$$</td>
</tr>
</tbody>
</table>

 = 0.1000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 1 \* 2<sup>1</sup>

 = 1.0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0001 \* 2<sup>0</sup>.  

When we subtract 1 from this amount, we are left with
*ϵ*<sub>*m**a**c**h*</sub>.

b\. Now $\left(\frac{1}{3}\right)\_{10} = 1.\overline{01}\*2^{-2}$, and
after rounding,
fl$\left(\frac{1}{3}\right) = 1.01 . . . 0101 0101 \* 2^{−2}$.
Subtracting this from $\frac{4}{3}$, we have:

  1.0101 0101 0101 0101 0101 0101 0101 0101 0101 0101 0101 0101 0101 00 \* 2<sup>0</sup>

 − 0.0101 0101 0101 0101 0101 0101 0101 0101 0101 0101 0101 0101 0101 01 \* 2<sup>0</sup>

<table>
<thead>
<tr class="header">
<th style="text-align: left;">$$</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">$$</td>
</tr>
</tbody>
</table>

 = 0.1111 1111 1111 1111 1111 1111 1111 1111 1111 1111 1111 1111 1111 11 \* 2<sup>0</sup>

 = 1.1111 1111 1111 1111 1111 1111 1111 1111 1111 1111 1111 1111 1111 1 \* 2<sup>−1</sup>,  

which rounds to 10 \* 2<sup>−1</sup> = 1.0 \* 2<sup>0</sup>. Subtracting
1 yields 0, not *ϵ*<sub>*m**a**c**h*</sub>.

Let’s test it out:

    options(digits=20)
    .Machine$double.eps

    ## [1] 2.2204460492503130808e-16

    (7/3-4/3)-1

    ## [1] 2.2204460492503130808e-16

    (4/3-1/3)-1

    ## [1] 0

### Problem 3

Using R, compute: 1000000000000000000 + 100 - 1000000000000000000. What
happens? Show some hand calculations to explain the result.

### Problem 3 Solution

    1000000000000000000 + 100 - 1000000000000000000

    ## [1] 128

In binary, 1000000000000000000 is equal to (use an online converter)

1101 1110 0000 1011 0110 1011 0011 1010 0111 0110 0100 0000 0000 0000 0000 = 1.10111100000101101101011001110100111011001000000000000000000 \* 2<sup>59</sup>

When we convert to double precision, we only keep 52 bits after the
decimal and the last 7 get rounded, but since they are all zeros, we can
exactly represent 1000000000000000000 as

1.1011110000010110110101100111010011101100100000000000 \* 2<sup>59</sup>

Meanwhile, in binary, 1000000000000000100 is equal to

1101 1110 0000 1011 0110 1011 0011 1010 0111 0110 0100 0000 0000 0110 0100 = 1.10111100000101101101011001110100111011001000000000001100100 \* 2<sup>59</sup>

When we convert this to double precision, the last seven digits
(1100100) get chopped off, and the 0 in the 8th to last position gets
rounded up to 1, so we have

1.1011110000010110110101100111010011101100100000000001 \* 2<sup>59</sup>

Now the difference between 1000000000000000100 and 1000000000000000000
is

0.000000000000000000000000000000000000000000000000001 \* 2<sup>59</sup> = (10000000)<sub>2</sub> = (128)<sub>10</sub>

This is another example of cancellation error, which occurs because the
function evaluation of *x*<sub>1</sub> − *x*<sub>2</sub> is bad when
*x*<sub>1</sub> and *x*<sub>2</sub> are relatively close.

### Problem 4

Find the smallest positive integer *i* such that *i* is not exactly
representable using the IEEE standard in double precision; i.e.,
$\hbox{fl}(i)\neq i$.

### Problem 4 Solution

Everything that has 52 bits in the mantissa can be represented exactly,
so we need 1.0\*2^53 at least, and then we need to put a 1 in the 53 bit
of the mantissa so that some rounding will occur:

1.0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 1 \* 2<sup>53</sup>,

which is equal to 9007199254740993. This gets rounded down to
9007199254740992 (recall that the rounding rule when the next digit is 1
followed by all zeros is to round so that the new final digit is 0,
which in this case is rounding down). You can try this out!

    print(9007199254740993,digits=16)

    ## [1] 9007199254740992

### Problem 5

Choose one section from 2.1 - 2.4 of [Bits and
Bugs](https://epubs.siam.org/doi/book/10.1137/1.9781611975567). Explain
the context of the problem you read about and give a mathematical
explanation of the error/challenge that occurred.

### Problem 5 Solution

Your solution goes here.
