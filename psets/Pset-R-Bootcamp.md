Pset - R Bootcamp
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

Load necessary packages.

``` r
library(pracma)
```

### Problem 1

In the 1999 movie Office Space, a character creates a program that takes
fractions of cents that are truncated in a bank’s transactions and
deposits them into his own account. This is not a new idea, and hackers
who have actually attempted it have been arrested. In this exercise you
will explore this scam. The point of this problem is to practice some R
skills.

Assume the following details:

- You have access to 10,000 bank accounts.

- Initially, the bank accounts have balances that are uniformly
  distributed between \$100 and \$100,000.

- The annual interest rate on the accounts is 5%.

- Interest is compounded each day and added to the accounts, except that
  fractions of a cent are truncated.

- The daily interest rate is thus $.05/365$.

- The truncated fractions are deposited into an illegal account that
  initially has a balance of \$0.

- The illegal account can hold fractional values and it also accrues
  daily interest.

- Assume that each day, interest is earned (on the main account and the
  illegal account) first, and then skimming from the main account into
  the illegal account happens.

Your job is to write an R script that simulates this situation and finds
how long it takes for the illegal account to reach a million dollars.

Here is some R help.

The following code generates the initial accounts:

``` r
accounts <- runif(10000,100,100000)
accounts <- floor(accounts*100)/100
```

The first line sets up 10,000 accounts with values uniformly between 100
and 100,000. The second line removes the fractions of cents (look at the
data before and after that line is applied).

For the main part of your code, I strongly recommend a `while` loop,
which would look something like:

``` r
while (illegal < 1000000) {
  do stuff here
}
```

In your response, before you provide your code and the final answer,
write out your algorithm in words. I want to read something like “Here’s
the procedure we’ll follow. First, we initialize 10,000 accounts each
with a random amount of money between \$100 and \$100,000, making sure
to round them so that there are no fractional cents. Next, we initialize
the illegal account to have zero dollars, and we set a day counter to
begin at 0. Then, we perform the following procedure until the illegal
account meets or exceeds the \$1,000,000 threshold:

- Apply interest to the illegal account
- Apply interest to the main accounts
- Etc. etc.

### Problem 1 Solution

Your solution goes here.

### Problem 2

The purpose of this problem is to practice some R skills and review a
critical topic in applied mathematics: Taylor series.

a\. Find the degree 5 Taylor polynomial $P(x)$ centered at $x = 0$ for
$f(x) = \cos x$. This is a paper-and-pencil calculation, but you should
write up the solution here in your markdown document.

b\. Find an upper bound for the error in approximating $f(x)$ by $P(x)$
on the interval $[-\pi/4,\pi/4]$. This is also a paper-and-pencil
calculation that you should write up here.

c\. Make a single plot with 3 different curves: $cos(x)$ and the 2nd and
4th degree Taylor approximations of $f(x)$ around $x=0$. Plot these on
the interval $[-\pi/2,\pi/2]$.

### Problem 2 Solution

a\. Your solution goes here.

b\. Your solution goes here.

c\. Your solution goes here.

### Problem 3

Consider the polynomial $f(x) = x^5 + x^4 - 2x^3 + 4x^2 - 8x + 5$. For a
vector of 5 million equally spaced values going from x = -100 to x =
100, evaluate the polynomial four different ways, time each one, and
compare. You should make your vector of x values using the `seq`
command.

- Regular polynomial evaluation using a loop to iterate over each value
  in the vector.
- Regular polynomial evaluation using the entire vector at once.
- Horner’s method using a loop to iterate over each value in the vector.
- Horner’s method using the entire vector at once.

In addition to stating the times for each calculation, identify the
slowest, the fastest, and the ratio of slowest to fastest.

### Problem 3 Solution

Your solution goes here.

### Problem 4

You read three articles on applications of linear algebra:

- [The Mathematics of Mass Testing for
  COVID-19](https://sinews.siam.org/Details-Page/the-mathematics-of-mass-testing-for-covid-19)
- [Machine Learning Has Been Used to Automatically Translate Long-Lost
  Languages](https://www.technologyreview.com/2019/07/01/65601/machine-learning-has-been-used-to-automatically-translate-long-lost-languages/)
- [The Best
  Bits](https://www.americanscientist.org/article/the-best-bits)

Choose one of these articles and write a brief paragraph about it in
which you explain the role that linear algebra plays. Strive to be
specific. For instance, rather than writing “linear algebra is used to
make pooled testing more efficient,” try to explain the actual
mathematical ideas used.

### Problem 4

Your solution goes here.
