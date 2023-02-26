Recap Exercise \#1
================

------------------------------------------------------------------------

This recap exercise provides a formative assessment of your
understanding of the course material. No collaboration is allowed.
Please work on your own. Once you complete your first attempt and
receive written feedback, you can seek support from the TAs and
professor if you want, but please continue not to collaborate with
peers. Any/all work you turn in should be your own work. Whether
intentional or unintentional, any potential violations of academic
integrity will be referred to the Honor Committee.

------------------------------------------------------------------------

Load necessary packages:

``` r
library(Matrix)
knitr::opts_chunk$set(options(digits=8))
```

### Problem 1

Assume that for all parts of this problem, I am using the same computer
and software. I have a $1000 \times 1000$ dense matrix $\mathbf{A}$
without special structure. My computer solves
$\mathbf{A} \mathbf{x}=\mathbf{b}$ using Gaussian elimination and back
substitution in 10 seconds.

a\. I’d like to solve another problem
$\mathbf{C} \mathbf{y}=\mathbf{g}$, where $\mathbf{C}$ is a dense
$10^5 \times 10^5$ matrix without special structure. How many hours will
it take to solve $\mathbf{C} \mathbf{y}=\mathbf{g}$ with the same method
as before?

b\. I have a third matrix $\mathbf{D}$ that is $1000 \times 1000$,
dense, and has no special structure. How many hours should it take to
solve

$$
\mathbf{D} \mathbf{x}_1 = \mathbf{f}_1, \ldots, \mathbf{D} \mathbf{x}_{100} = \mathbf{f}_{100}
$$

if I am being as efficient as possible?

### Problem 1 Solution

a\. Your solution goes here.

b\. Your solution goes here.

### Problem 2

\[Warning: Once, someone thought this problem was real, so let me assure
you that I am completely making this up.\] Before Apple manufactured the
iPhone, there was the hPhone, which, bizarrely, was a 24-bit machine: 1
bit for the sign, 16 for the mantissa, and 7 for the exponent. The
floating point representation corrected for exponential bias in a manner
analogous to the IEEE 64-bit standard. The largest exponent on the
hPhone was reserved to represent infinity.

a\. Which number is closest to the largest exactly representable number
on the hPhone? (Make sure to show/explain your work.)

A. $10^6$  
B. $10^{10}$  
C. $10^{19}$  
D. $10^{31}$  
E. $10^{63}$  
F. $10^{308}$

b\. On the hPhone, what number is the base-10 number 5901.546875
actually stored as? Give your answer in base-10.

### Problem 2 Solution

a\. Your solution goes here.

b\. Your solution goes here.

### Problem 3

The table below shows $||\mathbf{x}^i||$ and
$||\mathbf{A} \mathbf{x}^i||_{\infty}$ for four different vectors
$\mathbf{x}^i$, $i=1,\ldots,4$. $\mathbf{A}$ is a fixed but unknown
$n \times n$ invertible matrix. What are the best (largest) lower bounds
on $||\mathbf{A}||_\infty$, $||\mathbf{A}^{-1}||_\infty$, and
$\kappa_\infty(\mathbf{A})$?

| $i$ | $\|\|\mathbf{x}^i\|\|_{\infty}$ | $\|\|\mathbf{A} \mathbf{x}^i\|\|_{\infty}$ |
|:---:|:-------------------------------:|:------------------------------------------:|
| $1$ |               $1$               |                   $10^2$                   |
| $2$ |             $10^3$              |                   $10^4$                   |
| $3$ |             $10^2$              |                    $1$                     |
| $4$ |            $10^{-3}$            |                   $10^2$                   |

### Problem 3 Solution

Your solution goes here.

### Problem 4

Let’s learn about numerical calculation of derivatives and the
associated error. Write a function called `D2` that takes a function
$f(x)$, a value $x_0$ and a value of a quantity we’ll call $h$ as inputs
and approximates the second derivative $f^{\prime \prime}(x_0)$ as

$$
f^{\prime \prime}(x_0) \approx \frac{f(x_0-h)-2f(x_0)+f(x_0+h)}{h^2}
$$

Using this function, make a graph of the absolute error for $f(x)=e^x$
at $x_0=0$ as a function of h, and find the optimal size of $h$. You’ll
want to check many orders of magnitude for $h$, and so equally spaced
values are not appropriate. Instead, use values whose logs are equally
spaced by taking

``` r
h <- 10^seq(from = -8, to = -1, length = 501)
```

Then, in a few sentences, explain why it is not optimal to take $h$ as
small as possible for the numerical computation you carried out above.
Be specific.

### Problem 4 Solution

Your solution goes here.

### Problem 5

There are many different ways to rank and rate teams and individuals in
competitions. The following is the method behind one of the computer
rankings that used to be used to determine the college football playoff
teams.

For each team $i$, the algorithm uses the following information:

- $w_i=$ the number of wins for team $i$
- $l_i=$ the number of losses for team $i$
- $t_i=w_i+l_i=$ the total number of games played by team $i$

The goal is to determine $\mathbf{r}$, a vector of ratings, with the
i<sup>th</sup> element $r_i$ representing the rating for team $i$.
Rankings are just the order of the ratings (the highest $r_i$ is given
the \#1 ranking spot, etc.).

The ratings are computed as follows:

$$
r_i = \frac{1}{2+t_i}\left(1+\frac{w_i-l_i}{2}+\sum\limits_{j \in \{i\hbox{'s opponents}\}}r_j\right)
$$

where “opponent” means a team that a given team played (*not* all the
other teams in the league).

Try out this rating system on some of the local college football teams
from my old neighborhood in Minnesota, based on actual results from the
schedule some years back. Here are the results of 12 games from 6 local
teams: Macalester (M), Carleton (C), Hamline (H) , St. Thomas (T),
Gustavus Adolphus (G), and St. Olaf (O):

$$
\begin{array}{lllllll}
\text{M beat C}&\hskip.2in&\text{T beat C}&\hskip.2in&\text{T beat H} &\hskip.2in& \text{G beat H} \\
\text{H beat M}&&\text{G beat C}&&\text{T beat G} && \text{H beat O} \\
\text{C beat H}&&\text{O beat C}&&\text{T beat O} && \text{G beat O} \\
\end{array}
$$

By solving a linear system $\mathbf{A}\,\mathbf{r}=\mathbf{b}$, find the
ratings vector $\mathbf{r}$ and use it to give the ranking of teams from
first place to last place. In answering this problem, let
$i=(1,2,3,4,5,6)$ correspond to $(M,C,H,T,G,O)$ respectively. If you
find ties between teams you may choose their relative ordering.

### Problem 5 solution

Your solution goes here.

### Problem 6

a\. Choose a data set describing a large network (at least hundreds of
nodes) and represent it as an adjacency matrix in R. Explain the context
of the problem (that is, tell me what the nodes and edges in the network
mean). Tell me how many nodes and edges there are. Finally, tell me why
you chose this particular network — what about it seems important or
interesting? Hint: if the task of loading up a network seems daunting to
you, you can always adapt the code I provided on the Solving Linear
Systems Pset by providing a link to any other network in the Stanford
Large Network Data Set Collection. Of course, you can certainly branch
out to other sources if you prefer!

b\. Calculate the most central node according to Katz centrality.

c\. Visualize the network and make the most central node appear in red.

d\. Discuss any reasons that knowing the most central node is useful or
important (this is a question about the applied context of the problem,
not about the mathematics of it).

### Problem 6 Solution

a\. Your solution goes here.

b\. Your solution goes here.

c\. Your solution goes here.

d\. Your solution goes here.
