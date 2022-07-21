    library(Matrix)
    knitr::opts_chunk$set(options(digits=8))

### Problem 1

Assume that for all parts of this problem, I am using the same computer
and software. I have a 1000 × 1000 dense matrix **A** without special
structure. My computer solves **A** **x** = **b** using Gaussian
elimination and back substitution in 10 seconds. Show your work before
the code block and please give final answers in the code block.

a\. I’d like to solve another problem **C** **y** = **g**, where **C**
is a dense 10<sup>5</sup> × 10<sup>5</sup> matrix without special
structure. How many hours will it take to solve **C** **y** = **g** with
the same method as before?

b\. I have a third matrix **D** that is 1000 × 1000, dense, and has no
special structure. How many hours should it take to solve
**D** **x**<sub>1</sub> = **f**<sub>1</sub>, …, **D** **x**<sub>100</sub> = **f**<sub>100</sub>
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

A. 10<sup>6</sup>  
B. 10<sup>10</sup>  
C. 10<sup>19</sup>  
D. 10<sup>31</sup>  
E. 10<sup>63</sup>  
F. 10<sup>308</sup>

b\. On the hPhone, what number is the base-10 number 5901.546875
actually stored as? Give your answer in base-10.

### Problem 2 Solution

a\. Your solution goes here.

b\. Your solution goes here.

### Problem 3

The table below shows ||**x**<sup>*i*</sup>||<sub>∞</sub> and
||**A** **x**<sup>*i*</sup>||<sub>∞</sub> for four different vectors
**x**<sup>*i*</sup>, *i* = 1, …, 4. **A** is a fixed but unknown
*n* × *n* invertible matrix. What are the best (largest) lower bounds on
||**A**||<sub>∞</sub>, ||**A**<sup>−1</sup>||<sub>∞</sub>, and
*κ*<sub>∞</sub>(**A**)?

<table>
<colgroup>
<col style="width: 10%" />
<col style="width: 40%" />
<col style="width: 50%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: center;"><span
class="math inline"><em>i</em></span></th>
<th style="text-align: center;"><span
class="math inline">||<strong>x</strong><sup><em>i</em></sup>||<sub>∞</sub></span></th>
<th style="text-align: center;"><span
class="math inline">||<strong>A</strong><strong>x</strong><sup><em>i</em></sup>||<sub>∞</sub></span></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: center;"><span class="math inline">1</span></td>
<td style="text-align: center;"><span class="math inline">1</span></td>
<td style="text-align: center;"><span
class="math inline">10<sup>2</sup></span></td>
</tr>
<tr class="even">
<td style="text-align: center;"><span class="math inline">2</span></td>
<td style="text-align: center;"><span
class="math inline">10<sup>3</sup></span></td>
<td style="text-align: center;"><span
class="math inline">10<sup>4</sup></span></td>
</tr>
<tr class="odd">
<td style="text-align: center;"><span class="math inline">3</span></td>
<td style="text-align: center;"><span
class="math inline">10<sup>2</sup></span></td>
<td style="text-align: center;"><span class="math inline">1</span></td>
</tr>
<tr class="even">
<td style="text-align: center;"><span class="math inline">4</span></td>
<td style="text-align: center;"><span
class="math inline">10<sup>−3</sup></span></td>
<td style="text-align: center;"><span
class="math inline">10<sup>2</sup></span></td>
</tr>
</tbody>
</table>

### Problem 3 Solution

To think about ||**A**||<sub>∞</sub>, we think about the maximum
expansion achieved on a vector when multiplied by **A**, which here is
10<sup>5</sup> for *i* = 4. Similarly, for
||**A**<sup>−</sup>1||<sub>∞</sub>, the maximum expansion is
10<sup>2</sup> for *i* = 3. Therefore, the condition number is at least
10<sup>5</sup> ⋅ 10<sup>2</sup> = 10<sup>7</sup>.

### Problem 4

Let’s learn about numerical calculation of derivatives and the
associated error. Write a function called `D2` that takes a function
*f*(*x*), a value *x*<sub>0</sub> and a value of a quantity we’ll call
*h* as inputs and approximates the second derivative
*f*<sup>′′</sup>(*x*<sub>0</sub>) as

$$
f^{\prime \prime}(x\_0) \approx \frac{f(x\_0-h)-2f(x\_0)+f(x\_0+h)}{h^2}
$$

Using this function, make a graph of the absolute error for
*f*(*x*) = *e*<sup>*x*</sup> at *x*<sub>0</sub> = 0 and find the optimal
size of *h*. Then, in a few sentences, explain why it is not optimal to
take *h* as small as possible for the numerical computation you carried
out above. Be specific.

### Problem 4 Solution

Your solution goes here.

### Problem 5

There are many different ways to rank and rate teams and individuals in
competitions. The following is the method behind one of the computer
rankings that used to be used to determine the college football playoff
teams.

For each team *i*, the algorithm uses the following information:

-   *w*<sub>*i*</sub>= the number of wins for team *i*
-   *l*<sub>*i*</sub>= the number of losses for team *i*
-   *t*<sub>*i*</sub> = *w*<sub>*i*</sub> + *l*<sub>*i*</sub>= the total
    number of games played by team *i*

The goal is to determine **r**, a vector of ratings, with the
*i*<sup>*t**h*</sup> element *r*<sub>*i*</sub> representing the rating
for team *i*. Rankings are just the order of the ratings (the highest
*r*<sub>*i*</sub> is given the \#1 ranking spot, etc.).

The ratings are computed as follows:

$$
r\_i = \frac{1}{2+t\_i}\left(1+\frac{w\_i-l\_i}{2}+\sum\limits\_{j \in \\{i\hbox{'s opponents}\\}}r\_j\right)
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
\text{M beat C}&\hskip.2in&\text{T beat C}&\hskip.2in&\text{T beat H} &\hskip.2in& \text{G beat H} \\\\
\text{H beat M}&&\text{G beat C}&&\text{T beat G} && \text{H beat O} \\\\
\text{C beat H}&&\text{O beat C}&&\text{T beat O} && \text{G beat O} \\\\
\end{array}
$$

By solving a linear system **A** **r** = **b**, find the ratings vector
**r** and use it to give the ranking of teams from first place to last
place. In answering this problem, let *i* = (1,2,3,4,5,6) correspond to
(*M*,*C*,*H*,*T*,*G*,*O*) respectively. If you find ties between teams
you may choose their relative ordering.

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

1.  Your solution goes here.
