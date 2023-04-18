------------------------------------------------------------------------

When you begin work during class, work with your assigned partner.
Please have only one electronic device open and work on it jointly. When
writing up this assignment, please remember that showing all of your
work and giving your reasoning are critical parts of achieving mastery.
If the course staff cannot tell how you solved a problem or finds leaps
in explanation or logic, the problem is not mastered. Finally, as a
matter of academic integrity, please make sure that you are positioned
to honestly answer yes to these questions:

-   Have I disclosed everyone with whom I collaborated on this work?
    (Even if it is only my assigned partner.)

-   Have I made a substantive intellectual contribution to the solution
    of every problem?

-   Am I making sure not to pass off as my own work any work that
    belongs to someone else?

Whether intentional or unintentional, any potential violations of
academic integrity will be referred to the Honor Committee.

------------------------------------------------------------------------

Load necessary packages:

    library(pracma)

### Problem 1

Let **A** be an *n* × *n* matrix. As a reminder, to perform power
iteration, we select any vector **x**<sup>(0)</sup> in ℝ<sup>*n*</sup>,
then compute
$\mathbf{x}^{(1)}=\frac{\mathbf{A}\mathbf{x}^{(0)}}{||\mathbf{A}\mathbf{x}^{(0)}||}$,
then compute
$\mathbf{x}^{(2)}=\frac{\mathbf{A}\mathbf{x}^{(1)}}{||\mathbf{A}\mathbf{x}^{(1)}||}$,
and so forth. From direct (repeated) substitution into the definition of
power iteration, one can find that the iterates satisfy
**x**<sup>(*k*)</sup> = *c*<sub>*k*</sub>**A**<sup>*k*</sup>**x**<sup>(0)</sup>,
where the sequence of constants {*c*<sub>*k*</sub>}, *k* = 0, 1, …, is
defined as

$$
c\_k=\prod\_{i=0}^{k-1} \frac{1}{ ||\mathbf{A}\mathbf{x}^{(i)}||}.
$$

(You can check this if you want.)

a\. Let **A****v** = *λ***v**. For *k* ∈ ℤ<sup>+</sup>, what is
**A**<sup>*k*</sup>**v**?

b\. An orthonormal basis for ℝ<sup>3</sup> is

$$
\mathbf{v}\_1 = \begin{pmatrix}-\frac{4}{9}, -\frac{7}{9}, \frac{4}{9} \end{pmatrix}^T,\quad
\mathbf{v}\_2 = \begin{pmatrix} \frac{1}{9},\frac{4}{9}, \frac{8}{9} \end{pmatrix}^T,\quad
\mathbf{v}\_3 = \begin{pmatrix} \frac{8}{9},-\frac{4}{9}, \frac{1}{9} \end{pmatrix}^T.
$$

Consider a starting vector
**x**<sup>(0)</sup> = (−1,−2,−3)<sup>*T*</sup>. Find *a*<sub>1</sub>,
*a*<sub>2</sub>, and *a*<sub>3</sub> such that
**x**<sup>(0)</sup> = *a*<sub>1</sub>**v**<sub>1</sub> + *a*<sub>2</sub>**v**<sub>2</sub> + *a*<sub>3</sub>**v**<sub>3</sub>.

c\. The orthonormal vectors **v**<sub>1</sub>, **v**<sub>2</sub>, and
**v**<sub>3</sub> are actually eigenvectors of

$$
\mathbf{A}=\frac{1}{9} \begin{pmatrix}
83 & 296 & -128 \\\\
296 & 473 & -152 \\\\
-128 & -152 & 335
\end{pmatrix},
$$

and they are ordered so that their associated eigenvalues are in order
of decreasing magnitude:
|*λ*<sub>1</sub>| &gt; |*λ*<sub>2</sub>| ≥ |*λ*<sub>3</sub>|. Let’s use
this and put together with what we know from having worked through this
problem this far:

$$
\begin{align\*}
\mathbf{x}^{(k)}&=c\_k \mathbf{A}^k \mathbf{x}^{(0)} \\\\
&=c\_k \mathbf{A}^k (a\_1 \mathbf{v}\_1 + a\_2 \mathbf{v}\_2 + a\_3 \mathbf{v}\_3) \\\\
&=c\_k (a\_1 \mathbf{A}^k \mathbf{v}\_1 + a\_2 \mathbf{A}^k \mathbf{v}\_2 + a\_3 \mathbf{A}^k \mathbf{v}\_3)\\\\
&=c\_k (a\_1 \lambda\_1^k \mathbf{v}\_1 + a\_2 \lambda\_2^k \mathbf{v}\_2 + a\_3 \lambda\_3^k \mathbf{v}\_3)\\\\
&=c\_k \lambda\_1^k \left\[a\_1 \mathbf{v}\_1 + a\_2 \left(\frac{\lambda\_2}{\lambda\_1}\right)^k \mathbf{v}\_2 + a\_3 \left(\frac{\lambda\_3}{\lambda\_1}\right)^k \mathbf{v}\_3\right\]
\end{align\*}
$$

Explain what happens as *k* increases, meaning, as we continue
performing power iteration Note: we need *a*<sub>1</sub> ≠ 0. Give a
geometric explanation of when that condition is met.

d. Now write code to perform power iteration using the
**x**<sup>(0)</sup> and **A** defined above. In each iteration of your
loop, also compute and print out

$$
\lambda^{(k)}=\frac{\left(\mathbf{x}^{(k)}\right)^T \mathbf{A} \mathbf{x}^{(k)}}{\left(\mathbf{x}^{(k)}\right)^T \mathbf{x}^{(k)}}.
$$

This is called the Rayleigh quotient, and it is the estimate for the
eigenvalue. What vector does your power iteration converge to? What
value does your Rayleigh quotient converge to?

### Problem 1 Solution

a\. Your solution goes here.

b\. Your solution goes here.

c\. Your solution goes here.

d\. Your solution goes here.
