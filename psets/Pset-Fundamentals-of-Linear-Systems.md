Pset - How Computers Store Numbers
================

Load necessary packages:

``` r
library("Matrix")
library("igraph")
library("pracma")
```

### Problem 1

A real
![n \\times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%20n "n \times n")
matrix
![\\mathbf{Q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BQ%7D "\mathbf{Q}")
is *orthonormal* if

![ \\mathbf{Q}^{T}\\mathbf{Q}=\\mathbf{Q}\\mathbf{Q}^T=\\mathbf{I};\~i.e., \\mathbf{Q}^{-1}=\\mathbf{Q}^T.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%20%5Cmathbf%7BQ%7D%5E%7BT%7D%5Cmathbf%7BQ%7D%3D%5Cmathbf%7BQ%7D%5Cmathbf%7BQ%7D%5ET%3D%5Cmathbf%7BI%7D%3B~i.e.%2C%20%5Cmathbf%7BQ%7D%5E%7B-1%7D%3D%5Cmathbf%7BQ%7D%5ET. " \mathbf{Q}^{T}\mathbf{Q}=\mathbf{Q}\mathbf{Q}^T=\mathbf{I};~i.e., \mathbf{Q}^{-1}=\mathbf{Q}^T.")

Note: sometimes you will see these matrices just called *orthogonal
matrices*.

Show that
![\|\|\\mathbf{Q}\|\|\_2=\|\|\\mathbf{Q}^{-1}\|\|\_2=1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7C%7C%5Cmathbf%7BQ%7D%7C%7C_2%3D%7C%7C%5Cmathbf%7BQ%7D%5E%7B-1%7D%7C%7C_2%3D1 "||\mathbf{Q}||_2=||\mathbf{Q}^{-1}||_2=1"),
and therefore the 2-norm condition number of any orthonormal matrix
![\\mathbf{Q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BQ%7D "\mathbf{Q}")
is
![\\kappa\_2(\\mathbf{Q})=\|\|\\mathbf{Q}\|\|\_2 \|\|\\mathbf{Q}^{-1}\|\|\_2=1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ckappa_2%28%5Cmathbf%7BQ%7D%29%3D%7C%7C%5Cmathbf%7BQ%7D%7C%7C_2%20%7C%7C%5Cmathbf%7BQ%7D%5E%7B-1%7D%7C%7C_2%3D1 "\kappa_2(\mathbf{Q})=||\mathbf{Q}||_2 ||\mathbf{Q}^{-1}||_2=1").

### Problem 1 Solution

Your solution goes here.

### Problem 2

Consider the
![n \\times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%20n "n \times n")
square matrix

![A=\\begin{bmatrix} 1 & -1 & -1 & -1 & \\ldots & -1 \\\\ 0 & 1 & -1 & -1 & \\ldots & -1 \\\\ 0 & 0 & 1 & -1 & \\ldots & -1 \\\\ \\vdots & \\vdots & \\vdots & \\vdots & \\vdots & \\vdots \\\\ 0 & 0 & 0 & 0 & 0 & 1 \\end{bmatrix}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A%3D%5Cbegin%7Bbmatrix%7D%201%20%26%20-1%20%26%20-1%20%26%20-1%20%26%20%5Cldots%20%26%20-1%20%5C%5C%200%20%26%201%20%26%20-1%20%26%20-1%20%26%20%5Cldots%20%26%20-1%20%5C%5C%200%20%26%200%20%26%201%20%26%20-1%20%26%20%5Cldots%20%26%20-1%20%5C%5C%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%5C%5C%200%20%26%200%20%26%200%20%26%200%20%26%200%20%26%201%20%5Cend%7Bbmatrix%7D. "A=\begin{bmatrix} 1 & -1 & -1 & -1 & \ldots & -1 \\ 0 & 1 & -1 & -1 & \ldots & -1 \\ 0 & 0 & 1 & -1 & \ldots & -1 \\ \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\ 0 & 0 & 0 & 0 & 0 & 1 \end{bmatrix}.")

1.  Write a function that takes a value of
    ![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
    as input and outputs the matrix A above. Challenge yourself to be
    efficient in your coding when you create the matrix. Looking at
    patterns and Googling how to create them helps. I managed to create
    the matrix with no loops, in two lines of code.

2.  Use your function and the R command `kappa` to calculate the
    approximate condition number
    ![\\kappa(\\mathbf{A})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ckappa%28%5Cmathbf%7BA%7D%29 "\kappa(\mathbf{A})")
    for
    ![n=1,\\ldots,30](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%3D1%2C%5Cldots%2C30 "n=1,\ldots,30").
    Plot
    ![\\log\_{10}\[\\kappa(\\mathbf{A})\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clog_%7B10%7D%5B%5Ckappa%28%5Cmathbf%7BA%7D%29%5D "\log_{10}[\kappa(\mathbf{A})]")
    as a function of
    ![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
    and use the commands `lm` and `abline` to plot a best fit line. How
    would you describe the conditioning of the matrix as
    ![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
    increases?

3.  Now choose
    ![n=30](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%3D30 "n=30").
    Generate the matrix
    ![\\mathbf{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D "\mathbf{A}")
    and let
    ![\\mathbf{b\_1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bb_1%7D "\mathbf{b_1}")
    be an
    ![n \\times 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%201 "n \times 1")
    vector of random numbers chosen uniformly from 0 to 1. Solve
    ![\\mathbf{A} \\mathbf{x}\_1 = \\mathbf{b}\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D%20%5Cmathbf%7Bx%7D_1%20%3D%20%5Cmathbf%7Bb%7D_1 "\mathbf{A} \mathbf{x}_1 = \mathbf{b}_1")
    using any appropriate method that you want (including R’s built-in
    capabilities). Now let
    ![\\mathbf{b}\_2 = \\mathbf{b}\_1 + (0,\\ldots,0,0.001)^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bb%7D_2%20%3D%20%5Cmathbf%7Bb%7D_1%20%2B%20%280%2C%5Cldots%2C0%2C0.001%29%5ET "\mathbf{b}_2 = \mathbf{b}_1 + (0,\ldots,0,0.001)^T"),
    that is, you leave the first 29 elements the same as in
    ![\\mathbf{b}\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bb%7D_1 "\mathbf{b}_1")
    but add 0.001 to the last element. Solve
    ![\\mathbf{A} \\mathbf{x}\_2 = \\mathbf{b}\_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BA%7D%20%5Cmathbf%7Bx%7D_2%20%3D%20%5Cmathbf%7Bb%7D_2 "\mathbf{A} \mathbf{x}_2 = \mathbf{b}_2").
    Use the command `Norm` to find the (approximate)2-norm of
    ![\\mathbf{x\_1}-\\mathbf{x\_2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bx_1%7D-%5Cmathbf%7Bx_2%7D "\mathbf{x_1}-\mathbf{x_2}")
    and discuss vis-a-vis your result from part b. Also, to help build
    your intuition, find the magnitude of the difference between the
    first coordinate of
    ![\\mathbf{x}\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bx%7D_1 "\mathbf{x}_1")
    and that of
    ![\\mathbf{x}\_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bx%7D_2 "\mathbf{x}_2").

### Problem 2 Solution

Your solution goes here.

### Problem 3

Note: This problem is taken from Linear Algebra and its Applications, by
Lay et al.

<img src="https://github.com/chadtopaz/computationallineaaralgebra/raw/main/psets/heat.png" width="400">

An important concern in the study of heat transfer is to determine the
steady-state temperature distribution of a thin plate when the
temperature around the boundary is known. Assume the plate shown in the
figure above represents a cross section of a metal beam, with negligible
heat flow in the direction perpendicular to the plate. Let the variables
![x\_1, x\_2, \\ldots, x\_8](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_1%2C%20x_2%2C%20%5Cldots%2C%20x_8 "x_1, x_2, \ldots, x_8")
denote the temperatures at nodes 1 through 8 in the picture. In steady
state, the temperature at a node is approximately equal to the average
of the four nearest nodes (to the left, above, right, below).

1.  The solution to the approximate steady-state heat flow problem for
    this plate can be written as a system of linear equations
    ![Ax=b](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Ax%3Db "Ax=b"),
    where
    ![x=\[x\_1, x\_2, \\ldots, x\_8\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%3D%5Bx_1%2C%20x_2%2C%20%5Cldots%2C%20x_8%5D "x=[x_1, x_2, \ldots, x_8]")
    is the vector of temperatures at nodes 1 through 8. Find the
    ![8 \\times 8](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;8%20%5Ctimes%208 "8 \times 8")
    matrix
    ![A](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A "A")
    and the vector
    ![b](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b "b").

2.  Solve the system any way you want (using R) to find the unknown
    temperatures.

### Problem 3 Solution

Your solution goes here.

### Problem 4

The following command loads a matrix
![\\mathbf{M}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BM%7D "\mathbf{M}")
describing the western power network of the United States. Each element
is zero or one. A one in row
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i")
and column
![j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;j "j")
means that component
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i")
of the network is connected to component
![j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;j "j").
This type of matrix is called a network’s *adjacency matrix*.

``` r
M <- as.matrix(readMM(gzcon(url("https://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/bcspwr/bcspwr10.mtx.gz")))+0)
```

These commands visualize the network for you using a package called
`igraph`. It may take a little bit of time (several minutes) to run
because the graph is really large!

``` r
knitr::knit_hooks$set(crop = knitr::hook_pdfcrop)
g <- graph_from_adjacency_matrix(M, mode="undirected")
V(g)$color <- "grey"
V(g)$size <- 2
myLayout <- layout_nicely(g)
plot(simplify(g), layout=myLayout, vertex.color=V(g)$color, vertex.size=V(g)$size, vertex.label="")
```

<img src="Pset-Fundamentals-of-Linear-Systems_files/figure-gfm/unnamed-chunk-3-1.png" width="70%" />

1.  In the analysis of networks, one is often concerned with finding the
    most important (most central, in some sense) component. One measure
    of importance is called Katz centrality. The Katz centrality of the
    nodes,
    ![\\mathbf{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bx%7D "\mathbf{x}"),
    satisfies

    ![
    \\bigl(\\mathbf{I} - \\alpha M^T)(\\mathbf{x} + \\mathbf{1}) = \\mathbf{I}
    ](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbigl%28%5Cmathbf%7BI%7D%20-%20%5Calpha%20M%5ET%29%28%5Cmathbf%7Bx%7D%20%2B%20%5Cmathbf%7B1%7D%29%20%3D%20%5Cmathbf%7BI%7D%0A "
    \bigl(\mathbf{I} - \alpha M^T)(\mathbf{x} + \mathbf{1}) = \mathbf{I}
    ")

    where
    ![\\mathbf{I}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BI%7D "\mathbf{I}")
    is the identity matrix,
    ![\\alpha](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha "\alpha")
    is a small number, and
    ![\\mathbf{1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7B1%7D "\mathbf{1}")
    is a vector of all ones. Take
    ![\\alpha = 0.05](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha%20%3D%200.05 "\alpha = 0.05")
    and solve for
    ![\\mathbf{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bx%7D "\mathbf{x}")
    using any method you want. Print out the largest element of
    ![\\mathbf{x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bx%7D "\mathbf{x}")
    so that I can verify your answer.

2.  Re-plot the network using your same graph layout as before, color
    the most important node red and plot it with size 5. This way, the
    most important node (at least, according to the Katz measure) will
    stand out on your plot.

### Problem 4 Solution

Your solution goes here.
