# Computational Linear Algebra

## Table of Contents

-   [Introduction](#introduction)
-   [Citation](#citation)
-   [In Development](#in-development)
-   [Errata](#errata)
-   [Consulting/Speaking](#consulting/speaking)
-   [Course Syllabus](#course-syllabus)
-   [R Resources](#r-resources)
-   [Course Details](#course-details)

## Course Units

-   [Getting Started](#getting-started)
-   [R Bootcamp](#r-bootcamp)
-   [How Computers Store Numbers](#how-computers-store-numbers)
-   [Linear Systems](#linear-systems)
-   [Solving Linear Systems](#solving-linear-systems)
-   [Recap Exercise I](#recap-exercise-I)
-   [Interpolation](#interpolation)
-   [Least Squares](#least-squares)
-   [Eigenvalues](#eigenvalues)
-   [Principal Component Analysis](#principal-component-analysis)
-   [Topological Data Analysis](#topological-data-analysis)
-   [Recap Exercise II](#recap-exercise-II)

## Introduction

I am [Prof. Chad Topaz](http://www.chadtopaz.com). Welcome to the public facing version of my Computational Linear Algebra course. This course is a minimalist recreation of the one I use on my own campus at Williams College. You can access most of the basic course materials here. Since this site is for public sharing, I have eliminated certain items related to surveying my students, performing course evaluation, student assessment, and so forth. I enthusiastically credit colleagues from my time at Macalester College: Tom Halverson, Danny Kaplan, David Shuman, and Lori Ziegelmeier. Their versions of a Computational Linear Algebra course have influenced my own course, and I'm grateful to them for sharing wisdom and materials.

[Back to Table of Contents](#table-of-contents)

## Citation

If you use this site as a student (for your own learning) or as an instructor (for teaching your own course) I would appreciate you [letting me know](mailto:chad.topaz+impact@gmail.com) so that I can track impact of the work I have done to create it. If you do adopt materials from this site, please make sure to credit me.

[Back to Table of Contents](#table-of-contents)

## In Development

This course site is perpetually in development. I regularly revisit my course materials in order to improve them.

[Back to Table of Contents](#table-of-contents)

## Errata

If you find any problems with course materials (a typo, code that doesn't work, a broken link – anything at all) please let me know by [contacting me here](mailto:chad.topaz+errata@gmail.com).

[Back to Table of Contents](#table-of-contents)

## Consulting/Speaking

I have built this public facing site as a free, shareable public good. If you are interested in assistance with technical/mathematical aspects of this course or related consulting, please [contact me here](mailto:chad.topaz+consulting@gmail.com). If you would like me to speak to your department/institution/organization about curriculum and/or pedagogy (of this course, or of anything else) then please see my [Speaking](http://www.chadtopaz.com/speaking) page and contact me [here](mailto:chad.topaz+speaking@gmail.com).

[Back to Table of Contents](#table-of-contents)

## Course Syllabus

Rather than using a traditional text-based course syllabus, I use a rather concise graphical syllabus; see below. I designed this syllabus using Piktochart. A similar tool you could consider is Canva. Some details of this syllabus are not relevant for a broad public audience, but I include the full syllabus nonetheless in order to give a flavor of the course and the philosophy I bring to it.

![Computational Linear Algebra Course Syllabus](https://github.com/chadtopaz/computationallineaaralgebra/raw/main/syllabus/Computational%20Linear%20Algebra%20Syllabus.jpg)

[Back to Table of Contents](#table-of-contents)

## R Resources

This course uses the R programming language in the RStudio integrated development environment. [RStudio Cloud](https://rstudio.cloud) is a free version of RStudio that you can run on any device that has an internet connection and a web browser.

Here are some help resources for R.

-   Learning R

    -   [DataCamp](https://www.datacamp.com/)
    -   [DataQuest](http://dataquest.io/)
    -   [Learn R in Y Minutes](https://learnxinyminutes.com/docs/r/)
    -   [RStudio Primers](https://rstudio.cloud/learn/primers)
    -   [R for Data Science](https://r4ds.had.co.nz/index.html)

-   Cheat Sheets

    -   [RStudio Environment](https://raw.githubusercontent.com/rstudio/cheatsheets/main/rstudio-ide.pdf)
    -   [Base R](http://github.com/rstudio/cheatsheets/blob/main/base-r.pdf)
    -   [R Markdown](https://raw.githubusercontent.com/rstudio/cheatsheets/main/rmarkdown.pdf)
    -   [Data Visualization](https://raw.githubusercontent.com/rstudio/cheatsheets/main/data-visualization.pdf)

-   Help Resources

    -   [Stack Overflow](http://stackoverflow.com)
    -   [My Favorite Choice](https://letmegooglethat.com/?q=for+loop+in+R)

[Back to Table of Contents](#table-of-contents)

## Course Details

A few characteristics of my on-the-ground course include:

-   Ungrading
-   Flipped classroom style, including pre-class quizzes
-   Introductory activities that get students oriented to the course
-   Organization around weekly blocks of material
-   Substantial in-class small-group discussion

To enable you to make easy and flexible use of the resources here, I have eliminated most of these elements from this public facing version of the course. Essentially, what I am providing here is a core dump of materials I use. The materials are linked throughout the various course units. Additionally, you can access the raw mateerials in this respository's subdirectories:

-   [/coursenotes](/coursenotes)
-   [/demos](/demos)
-   [/activities](/activities)
-   [/psets](psets)

[Back to Table of Contents](#table-of-contents)

## Getting Started

-   Videos

    -   [Why is Linear Algebra Useful (9:56)](https://www.youtube.com/watch?v=X0HXnHKPXSo&feature=youtu.be)

-   Reading

    -   [Why I Don't Grade, by Jesse Stommel](https://www.jessestommel.com/why-i-dont-grade/)
    -   [The Mathematics of Mass Testing for COVID-19](https://sinews.siam.org/Details-Page/the-mathematics-of-mass-testing-for-covid-19)
    -   [DataCamp Hid a Sexual Harassment Incident by its CEO](https://www.buzzfeednews.com/article/daveyalba/datacamp-sexual-harassment-metoo-tech-startup) [CW/TW: Sexual Harassment]. DataCamp is this course's default resource for learning R because it is the best-designed and most effective resource I know of. However, **DataCamp has a problematic past**, and I feel an ethical imperative to be transparent about this, which is why I am having you read the article linked above. If you feel you cannot or should not use DataCamp, I encourage you use an alternative option, including one of the choices I have listed in the [R Resources](#r-resources) section. In case it impacts your decision, please know that DataCamp does not gain financially from us using it, as instructors at colleges/universities can get six free months of DataCamp for their students. In any case, you are 100% free to make whatever choice you want and I will support you fully!

-   DataCamp Coding Practice

    -   [Introduction to R (course)](https://www.datacamp.com/courses/free-introduction-to-r)

## R Bootcamp

-   Videos

    -   [Horner's Method (17:16)](https://www.youtube.com/watch?v=RGrmEWj38bs) 
    -   [Intro to Taylor Series (12:42)](https://youtu.be/9YAaCEA08yM)
    -   [Why Taylor Series Actually Work (9:34)](https://youtu.be/Iub16Y1ZXRw)

-   Reading

    -   [Course Notes: R Bootcamp](/coursenotes/coursenotes.md#r-bootcamp)
    -   [Machine Learning Has Been Used to Automatically Translate Long-Lost Languages](https://www.technologyreview.com/2019/07/01/65601/machine-learning-has-been-used-to-automatically-translate-long-lost-languages/)
    -   [The Best Bits](https://www.americanscientist.org/article/the-best-bits)

-   DataCamp Coding Practice

    -   [Getting Started with R Markdown (chapter)](https://www.datacamp.com/courses/reporting-with-rmarkdown)
    -   [Introduction to Linear Algbra (chapter)](https://www.datacamp.com/courses/linear-algebra-for-data-science-in-r)
    -   [Conditionals and Control Flow (chapter)](https://www.datacamp.com/courses/intermediate-r)
    -   [Loops (chapter)](https://www.datacamp.com/courses/intermediate-r)
    -   [Functions (chapter)](https://www.datacamp.com/courses/intermediate-r)

-   Activity - R Bootcamp [[.md for viewing](activities/Activity-R-Bootcamp.md)] [[.Rmd for editing](activities/Activity-R-Bootcamp.Rmd)]

-   Pset - R Bootcamp [[.md for viewing](psets/Pset-R-Bootcamp.md)] [[.Rmd for editing](psets/Pset-R-Bootcamp.Rmd)]

-   Challenge Problems

    -   Many students take the formula for Taylor series as a given, but it's valuable to understand how it comes from calculus. Derive the Taylor series formula [hint: start with the Fundamental Theorem of Calculus].

-   Looking ahead

    -   8.2 - 7.2 - 1 = ?

[Back to Course Units](#course-units)

## How Computers Store Numbers

-   Videos

    -   [Floating Point Numbers (17:29)](https://youtu.be/gc1Nl3mmCuY)
    -   [Floating Point Representation and Rounding Error (16:26)](https://youtu.be/wbxSTxhTmrs)

-   Reading

    -   [Course Notes: How Computers Store Numbers](/coursenotes/coursenotes.md#how-computers-store-numbers)
    -   [Bits and Bugs](https://drive.google.com/file/d/1Ez-Vjk9e97o5lMMsSkFLFZF968q7rjyC/view?usp=sharing) 2.1 - 2.4

-   DataCamp Coding Practice

    -   [Quick Introduction to Base R Graphics (chapter)](https://www.datacamp.com/courses/data-visualization-in-r)

-   Activity - How Computers Store Numbers [[.md for viewing](activities/Activity-How-Computers-Store-Numbers.md)] [[.Rmd for editing](activities/Activity-How-Computers-Store-Numbers.Rmd)]

-   Pset - How Computers Store Numbers [[.md for viewing](psets/Pset-How-Computers-Store-Numbers.md)] [[.Rmd for editing](psets/Pset-How-Computers-Store-Numbers.Rmd)]

[Back to Course Units](#course-units)

## Linear Systems

[Back to Course Units](#course-units)

## Solving Linear Systems

[Back to Course Units](#course-units)

## Recap Exercise I

[Back to Course Units](#course-units)

## Interpolation

[Back to Course Units](#course-units)

## Least Squares I

[Back to Course Units](#course-units)

## Least Squares II

[Back to Course Units](#course-units)

## Eigenvalues

[Back to Course Units](#course-units)

## Principal Component Analysis

[Back to Course Units](#course-units)

## Topological Data Analysis

[Back to Course Units](#course-units)

## Recap Exercise II

[Back to Course Units](#course-units)
