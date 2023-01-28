Load necessary packages:

    library(Matrix)
    library(tidyverse)

### Problem 1

By completing this problem, you’ll hone your PCA skills in R. We’re
going to study nutritional information of food. Let’s load some data
from the USDA National Nutrient Database.

    food <- read.csv("https://query.data.world/s/ixxtapv4y3rbkb6xo5c32ipljmdkor")

a\. Perform a principal component analysis based on the macronutrient
info, that is, protein, fat, carb, sugar, and fiber. What percent of
variation in the data do the first two principal components capture?

b\. For ease of bookkeeping, create a new data frame that has columns
`FoodGroup` and `ShortDescrip` from the original data set, and `PC1` and
`PC2` from your principal component analysis. Plot the data in
two-dimensional principal component space.

c\. That’s a lot of data! Just for funsies, subset your data to include
only the “Beef Products”,“Fruits and Fruit Juices”, and “Fats and Oils”
food groups. Plot these in two-dimensional principal component space,
colored by food group. If your data frame was called `myfood` you would
use the command

    ggplot(myfood, aes(x = PC1, y = PC2, color = FoodGroup)) +
      geom_point()

Interpret what you see.

### Problem 1 Solution

a\. Your solution goes here.

b\. Your solution goes here.

c\. Your solution goes here.