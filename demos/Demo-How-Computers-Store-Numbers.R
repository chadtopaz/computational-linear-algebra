# Load libraries
library(kableExtra)

# Write numbers out as bit representation
x <- 97
numToBits(x)

# That's annoying! Let's write a little function to
# make more friendly
bitRep <- function(x) {
  x <- numToBits(x)
  x <- as.numeric(x)
  x <- rev(x)
  paste0(c(x[1], " | ", x[2:12], " | ", x[13:64]), collapse = "")
}

# Let's test it
bitRep(x)

# Check our work
(-1)^0 * 2^(1 + 4 + 1024 - 1023) * (1 + 1/2 + 1/64)

# Loss of significance
options(digits = 12)
x <- 10^(-(0:12))
E1 <- (1 - cos(x))/sin(x)^2
E2 <- 1/(1 + cos(x))
kable(cbind(x, E1, E2), format = "simple")