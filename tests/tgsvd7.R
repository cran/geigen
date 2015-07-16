
library("geigen")
source("testgsvd.R")

set.seed(11)
A <- matrix(round(runif(18),2), nrow=3, byrow=TRUE)

B <- matrix(round(runif(24),2),byrow=TRUE, ncol=6)

A
B

z <- gsvd(A,B)
testgsvd(z,A,B)
