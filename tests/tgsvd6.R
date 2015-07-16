
library("geigen")
source("testgsvd.R")

set.seed(11)
A <- matrix(round(runif(36),2),byrow=TRUE, ncol=6)

B <- matrix(round(runif(36),2),byrow=TRUE, ncol=6)

A
B

z <- gsvd(A,B)
testgsvd(z,A,B)
