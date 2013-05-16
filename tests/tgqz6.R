
library("geigen")
source("testqz.R")

set.seed(11)

n <- 5
A <- matrix(runif(n*n),nrow=n)+0i
B <- matrix(runif(n*n),nrow=n)+0i

# Test interface to zgges (QZ method)

z <- gqz(A, B,"-")
testqz(A,B,z)
gev <- z$alpha/z$beta
length(which(Re(gev)<0)) == z$sdim

z <- gqz(A, B,"+")
testqz(A,B,z)
gev <- z$alpha/z$beta
length(which(Re(gev)>0)) == z$sdim

z <- gqz(A, B,"S")
testqz(A,B,z)
gev <- z$alpha/z$beta
length(which(abs(gev)<1)) == z$sdim

z <- gqz(A, B,"B")
testqz(A,B,z)
gev <- z$alpha/z$beta
length(which(abs(gev)>1)) == z$sdim

z <- gqz(A, B,"R")
testqz(A,B,z)
gev <- z$alpha/z$beta
length(which(Im(gev)==0)) == z$sdim
