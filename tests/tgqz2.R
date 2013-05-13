
library("geigen")
source("testqz.R")

set.seed(11)

n <- 5
A <- matrix(runif(n*n),nrow=n)
B <- matrix(runif(n*n),nrow=n)

# Test interface to dgges (QZ method)

z <- gqz(A, B,"-")
testqz(A,B,z)
gev <- complex(real=z$alphar,imaginary=z$alphai)/z$beta
length(which(Re(gev)<0)) == z$sdim

z <- gqz(A, B,"+")
testqz(A,B,z)
gev <- complex(real=z$alphar,imaginary=z$alphai)/z$beta
length(which(Re(gev)>0)) == z$sdim

z <- gqz(A, B,"S")
testqz(A,B,z)
gev <- complex(real=z$alphar,imaginary=z$alphai)/z$beta
length(which(abs(gev)<1)) == z$sdim

z <- gqz(A, B,"B")
testqz(A,B,z)
gev <- complex(real=z$alphar,imaginary=z$alphai)/z$beta
length(which(abs(gev)>1)) == z$sdim

z <- gqz(A, B,"R")
testqz(A,B,z)
gev <- complex(real=z$alphar,imaginary=z$alphai)/z$beta
length(which(Im(gev)==0)) == z$sdim
