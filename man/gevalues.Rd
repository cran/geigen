\name{gevalues}
\alias{gevalues}
\encoding{UTF-8}
\title{Calculate generalized eigenvalues from a Generalized Schur decomposition}
\description{Computes the generalized eigenvalues from an object constructed with \code{gqz}.}
\usage{
gevalues(x)
}
\arguments{
 \item{x}{an object created with \code{gqz}.}
}
\details{The function calculates the generalized eigenvalues from elements of the object returned by
the function \code{gqz}.
The generalized eigenvalues are computed from a ratio where the denominator 
(the \eqn{\beta}{beta} component of the argument) may be zero. 
The function attempts to guard against nonsensical complex \code{NaN} values when dividing by zero 
which may happen on some systems.
}
\value{A vector containing the generalized eigenvalues.
The vector is numeric if the imaginary parts of the eigenvalues are all zero and complex otherwise.
} 
\seealso{
  \code{geigen}, \code{gqz}
}
\examples{
# Real matrices
# example from NAG: http://www.nag.com/lapack-ex/node116.html
# Find the generalized Schur decomposition with the real eigenvalues ordered to come first

A <- matrix(c(  3.9, 12.5,-34.5,-0.5,
                4.3, 21.5,-47.5, 7.5,
                4.3, 21.5,-43.5, 3.5,
                4.4, 26.0,-46.0, 6.0), nrow=4, byrow=TRUE)

B <- matrix(c( 1.0, 2.0, -3.0, 1.0,
               1.0, 3.0, -5.0, 4.0,
               1.0, 3.0, -4.0, 3.0,
               1.0, 3.0, -4.0, 4.0), nrow=4, byrow=TRUE)

z <- gqz(A, B,"R")
z
# compute the generalized eigenvalues
ger <- gevalues(z)
ger
}
\keyword{algebra}
\keyword{array}
