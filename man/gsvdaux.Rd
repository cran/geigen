\name{gsvd.Auxiliaries}
\alias{gsvd.aux}
\title{Extract the R, D1, D2 matrices from a gsvd object}
\description{Returns a component of the object as a matrix}
\usage{
gsvd.R(z)
gsvd.oR(z)
gsvd.D1(z)
gsvd.D2(z)
}
\alias{gsvd.R}
\alias{gsvd.oR}
\alias{gsvd.D1}
\alias{gsvd.D2}
\arguments{
 \item{z}{an object created with \code{gsvd}}
}
\value{
\code{gsvd.R} returns the \code{R} matrix implied by the GSVD.

\code{gsvd.oR} returns the matrix \eqn{[0\,R]}{[0 R]} implied by the GSVD.

\code{gsvd.D1} returns the matrix \code{D1} implied by the GSVD.

\code{gsvd.D2} returns the matrix \code{D2} implied by the GSVD.

%\code{gsvd.X} returns the matrix
%\eqn{X = \left( \begin{array}{cc} I & 0 \\0 & R \end{array} \right) Q^T}{blockdiagonal(I,R) t(Q)}
%where the identity matrix \eqn{I} has size \eqn{n-k-l} (see \code{\link{gsvd}} for the meaning of \eqn{n,k} and \eqn{l}).
%This matrix is used to construct an alternative formulation of the Generalized Singular Value Decomposition.
}
\seealso{
\code{\link{gsvd}}
}
\examples{
A <- matrix(c(1,2,3,3,2,1,4,5,6,7,8,8), nrow=2, byrow=TRUE)
B <- matrix(1:18,byrow=TRUE, ncol=6)
A
B

z <- gsvd(A,B)
z

R <- gsvd.R(z)
oR <- gsvd.oR(z)
D1 <- gsvd.D1(z); D2 <- gsvd.D2(z)
R;oR
D1;D2
}
\keyword{algebra}
\keyword{array}
