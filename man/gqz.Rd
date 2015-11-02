\name{gqz}
\alias{gqz}
\encoding{UTF-8}
\title{Generalized Schur decomposition}
\description{Computes the generalized eigenvalues and Schur form of a pair of matrices.}
\usage{
gqz(A, B, sort=c("N","-","+","S","B","R"))
}
\arguments{
 \item{A}{left hand side matrix.}
 \item{B}{right hand side matrix.}
 \item{sort}{how to sort the generalized eigenvalues. See \sQuote{Details}.}
}
\details{
Both matrices must be square.
This function provides the solution to the generalized eigenvalue problem defined by
\deqn{A x = \lambda Bx}{A*x = lambda B*x}
If either one of the matrices is complex the other matrix is coerced to be complex.
If the matrices are complex and Hermitian then the matrix \code{B} must be positive definite.

The \code{sort} argument specifies how to order the eigenvalues on the
diagonal of the generalized Schur form, where it is noted that non-finite eigenvalues never
satisfy any ordering condition (even in the case of a complex infinity).
Eigenvalues that are placed in the leading block of the Schur form
satisfy
\describe{
  \item{\code{N}}{unordered.}
  \item{\code{-}}{negative real part.}
  \item{\code{+}}{positive real part.}
  \item{\code{S}}{absolute value < 1.}
  \item{\code{B}}{absolute value > 1.}
  \item{\code{R}}{imaginary part identical to 0 with a tolerance of 100*machine_precision as determined by Lapack.}
}
}
\value{
The generalized Schur form for \emph{numeric} matrices is
\deqn{ (A,B) = (Q S Z^T, Q T Z^T)}{(A,B) = (Q*S*Z^T, Q*T*Z^T)}
The matrices \eqn{Q}{Q} and \eqn{Z}{Z} are orthogonal. The matrix \eqn{S}{S} is quasi-upper triangular and
the matrix \eqn{T}{T} is upper triangular.
The return value is a list containing the following components
\describe{
  \item{\code{S}}{generalized Schur form of A.}
  \item{\code{T}}{generalized Schur form of B.}
  \item{\code{sdim}}{the number of eigenvalues (after sorting) for which the sorting condition is true.}
  \item{\code{alphar}}{numerator of the real parts of the eigenvalues (numeric).}
  \item{\code{alphai}}{numerator of the imaginary parts of the eigenvalues (numeric).}
  \item{\code{beta}}{denominator of the expression for the eigenvalues (numeric).}
  \item{\code{Q}}{matrix of left Schur vectors (matrix Q).}
  \item{\code{Z}}{matrix of right Schur vectors (matrix Z).}
}

The generalized Schur form for \emph{complex} matrices is
\deqn{(A,B) = (Q S Z^H, Q T Z^H)}{(A,B) = (Q*S*Z^H, Q*T*Z^H)}
The matrices \eqn{Q}{Q} and \eqn{Z}{Z} are unitary and the matrices \eqn{S}{S} and
\eqn{T}{T} are upper triangular.
The return value is a list containing the following components
\describe{
  \item{\code{S}}{generalized Schur form of A.}
  \item{\code{T}}{generalized Schur form of B.}
  \item{\code{sdim}}{the number of eigenvalues. (after sorting) for which the sorting condition is true.}
  \item{\code{alpha}}{numerator of the eigenvalues (complex).}
  \item{\code{beta}}{denominator of the eigenvalues (complex).}
  \item{\code{Q}}{matrix of left Schur vectors (matrix Q).}
  \item{\code{Z}}{matrix of right Schur vectors (matrix Z).}
}

The generalized eigenvalues can be computed by calling function \code{gevalues}.
}
\source{
  \code{gqz} uses the LAPACK routines \code{DGGES} and \code{ZGGES}.
  LAPACK is from \url{http://www.netlib.org/lapack}.
  The complex routines used by the package come from LAPACK version 3.4.2.
}
\references{
  Anderson. E. and ten others (1999)
  \emph{LAPACK Users' Guide}.  Third Edition.  SIAM.\cr%'
  Available on-line at
  \url{http://www.netlib.org/lapack/lug/lapack_lug.html}.
  See the section \emph{Eigenvalues, Eigenvectors and Generalized Schur Decomposition}
  (\url{http://www.netlib.org/lapack/lug/node56.html}).
}
\seealso{
  \code{geigen}, \code{gevalues}
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

# complexify
A <- A+0i
B <- B+0i
z <- gqz(A, B,"R")
z
}
\keyword{algebra}
\keyword{array}
