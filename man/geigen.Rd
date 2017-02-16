\name{geigen}
\alias{geigen}
\encoding{UTF-8}
\title{Generalized Eigenvalues}
\description{Computes generalized eigenvalues and eigenvectors of a pair of matrices.}
\usage{
geigen(A, B, symmetric, only.values=FALSE)
}
\arguments{
 \item{A}{left hand side matrix.}
 \item{B}{right hand side matrix.}
 \item{symmetric}{if \code{TRUE}, both matrices are assumed to be symmetric (or Hermitian if complex)
 and only their lower triangle (diagonal included) is used.
 If \code{symmetric} is not specified, the matrices are inspected for symmetry.}
 \item{only.values}{if \code{TRUE} only eigenvalues are computed otherwise both eigenvalues and eigenvctors are
returned.}
}
\details{
If the argument \code{symmetric} is missing, the function
will try to determine if the matrices are symmetric with the function \code{isSymmetric} from
the \pkg{base} package. It is faster to specify the argument.

Both matrices must be square.
This function provides the solution to the generalized eigenvalue problem defined by
\deqn{A x = \lambda Bx}{A*x = lambda B*x}
If either one of the matrices is complex the other matrix is coerced to be complex.
If the matrices are complex and Hermitian then the matrix \code{B} must be positive definite.

For general matrices the generalized eigenvalues \eqn{\lambda}{lambda}
are calculated as the ratio \eqn{\alpha / \beta}{alpha/beta}
where \eqn{\beta}{beta} may be zero or very small leading
to non finite or very large values for the eigenvalues.
Therefore the values for \eqn{\alpha}{alpha} and \eqn{\beta}{beta} are also included in the return value
of the function.
When both matrices are complex (or coerced to be so) the generalized eigenvalues,
\eqn{\alpha}{alpha} and \eqn{\beta}{beta} are complex.
When both matrices are numeric \eqn{\alpha}{alpha} may be numeric or complex and
\eqn{\beta}{beta} is numeric.

When both matrices are symmetric (or Hermitian) the generalized eigenvalues are numeric and
no components \eqn{\alpha}{alpha} and \eqn{\beta}{beta} are available.}
\value{
A list containing components
  \item{values}{a vector containing the \eqn{n} generalized eigenvalues.}
  \item{vectors}{an \eqn{n\times n}{n-by-n} matrix containing the generalized eigenvectors or NULL
  if \code{only.values} is \code{TRUE}.}
  \item{alpha}{the numerator of the generalized eigenvalues and may be NULL if not applicable.}
  \item{beta}{the denominator of the generalized eigenvalues and may be NULL if not applicable.}
}

\source{
  \code{geigen} uses the LAPACK routines \code{DGGEV},
  \code{DSYGV}, \code{ZHEGV} and \code{ZGGEV}.
  LAPACK is from \url{http://www.netlib.org/lapack}.
  The complex routines used by the package come from LAPACK 3.7.0.
}
\references{
  Anderson. E. and ten others (1999)
  \emph{LAPACK Users' Guide}.  Third Edition.  SIAM.\cr%'
  Available on-line at
  \url{http://www.netlib.org/lapack/lug/lapack_lug.html}.
  See the section \emph{Generalized Eigenvalue and Singular Value Problems}
  (\url{http://www.netlib.org/lapack/lug/node33.html}).
}
\seealso{
  \code{eigen}
}

\examples{
A <- matrix(c(14, 10, 12,
              10, 12, 13,
              12, 13, 14), nrow=3, byrow=TRUE)

B <- matrix(c(48, 17, 26,
              17, 33, 32,
              26, 32, 34), nrow=3, byrow=TRUE)

z1 <- geigen(A, B, symmetric=FALSE, only.values=TRUE)
z2 <- geigen(A, B, symmetric=FALSE, only.values=FALSE )
z2

# geigen(A, B)
z1 <- geigen(A, B, only.values=TRUE)
z2 <- geigen(A, B, only.values=FALSE)
z1;z2

A.c <- A + 1i
B.c <- B + 1i

A[upper.tri(A)] <- A[upper.tri(A)] + 1i
A[lower.tri(A)] <- Conj(t(A[upper.tri(A)]))

B[upper.tri(B)] <- B[upper.tri(B)] + 1i
B[lower.tri(B)] <- Conj(t(B[upper.tri(B)]))

isSymmetric(A)
isSymmetric(B)

z1 <- geigen(A, B, only.values=TRUE)
z2 <- geigen(A, B, only.values=FALSE)
z1;z2
}
\keyword{algebra}
\keyword{array}
