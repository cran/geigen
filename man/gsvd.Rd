\name{gsvd}
\alias{gsvd}
\encoding{UTF-8}
\title{Generalized Singular Value Decomposition}
\description{Computes the generalized singular value decomposition of a pair of matrices.}
\usage{
gsvd(A,B)
}
\arguments{
 \item{A}{a numeric matrix with \eqn{m}{m} rows and \eqn{n}{n} columns.}
 \item{B}{a numeric matrix with \eqn{p}{p} rows and \eqn{n}{n} columns.}
}
\details{
The matrix A is a \eqn{m}-by-\eqn{n} matrix and the matrix B is a \eqn{p}-by-\eqn{n} matrix.
This function decomposes both matrices.
The Generalized Singular Value Decomposition of the matrices \eqn{A} and \eqn{B}
is given as
\deqn{A = U D_1 [0\, R] Q^T}{A = U * D_1 * [0 R] * t(Q)} and
\deqn{B = V D_2 [0\, R] Q^T}{B = V * D_2 * [0 R] * t(Q)}
where
\describe{
\item{\eqn{U}}{an \eqn{m\times m}{m-by-m} orthogonal matrix.}
\item{\eqn{V}}{a \eqn{p\times p}{p-by-p} orthogonal matrix.}
\item{\eqn{Q}}{an \eqn{n\times n}{n-by-n} orthogonal matrix.}
\item{\eqn{R}}{an \eqn{r}-by-\eqn{r} upper triangular non singular matrix and the matrix
\eqn{[0\, R]} is an \eqn{r}-by-\eqn{n} matrix. The quantity \eqn{r} is the rank of
the matrix \eqn{\left( \begin{array}{c} A \\B \end{array} \right)}{rbind(A,B)} with \eqn{r \le n}{r <= n}.}
\item{\eqn{D_1},\eqn{D_2}}{are quasi diagonal matrices and nonnegative and satisfy
\eqn{D_1^T D_1 + D_2^T D_2 = I}{t(D_1)\%*\% D_1 + t(D_2) \%*\% D_2}.
\eqn{D_1} is an \eqn{m}-by-\eqn{r} matrix and \eqn{D_2} is a \eqn{p}-by-\eqn{r} matrix.}
}
For details on this decomposition and the structure of the matrices \eqn{D_1} and \eqn{D_2}
see \url{http://www.netlib.org/lapack/lug/node36.html}.
}
\value{
The return value is a list containing the following components
\describe{
  \item{\code{A}}{the upper triangular matrix or a part of \eqn{R}.}
  \item{\code{B}}{lower part of the triangular matrix \eqn{R} if
            \eqn{k+l>m} (see below).}
  \item{\code{m}}{number of rows of \eqn{A}.}
  \item{\code{k}}{\eqn{r{-}l}{r-l}. The number of rows of the matrix \eqn{R} is \eqn{k{+}l}{k+l}.
                  The first \eqn{k} generalized singular values are infinite.}
  \item{\code{l}}{effective rank of the input matrix \eqn{B}.
                  The number of finite generalized singular values after the first \eqn{k} infinite ones.}
  \item{\code{alpha}}{a numeric vector with length \eqn{n} containing
                      the numerators of the generalized singular values in entries
                      \eqn{1} \ldots \eqn{k{+}l}{k+l}.}
  \item{\code{beta}}{a numeric vector with length \eqn{n} containing
                     the denominators of the generalized singular value in entries
                     \eqn{1} \ldots \eqn{k{+}l}{k+l}.}
  \item{\code{U}}{the orthogonal matrix \eqn{U}.}
  \item{\code{V}}{the orthogonal matrix \eqn{V}.}
  \item{\code{Q}}{the orthogonal matrix \eqn{Q}.}
}
For a detailed description of these items see \url{http://www.netlib.org/lapack/lug/node36.html}.
Auxiliary functions are provided for extraction and manipulation of the various items.
}
\source{
  \code{gsvd} uses the LAPACK routine \code{DGGSVD}.
  LAPACK is from \url{http://www.netlib.org/lapack}.
  The decomposition is fully explained in \url{http://www.netlib.org/lapack/lug/node36.html}.
}
\references{
  Anderson. E. and ten others (1999)
  \emph{LAPACK Users Guide}. Third Edition. SIAM.\cr
  Available on-line at
  \url{http://www.netlib.org/lapack/lug/lapack_lug.html}.
  See the section \emph{Generalized Eigenvalue and Singular Value Problems}
  (\url{http://www.netlib.org/lapack/lug/node33.html}) and
  and the section \emph{Generalized Singular Value Decomposition (GSVD)}
  (\url{http://www.netlib.org/lapack/lug/node36.html}).
}
\seealso{
\link{gsvd.aux}
}
\examples{
A <- matrix(c(1,2,3,3,2,1,4,5,6,7,8,8), nrow=2, byrow=TRUE)
B <- matrix(1:18,byrow=TRUE, ncol=6)
A
B

z <- gsvd(A,B)
z
}
\keyword{algebra}
\keyword{array}
