testgv <- function(A,B,z) {
    V <- z$vectors
    ret <- A %*% V - B %*% V %*% diag(z$values)
    tol <- 100 * sqrt(.Machine$double.eps)
    all(abs(ret)<=tol)
}
