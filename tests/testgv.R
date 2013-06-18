testgv <- function(A,B,z) {
    V <- z$vectors 
    ret <- A %*% V - B %*% V %*% diag(z$values)
    if( !is.null(z$alpha) )
        ret <- c(ret, A %*% V %*% diag(z$beta) - B %*% V %*% diag(z$alpha))
    tol <- 100 * sqrt(.Machine$double.eps)
    all(abs(ret)<=tol)
}
