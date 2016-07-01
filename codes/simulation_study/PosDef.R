# This function for generating positive definite matrices was taken from
#https://stat.ethz.ch/pipermail/r-help/2008-February/153708.html
# copyright Ravi Varadhan
Posdef <- function (n, ev = runif(n, 0, 10)) {
     Z <- matrix(ncol=n, rnorm(n^2))
     decomp <- qr(Z)
     Q <- qr.Q(decomp) 
     R <- qr.R(decomp)
     d <- diag(R)
     ph <- d / abs(d)
     O <- Q %*% diag(ph)
     Z <- t(O) %*% diag(ev) %*% O
     return(Z)
}
