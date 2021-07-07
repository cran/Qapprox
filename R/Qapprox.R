#' Right-tail probability of quadratic forms of centered Gaussian variables.
#' @param q - quantile, could be a vector.
#' @param Sigma - covariance matrix of Gaussian variables.
#' @param A - a positive-semi-definite matrix that defines the quadratic form.
#' @param method - "MR": moment-ratio (skewness-kurtosis) matching method; "SW": Satterthwaite-Welch method that matches mean and variance; "LTZ4": Liu-Tang-Zhang method that matches the kurtosis.
#' @return The right-tail probability of a quadratic form (Q = X'AX) of centered Gaussian variables.
#' @references 1. Hong Zhang, Judong Shen and Zheyang Wu. "An efficient and accurate approximation to the distribution of quadratic forms of Gaussian variables", arXiv:2005.00905.
#' @examples
#' n <- 100
#' Sigma <- toeplitz(1/(1:n))
#' thr <- 180
#' Qapprox(thr, Sigma, method="SW")
#' Qapprox(thr, Sigma, method="LTZ4")
#' Qapprox(thr, Sigma, method="MR")
#' @export
#' @importFrom stats pchisq pgamma

Qapprox <- function(q, Sigma, A = NULL, method="MR") {
  if(!is.null(A)){
    Sigma <- A%*%Sigma
  }

  c1 <- sum(diag(Sigma)) # mean
  c2 <- 2*sum(Sigma*Sigma) # variance

  if(method == "SW"){
    a <- c1^2/c2
    s <- c1/a
    pval <- pgamma(q, shape = a, scale = s, lower.tail = FALSE)
  }else{
    Sigma2 <- Sigma%*%Sigma
    c3 <- 8*sum(Sigma*Sigma2)
    c4 <- 48*sum(Sigma2*Sigma2)
    gm <- c3/sqrt(c2)^3 # skewness
    kp <- c4/c2^2 + 3 # kurtosis
    if(method == "MR"){
      a <- 9*gm^2/(kp-3)^2
      pval <- pgamma((q-c1)/sqrt(c2)*sqrt(a)+a, shape=a,  scale = 1, lower.tail = FALSE)
    }else if(method == "LTZ4"){
      s1 = gm/sqrt(8); s2 = (kp - 3)/12
      D = s1^2 - s2
      if(D > 0){
        a <- 1/(s1-sqrt(D))
        ncp <- s1*a^3-a^2
        df <- a^2 - 2*ncp
        pval <- pchisq((q-c1)/sqrt(c2)*sqrt(2)*a+df+ncp, df=df, ncp=ncp, lower.tail = FALSE)
      }else{
        a <- 1/sqrt(s2)
        ncp <- 0
        df <- a^2
        pval <- pchisq((q-c1)/sqrt(c2)*sqrt(2)*a+df, df=df, lower.tail = FALSE)
      }
    }
  }
  pval
}
