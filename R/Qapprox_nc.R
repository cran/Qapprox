#' Right-tail probability of quadratic forms (Q = X'AX) of noncentral Gaussian variables.
#' @param q - quantile, could be a vector.
#' @param mu - mean vector of Gaussian variables.
#' @param Sigma - covariance matrix of Gaussian variables.
#' @param A - a positive-semi-definite matrix that defines the quadratic form.
#' @param method - "MR": moment-ratio (skewness-kurtosis) matching method; "SW": Satterthwaite-Welch method that matches mean and variance; "LTZ4": Liu-Tang-Zhang method that matches the kurtosis.
#' @return The right-tail probability of a quadratic form (Q = X'AX) of noncentral Gaussian variables.
#' @references 1. Hong Zhang, Judong Shen and Zheyang Wu. "An efficient and accurate approximation to the distribution of quadratic forms of Gaussian variables", arXiv:2005.00905.
#' @examples
#' n = 100
#' Sigma = toeplitz(1/(1:n))
#' mu = rep(1, n)
#' thr = 500
#' Qapprox_nc(thr, mu, Sigma, method="SW")
#' Qapprox_nc(thr, mu, Sigma, method="LTZ4")
#' Qapprox_nc(thr, mu, Sigma, method="MR")
#' @export
#' @importFrom stats pchisq pgamma

Qapprox_nc <- function(q, mu, Sigma, A = NULL, method="MR"){
  if(is.null(A)){
    A = diag(1, dim(Sigma))
  }
  Sigma_chol = chol(Sigma)
  ei = eigen(Sigma_chol%*%A%*%t(Sigma_chol))
  lam = ei$values; P = t(ei$vectors)

  muy = P%*%solve(t(Sigma_chol))%*%mu
  delta = muy^2

  CM = sapply(1:4, function(k)sum(lam^k)+k*sum(lam^k*delta))
  c1 = CM[1]; c2 = 2*CM[2];
  s1 = CM[3]/(CM[2]^(3/2)); s2 = CM[4]/CM[2]^2
  gm = sqrt(8)*s1; kp = 12*s2 + 3
  D = s1^2-s2

  if(method=="SW"){
    a = c1^2/c2
    s = c1/a
    pval = pgamma(q, shape = a, scale = s, lower.tail = FALSE)
  }else if(method=="MR"){
    a = 9*gm^2/(kp-3)^2
    pval  = pgamma((q-c1)/sqrt(c2)*sqrt(a)+a, shape=a,  scale = 1, lower.tail = FALSE)
  }else if(method=="LTZ4"){
    if(D>0){
      a = 1/(s1-sqrt(D)); ncp = s1*a^3-a^2; df = a^2 - 2*ncp
      pval = pchisq((q-c1)/sqrt(c2)*sqrt(2)*a+df+ncp, df=df, ncp=ncp, lower.tail = FALSE)
    }else{
      a = 1/sqrt(s2); ncp = 0; df = a^2
      pval = pchisq((q-c1)/sqrt(c2)*sqrt(2)*a+df, df=df, lower.tail = FALSE)
    }
  }
  return(pval)
}
