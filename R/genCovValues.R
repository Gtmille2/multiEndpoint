#' Generate covariate Values
#'
#' This is a function generate covariate values to use before an allocation scheme
#' @param p1 The probability vector for covariate 1
#' @param p2 The probability vector for in covariate 2
#' @param N The total sample size
#' @export
#' @examples
#' genCovValues()

genCovValues = function(p = c(0.5,0.5), N = 200) {
  full = NULL
  for (i in 1:length(p)) full = cbind(full,sample(2,N,p[i]))
  full
  # matrix(c(z1,z2),ncol=2)
}
