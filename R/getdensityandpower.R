#'get_density_first
#'
#'This function gets the density of the test statistic z at the first analysis
#'@param ntmt Ntmt is the number of treatments
#'@param theta Theta is the treatment effect variable
#'@param v V is the value of the fisher information statistic
#'@param z Z is the current value for the Efficient Score statistic
get_density_first = function(ntmt,theta,v,z) {
  lowerLimit = -20*sqrt(v)
  upperLimit = 20*sqrt(v)

  xvalweights = gaussLegendre(INTEGR_POINTS_First,lowerLimit,upperLimit)
  integral = 0
  for ( j in 0 :(INTEGR_POINTS_First-1)){
    integral = integral + get_integrand(ntmt,theta,v,z,xvalweights$x[j+1])*xvalweights$w[j+1]
  }
  return(integral)
}

#'get_integrand
#'
#'This is a function that returns the integrand for the
#'@param k K is the number of treatments used in this trial before the trial
#'@param theta Theta is the treament effect
#'@param v V is the current fisher's information
#'@param z Z is the current
#'@param x X is used in the integrand
get_integrand = function(k,theta,v,z,x) {
  root = sqrt(v/2)
  product = 1
  for( i in 1:k ) {
    product = product*pnorm((x-theta1[i]*v)/root)

  }
  integrand = product*2*(x-z)*dnorm((x-z)/root)/(v*root)
  return(integrand)
}

#' get_power_first
#'
#' This function gets the power for the first test statistic
#' @param ntmt Ntmt is the number of treatments
#' @param theta Theta is the treatment effect
#' @param v V is the Fisher information test statistic
#' @param z Z is the current Efficient Score test statistic
get_power_first = function(ntmt,theta,v,z) {
  lowerLimit = -20*sqrt(v)
  upperLimit = 20*sqrt(v)

  xvalweights = gaussLegendre(INTEGR_POINTS_First,lowerLimit,upperLimit)

  integral = 0
  for ( j in 0:(INTEGR_POINTS_First-1)) {
    integral = integral + get_integrand_power(ntmt,theta,v,z,xvalweights$x[j+1])*xvalweights$w[j+1]
  }
  return(integral)
}

#' get_integrand_power
#'
#' This function gets the power for the first test statistic
#' @param ntmt Ntmt is the number of treatments
#' @param theta Theta is the treatment effect
#' @param v V is the Fisher information test statistic
#' @param z Z is the current Efficient Score test statistic
#' @param x X is the value used to find power
get_integrand_power = function(k,theta,v,z,x) {
  root = sqrt(v/2)
  product = 1
  for ( i in 2:k) {
    product = product*pnorm((x-theta1[i]*v)/root)
  }
  integrand = product*dnorm((x-theta1[1]*v)/root)*dnorm((x-z)/root)*2/v
  return(integrand)
}
