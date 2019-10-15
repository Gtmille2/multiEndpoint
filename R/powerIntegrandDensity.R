#' getErrorRates
#'
#' This function gets the error rates for trial
#' @param get_power Get_power determines whether or not this is returning the power of the trial
#' @param numberInspections The number of inspections in the trial
#' @param print Print is an option to print the error rates or not
getErrorRates = function(get_power,numberInspections,print) {
  alpha_u = 0
  alpha_l = 0

  densityFunctionDiscreteArray <<- densityArray(get_power,numberInspections,INTEGR_POINTS)

  for ( inspection in 1:(numberInspections+1)) {
    alpha_u = alpha_u + probCrossUpperDiscrete(get_power,inspection)
    if ( get_power == 0) alpha_l = alpha_l + probCrossLowerDiscrete(get_power,inspection)
    if ( print == 1 & get_power == 0) cat("alpha_u: ", alpha_u, "alpha_l", alpha_l,"\n")
    if ( print == 1 & get_power == 1) cat("alpha_u: ",alpha_u,"\n")
  }
  if (print==1 & get_power == 0) cat("alpha_u: ",alpha_u, "alpha_l", alpha_l,"\n")
  if (print==1 & get_power == 1) cat("alpha_u: ",alpha_u,"\n")
  return(alpha_u)

}

#' probCrossUpperDiscrete
#'
#' This function is used to determine the probability of crossing the upper boundary
#' @param get_power get_power is the arugment to calculate the power in the trial
#' @param inspection Inspection is the current inspection in the analysis
probCrossUpperDiscrete = function(get_power,inspection) {
  upperBound = monitorArray[[inspection+1]][UPPER_BOUNDI]
  vStep = monitorArray[[inspection+1]][VI] - monitorArray[[inspection]][VI]
  if(inspection == 1) {
    return(0)
  }
  else
  {
    return(integrateNormalDistributionDensityFunctionDiscrete(inspection))
  }
}

#' integrateNormalDistributionDensityFunctionDiscrete
#'
#' This function is used to compute the integral in the probCrossUpperDiscrete function
#' @param inspection Inspection is the current inspection in the trial
integrateNormalDistributionDensityFunctionDiscrete = function(inspection) {
  integral = 0
  vStep = monitorArray[[inspection+1]][VI] - monitorArray[[inspection]][VI]
  cat("vStep in integrate function: ",vStep,"inspection: ",inspection, "\n")
  for ( i in 0:(INTEGR_POINTS-1)) {
    s = densityFunctionDiscreteArray[[inspection-1]][ARRVALUE,i+1]
    if (inspection == 2) {
      mu_second = theta[1]*monitorArray[[3]][VI] - rho*sqrt(monitorArray[[3]][VI]/monitorArray[[2]][VI])*
        (theta[1]*monitorArray[[2]][VI]-s)
      var_second = monitorArray[[3]][VI]*(1-rho*rho)
      normalDistribution = 1 - pnorm((monitorArray[[inspection+1]][UPPER_BOUNDI] - mu_second)/sqrt(var_second))
    }
    else {
      normalDistribution = 1 - pnorm((monitorArray[[inspection+1]][UPPER_BOUNDI] - s - theta[1]*vStep)/sqrt(vStep))
    }
    #cat("Test 1: ",densityFunctionDiscreteArray[[inspection-1]][ARRWEIGHT,i + 1]," Test 2: ",
    #densityFunctionDiscreteArray[[inspection -1]][ARRDENSITY,i+1],"\n")
    integral = integral + densityFunctionDiscreteArray[[inspection-1]][ARRWEIGHT,i + 1]*
      densityFunctionDiscreteArray[[inspection -1]][ARRDENSITY,i+1]*
      normalDistribution
  }
  #cat("Integral result Here: ", integral, "\n")

  return(integral)
}

#' probCrossLowerDiscrete
#'
#' This function is used to compute the probability of crossing the lower boundary
#' @param get_power Get_power determines the power in the trial
#' @param inspection Inpsection is current inspection in the trial
probCrossLowerDiscrete = function(get_power, inspection) {
  lowerBound = monitorArray[[inspection+1]][LOWER_BOUNDI]
  vStep = monitorArray[[inspection+1]][VI] - monitorArray[[inspection]][VI]

  if(inspection == 1) return(0)
  else return (integrateNormalDistributionDensityFunctionDiscreteLower(inspection))
}


#' integrateNormalDistributionDensityFunctionDiscreteLower
#'
#' This function is used to compute the integral in the probCrossLowerDiscrete function
#' @param inspection Inspection is the current inspection in the trial
integrateNormalDistributionDensityFunctionDiscreteLower = function(inspection) {
  integral = 0
  vStep = monitorArray[[inspection+1]][VI] - monitorArray[[inspection]][VI]
  cat("vStep in integrate function: ",vStep,"inspection: ",inspection, "\n")

  for ( i in 0:(INTEGR_POINTS - 1)) {
    s = densityFunctionDiscreteArray[[inspection-1]][ARRVALUE,i+1]
    if ( inspection == 2) {
      mu_second = theta[1]*monitorArray[[3]][VI] - rho*sqrt(monitorArray[[3]][VI]/monitorArray[[2]][VI])*
        (theta[1]*monitorArray[[2]][VI]-s)
      var_second = monitorArray[[3]][VI]*(1-rho*rho)
      # #cat("mu_second: ",mu_second,"var_second: ",var_second,"s: ",s,"test: ",monitorArray[[inspection+1]][UPPER_BOUNDI] ,"\n")

      normalDistribution = pnorm((monitorArray[[inspection+1]][LOWER_BOUNDI] - mu_second)/sqrt(var_second))
      # #cat("normalDistribution: ",normalDistribution,"\n")

    }
    else {
      normalDistribution = pnorm((monitorArray[[inspection+1]][LOWER_BOUNDI] - s - theta[1]*vStep)/sqrt(vStep))
    }
    integral = integral + densityFunctionDiscreteArray[[inspection-1]][ARRWEIGHT,i+1]*
      densityFunctionDiscreteArray[[inspection-1]][ARRDENSITY,i+1]*
      normalDistribution

  }
  #cat("Integral result Here: ", integral, "\n")

  return(integral)
}

#' densityArray
#'
#' This function is used to compute the densities for each analysis in the trial in order to recursively find the stopping boundaries
#' @param get_power Get_power is an argument to determine the power of the trial
#' @param numberInspections Number of inspections in the current analysis
densityArray = function(get_power,numberInspections) {
  array = list(c(0))
  for ( i in 1:numberInspections) array[[i]] = matrix(rep(0,3*INTEGR_POINTS),nrow = 3)
  for ( i in 0:(numberInspections-1)) {

    lowerLimit = monitorArray[[i+2]][LOWER_BOUNDI]
    upperLimit = monitorArray[[i+2]][UPPER_BOUNDI]
    if (i == 0) lowerLimit = -10*sqrt(monitorArray[[2]][VI])
    if (i == 0) upperLimit = 10*sqrt(monitorArray[[2]][VI])
    cat("lowerLimit: ",lowerLimit,"upperLimit: ",upperLimit,"\n")

    arrayvalweights = gaussLegendre(n = INTEGR_POINTS,lowerLimit,upperLimit)

    array[[i+1]][ARRVALUE,] = arrayvalweights$x
    array[[i+1]][ARRWEIGHT,] = arrayvalweights$w

    vStep = monitorArray[[i+2]][VI] - monitorArray[[i+1]][VI]
    cat("vStep: ",vStep,"\n")
    if ( i ==0 ) {
      for( j in 0:(INTEGR_POINTS-1)){
        # z = array[[i+1]][ARRVALUE,j+1]
        z = monitorArray[[i+1]][ZI]

        # cat("z: ",z,"\n")
        if(get_power==0) {
          array[[i+1]][ARRDENSITY,j+1] = get_density_first(ntmt,theta,vStep,z)
        }
        if (get_power == 1) {
          array[[i+1]][ARRDENSITY,j+1] = get_power_first(ntmt,theta,vStep,z)

        }
      }
      # for ( r in 1:11 ) cat("BLAH HERE: ",array[[i+1]][ARRDENSITY,r],"\n")
    }
    if ( i ==1 ) {
      for ( j in 0:(INTEGR_POINTS-1)) {
        density = 0.0
        z = array[[i+1]][ARRVALUE,j+1]

        for( k  in 0:(INTEGR_POINTS-1)){
          s  = array[[i]][ARRVALUE,k+1]
          mu_second = theta[1]*monitorArray[[3]][VI] - rho*sqrt(monitorArray[[3]][VI]/monitorArray[[2]][VI])*
            (theta[1]*monitorArray[[2]][VI] - s)
          var_second = monitorArray[[3]][VI]*(1-rho*rho)
          density = density + array[[i]][ARRWEIGHT,k+1]*
            (1/sqrt(var_second))*dnorm((z-mu_second)/sqrt(var_second))*(array[[i]][ARRDENSITY,k+1])

        }
        #cat("density :",density,"i: ",i,"\n")

        array[[i+1]][ARRDENSITY,j+1] = density
      }
    }
    if ( i >= 2) {
      #print(TRUE)
      for(j in 0:(INTEGR_POINTS-1) )
      {
        density = 0
        z = array[[i+1]][ARRVALUE,j+1]

        for (k in 0:(INTEGR_POINTS-1)) {
          s = array[[i]][ARRVALUE,k+1]

          density = density + array[[i]][ARRWEIGHT,k+1]*(1/sqrt(vStep))*dnorm((z-s-theta[1]*vStep)/sqrt(vStep))*array[[i]][ARRDENSITY,k+1]
        }
        # cat("density :",density,"\n")

        array[[i+1]][ARRDENSITY,j+1] = density
      }

    }
  }
  return(array)
}
