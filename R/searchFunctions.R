#' V Search function
#'
#' Function to search for solution to V distribution
#' @param vmax Vmax is the max value passed to this function
v_search_fn = function(vmax) {
monitorArray[[1]][VI] <<- 0
# monitorArray[[2]][VI] <<- 12.5
print(vmax)
for ( j in 3:(numberInspectionsMax+1)) {
  monitorArray[[j]][VI] <<- vmax*(j-2)/(numberInspectionsMax-1)
  cat("Monitor Array at Point: ",j," ",monitorArray[[j]][VI]," \n")
}
for ( j in 1:ntmt) theta[j] <<- 0
for ( j in 1:ntmt) theta1[j] <<- 0
for ( numberInspections in 1:(numberInspectionsMax) ) {
  print(numberInspections)
  if(numberInspections == 1) {
    monitorArray[[numberInspections]][UPPER_BOUNDI]  <<- 0
    monitorArray[[numberInspections]][LOWER_BOUNDI] <<- 0
  }
  print(TRUE)
  if ( numberInspections > 1) {
    #cat("numberInspections here: ",numberInspections,"\n")
    densityFunctionDiscreteArray <<- densityArray(0,numberInspections-1,INTEGR_POINTS)

    upperBoundary = binarySearchUpper(0,5*sqrt(monitorArray[[numberInspections+1]][VI]),numberInspections)
    lowerBoundary = binarySearchLower(-5*sqrt(monitorArray[[numberInspections+1]][VI]),5*sqrt(monitorArray[[numberInspections+1]][VI]),numberInspections)
    #cat("lowerBoundary here: ", lowerBoundary,"\n")

    #Think I have to use two different functions for upper and lower since R can't do pointer functions
    monitorArray[[numberInspections+1]][UPPER_BOUNDI] <<- upperBoundary
    monitorArray[[numberInspections+1]][LOWER_BOUNDI] <<- lowerBoundary
    cat(monitorArray[[numberInspections + 1]][ZI],"upperBoundary: ",upperBoundary," lowerBoundary: ",lowerBoundary,"\n")
    if ( monitorArray[[numberInspections + 1]][ZI]  > upperBoundary || monitorArray[[numberInspections + 1]][ZI] < lowerBoundary ) {
      print("BLAHHHH")
      break
    }
  }
}
theta[1]<<-0.315
theta[2]<<-0.24
theta[3]<<-0.24
theta1[1]<<-(1.08-0.52)/.6
theta1[2]<<-(0.68-0.52)/.6
theta1[3]<<-(0.68-0.52)/.6
difference = getErrorRates(1,numberInspectionsMax-1,0)-.9
if (numberInspections > 1)
  cat("vmax: ",vmax, "difference: ",difference,"\n")
return(difference)
}

#'Upper search function
#'
#'This is a function to find the upper boundary
#'@param upperBoundary The upperBoundary passed to upper_search_fn
#'@param numberInspections The current number of numberInspections
upper_search_fn = function(upperBoundary,numberInspections) {
  target = alpha_star_u[numberInspections+1] - alpha_star_u[numberInspections]
  monitorArray[[numberInspections+1]][UPPER_BOUNDI] <<- upperBoundary
  monitorArray[[numberInspections+1]][LOWER_BOUNDI] <<- 0
  cat("target: ",target,"UpperBoundary: ",monitorArray[[numberInspections+1]][[UPPER_BOUNDI]],"\n")
  return(probCrossUpperDiscrete(0,numberInspections)-target)
}

#'Lower search function
#'
#'This is a function to find the upper boundary
#'@param lowerBoundary The upperBoundary passed to upper_search_fn
#'@param numberInspections The current number of numberInspections
lower_search_fn = function(lowerBoundary,numberInspections) {
  target = alpha_star_l[numberInspections + 1] - alpha_star_l[numberInspections]
  monitorArray[[numberInspections+1]][UPPER_BOUNDI] <<- 0
  monitorArray[[numberInspections+1]][LOWER_BOUNDI] <<- lowerBoundary
  return(probCrossLowerDiscrete(0,numberInspections) - target)

}


#' binarySearchUpper
#'
#' A general search function to find the root for the upper search
#' @param lowerLimit The lowerlimit of the search passed to the function
#' @param upperLimit The upperLimit of the search passed to the function
#' @param numberInspections The current number of inspection
binarySearchUpper = function(lowerLimit,upperLimit,numberInspections) {
  loop = 0
  f_upr = upper_search_fn(upperLimit,numberInspections)
  #cat("Number inspections: ", numberInspections, " f_upr: ",f_upr, "\n")
  while ( loop < BINARY_SEARCH_MAX_ITERATIONS) {
    mid_val = (upperLimit+lowerLimit)/2
    if( abs(mid_val) >= CLOSE_ZERO & abs((upperLimit - mid_val)/mid_val) < SEARCH_TOLERANCE)
      return(mid_val)

    loop = loop + 1
    f_mid = upper_search_fn(mid_val,numberInspections)
    #cat("f_mid: ", f_mid,"\n")
    if(f_mid*f_upr == 0 ){

      return(mid_val)
    }
    else{
      if(f_mid*f_upr > 0 )
      {
        upperLimit = mid_val
        f_upr = f_mid
      }
      else {
        lowerLimit = mid_val
      }
    }

  }
  return(mid_val)
}

#' binarySearchLower
#'
#' A general search function to find the root for the upper search
#' @param lowerLimit The lowerlimit of the search passed to the function
#' @param upperLimit The upperLimit of the search passed to the function
#' @param numberInspections The current number of inspection
binarySearchLower = function(lowerLimit,upperLimit,numberInspections) {
  loop = 0
  f_upr = lower_search_fn(upperLimit,numberInspections)

  while ( loop < BINARY_SEARCH_MAX_ITERATIONS) {
    mid_val = (upperLimit+lowerLimit)/2
    if( abs(mid_val) >= CLOSE_ZERO & abs((upperLimit - mid_val)/mid_val) < SEARCH_TOLERANCE)
      return(mid_val)

    loop = loop + 1
    f_mid = lower_search_fn(mid_val,numberInspections)
    if(f_mid*f_upr == 0 ){

      return(mid_val)
    }
    else{
      if(f_mid*f_upr > 0 )
      {
        upperLimit = mid_val
        f_upr = f_mid
      }
      else {
        lowerLimit = mid_val
      }
    }

  }
  #print(mid_val)
  return(mid_val)

}

#' binarySearchV
#'
#' THe binary search v function is a bsic search function to find the root for V
#' @param lowerLimit The lowerlimit of the search passed to the function
#' @param upperLimit The upperLimit of the search passed to the function
#' @param numberInspections The current number of inspection
binarySearchV = function(lowerLimit,upperLimit,numberInspections) {
  loop = 0
  f_upr = v_search_fn(upperLimit)
  #print(lowerLimit)
  #print(upperLimit)
  while ( loop < BINARY_SEARCH_MAX_ITERATIONS) {
    mid_val = (upperLimit+lowerLimit)/2
    if( abs(mid_val) >= CLOSE_ZERO & abs((upperLimit - mid_val)/mid_val) < SEARCH_TOLERANCE)
      return(mid_val)

    loop = loop + 1
    f_mid = v_search_fn(mid_val)
    if(f_mid*f_upr == 0 ){

      return(mid_val)
    }
    else{
      if(f_mid*f_upr > 0 )
      {
        upperLimit = mid_val
        f_upr = f_mid
      }
      else {
        lowerLimit = mid_val
      }
    }

  }
  return(mid_val)
}



