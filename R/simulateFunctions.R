#' Simulate Data For CAR Function
#'
#' Function to simulate data with treatment assignments from Covariate Adaptive Randomization Scheme
#' @param mean.s Mean of the short term response
#' @param sigma sigma in bivariate normal distribution
#' @param tau1 Tau1 is the weight on the first covariate
#' @param tau2 Tau2 is the weight on the second covariate
#' @param treat List containing the treatment assignments and the base matrix
#' @param covValues covValues is the dummy variable matrix from the covariate values
#' @param data Is the data from the previous analysis if applicable. The default is NULL
#' @export
simulatedata.car  <- function(mean.s,sigma = 1, tau1 = 1, tau2 = 1,treat,covValues,data=NULL,dist,inspection)
  # simulate multiple treatment groups
{

  trts = unique(treat)
  trts = trts[order(trts)]
  if (inspection == 1) keeps = rep(TRUE, nrow(covValues)) else keeps = data$treat %in% trts
  keeps = c(keeps,rep(TRUE,nrow(covValues)-length(keeps)))
  # if (inspection <= 1) covValues = covValues else covValues = covValues[(length(data$treat)+1):nrow(covValues),]
  if (inspection <= 1) covValues <<- covValues else covValues <<- covValues[keeps,]
  n.trt = length(unique(treat))-1
  if (inspection <= 2) n.s = table(treat) else n.s = table(treat)- table(data[data$treat %in% trts,]$treat) # n.s is the short term
  if (inspection == 2) data = NULL
  treat = tail(treat,sum(n.s))
  # if (dist == "normal") full.data = simnormal(mean.s = mean.s,sigma,trts,n.s,tau1,tau2,get_power = 0,covValues) else print(TRUE)
  print(n.s)
  if (dist == "normal") full.data = simnormal(mean.s = mean.s,sigma,trts,n.s,tau1,tau2,get_power = 0,covValues,treat) else full.data = simbin(mean.s,trts,n.s,treat)
  print(n.s)
  print(nrow(full.data))
  out = data.frame(s = full.data,treat = treat)
  out = rbind(data[data$treat %in% trts,],out)
  colnames(out) = c("s","treat")
  out
  }

#' Sim normal function
#'
#' Function to simulate normally distributed data based on given treatment groups
#' @param mean.s Mean.s is the mean of the outcome varaible
#' @param sigma Sigma is the variance
#' @param trts Trts is the possible treatment assignments
#' @param n.s N.s is the frequency of each treatment
#' @param tau1 Tau1 is the weight assigned to the first covariate value
#' @param tau2 Tau2 is the weight assigned to the second covariate value
#' @param get_power Get_power determines whether or not the power is being calculated in this simulation
simnormal = function (mean.s,sigma,trts,n.s, tau1, tau2, get_power = 0,covValues,treat) {
  mean.full.s = rep(rep(0,length(trts)),n.s)

  for (trt in trts)
  {
    mean.full.s[treat==trt] = rep(mean.s[which(trts==trt)], n.s[which(trts==trt)])

  }
  mean.full = cbind(mean.full.s)
  mean = c(0)
  var = c(sigma)
  dim(var) = c(1,1)
  full.error = mvtnorm::rmvnorm(nrow(mean.full),mean,var)
  covmatrix = matrix(c(tau1*covValues[,1],tau2*covValues[,2]),ncol=2)
  powermatrix = matrix(c(get_power*replace(treat,treat!=0,1)),ncol=1)
  full.data = full.error + mean.full + tau1*covValues[,1] + tau2*covmatrix[,2] + powermatrix
  full.data
}
#' simbin function
#'
#' This function is used to simulate binary endpoint data
#' @param mean.s mean.s is the mean of the binary endpoint
#' @param trts Trts is the vector of possible treatments
#' @param n.s N.s is the frequency of each treatment
simbin = function(mean.s,trts,n.s,treat) {
  print(TRUE)
  mean.full.s = rep(rep(0,length(trts)),n.s)
  for (trt in trts )
  {
    print(n.s[which(trts==trt)])
    print(mean.s[which(trts==trt)])

    mean.full.s[treat == trt] = rbinom(n.s[which(trts==trt)],1,mean.s[which(trts==trt)])

  }
  mean.full = cbind(mean.full.s)
  full.data = data.frame(s = mean.full.s)
  full.data

}
