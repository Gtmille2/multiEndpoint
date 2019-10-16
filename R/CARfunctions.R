#' Function for getting the treatment assignments
#'
#' Get treamtments points either to Pocock's design or Stratified Permuted Block Design
#' @param covValues Covariate Values for the number of patients in the trial
#' @param best Best is the selected best treatment after the first analysis
#' @param Tr Treatment assignments passed back to the treatment assignment function after the first interim analysis
#' @param Design Design is the variable determing whether there to use Pocock's design or Stratified Permuted Block Randomization
#' @export N.trt Number of treatments in the trial. After the first interim analysis, there will only be the control and the selected treatment
getTreat = function(covValues, best = 0, tr = NULL, n.trt, design) {
  if ( design == "Pocock" ) treat = psd(covValues, best = best, tr = tr,n.trt = n.trt) else treat = spbd(covValues,best = best, tr = tr, n.trt = n.trt)
  treat
  }



#' Function for returning treatment assignments in Pocock's Design
#'
#' This function takes a vector of factor level assignments, and returns treatment assignments
#' @param x Vector of factor level assignments
#' @param p Default is 3/4.
#' @param n.trt Number of treatments used in the trial
#' @param Best is the selected best treatment if past the first analysis. Default is 0
#' @param tr Tr is the vector of treamtnet assignments. Default is NULL
#' @param n.trt The number of treatments used in the analaysis if at the first stage in the analysis.
#' @examples
#' psd(x = seamlessTrials::all)
psd=function(x,p1=3/4,best = 0,tr = NULL,n.trt)
{
  if (best == 0) trts = seq(0,n.trt) else trts = c(0,best)
  if (is.null(tr)) keeps = rep(TRUE,nrow(x)) else keeps = c(tr %in% trts,rep(TRUE,nrow(x)-length(tr)))
  if (is.null(tr)) arrival = as.matrix(ade4::acm.disjonctif(x)) else arrival = as.matrix(ade4::acm.disjonctif(x[keeps,]))
  # arrival = as.matrix(ade4::acm.disjonctif(x))
  weight = rep(1,ncol(arrival))
  s = length(tr[tr %in% trts]) + 1
  n.trt = length(trts)
  N = nrow(arrival)
  base = matrix(rep(0,n.trt*ncol(arrival)),ncol=ncol(arrival))
  for ( i in 1:length(trts)) base[i,] = colSums(arrival[which(tr[tr %in% trts]==trts[i]),])
  if (is.null(tr)) tr=rep(NA,N) else tr = tr[keeps]
  out = list()
  # trin = which(trts==trts)
  delts = rep(0,(n.trt))
  while(s<=N)
  {
    for ( i in 1:length(trts)) {

      inc = base
      inc[i,] = base[i,] + t(arrival[s,])
      delts[i]=t(arrival[s,])%*%(Rfast::colrange(inc)*weight)
    }
    p=g(delts,p1,n.trt)
    tr[s] = sample(trts,1,p,replace=TRUE)
    base[which(trts==tr[s]),] = base[which(trts==tr[s]),] + t(arrival[s,])
    s=s+1
  }
  tr
  }

#' A function to calculate G in the Pocock Design and return a vector of p
#'
#' This function is not generally used. It is used in the calculation of the treatment assignments in psd
#' @param x x is the vector of deltas in Pocock's design
#' @param p1 = 3/4
g=function(x,p1=3/4,n.trt)
{
  p = NULL
  best = which(x == min(x))
  if (length(best) > 1) {
    p = rep(1/n.trt,n.trt)
  } else {
    p = rep(1-p1/(n.trt-1),n.trt)
    p[best] = 3/4

  }
  p

}

#' spbd Function
#'
#' Function to determine treatment assignments in SPB Design
#' @param x X vector of covariate values, ranging from 0 to total number of factor levels
#' @param n Total sample size
#' @param m M is the level of discretization
#' @param n.trt n.trt is the number of treatments in the trial
#' @export
#' @examples
#' spbd()
spbd=function(covValues,m=4, best = 0, tr = NULL, n.trt)
{
  if (best == 0) trts = seq(0,n.trt) else trts = c(0,best)
  if (is.null(tr)) keeps = rep(TRUE,nrow(covValues)) else keeps = c(tr %in% trts,rep(TRUE,nrow(covValues)-length(tr)))
  n.trt = length(trts)
  # if (is.null(tr)) arrival = as.matrix(ade4::acm.disjonctif(x)) else arrival = as.matrix(ade4::acm.disjonctif(x[keeps,]))
  s = length(tr) + 1
  z1 = covValues[,1]-1
  z2 = covValues[,2]-1
  n = length(z1)
  x=NULL  #covaiate info used for function spbd, values from 1 to 4
  #think I want this to be 100, 130-30
  while(s<=n)
  {
    if(z1[s]==1 & z2[s]==1) x[s]=1 else
      if(z1[s]==1 & z2[s]==0) x[s]=2 else
        if(z1[s]==0 & z2[s]==1) x[s]=3 else
          if(z1[s]==0 & z2[s]==0) x[s]=4

          s=s+1
  }
  x = na.omit(x)
  trkeeps = tr[tr %in% trts]
  i = 1
  tr=NULL
  while(i<=m)
  {
    # print(i)
    # print(length(x[x==i]))

    tr[x==i]=pbr(length(x[x==i]),n.trt = n.trt,trts=trts)
    # print(tr)
    i=i+1
  }
  tr = c(trkeeps,tr)
  return(tr)
}

#' PBR Function
#'
#' This function is not generally used. This is used in spbd to determine treatment assignments
#' @param n Total n size
#' @param block.size This is the number of subjects in each block, needs to be a multiple of the number of treatments
#' @param n.trt is the number of treatments in the trial, including the control
pbr=function(n,block.size=12,n.trt,trts) #block.size is the number of subjects in each block, assuming fixed
{
  # trts = trts[order(trts)]
  block.num=ceiling(n/block.size)
  # print(block.num)
  cards=NULL
  i=1
  while(i<=block.num)
  {

    full = matrix(sort(rep(trts,block.size/2),decreasing = TRUE),ncol=n.trt)
    # print(full)
    cards = c(cards, sample(full,block.size))
    # cards=c(cards,sample(cbind(rep(1,block.size/2),rep(0,block.size/2)),block.size))
    i=i+1
    # print(cards)
  }
  # print(n)
  cards=cards[1:n] #Why up to n? # Is this the step that deals with "a block size larger than the sum of the treatment group ratio"

  return(cards)
}



