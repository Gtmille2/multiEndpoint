#' get.z.v.full function
#'
#' This function is used to calculate the Z and V test statistics
#' @param data Data is data containing all the interim analysis data. In the form of a long-term endpoint followed by a short-term endpoint column.
#' @param inspection The current inspection in the trial
#' @param dist Dist the distribution of the endpoint at this inspection
#' @export
get.z.v.full = function(full.data,inspection,dist = "normal"){


  if (dist == "normal") z.v = get.z.v.norm(full.data) else z.v = get.z.v.bin(full.data)
  trts = unique(full.data$treat)
  trts = trts[order(trts)]
  trts = trts[trts>0]
  z = z.v[,1]
  v = z.v[,2]
  effect = z/v
  z.v
  best <<- trts[which.max(effect)]
  monitorArray[[inspection +1]][ZI] <<- z[which(trts==best)]
  monitorArray[[inspection +1]][VI] <<- v[which(trts==best)]
  print(z[which(trts==best)])
  print(v[which(trts==best)])
}

#' get.z.v.norm function
#'
#' This function computes the Z and V test statistics for a normally distributed endpoint
#' @param data The data containing a noramlly distributed endpoint and treatment assignments
get.z.v.norm = function(full.data) {
  n.trt = length(unique(full.data$treat))-1
  mean.full = rep(0,N)
  z.t = rep(0,n.trt)
  v.t = rep(0,n.trt)
  theta.t = rep(0,n.trt)
  for ( trt in 1:max(n.trt)){
    n=(length(full.data[treat==trt,]$s)+length(full.data[treat==0,]$s))
    n.e = length(full.data[treat==trt,]$s)
    n.c = length(full.data[treat==0,]$s)
    x.bar.e = mean(full.data[treat==trt,]$s)
    x.bar.c = mean(full.data[treat==0,]$s)
    table(treat)

    #
    # z.bin[trt] = length(full.data[treat==trt,]$s)*length(full.data[treat==0,]$s)*(props[trt+1]-props[1])/n
    # v.bin[trt] = length(full.data[treat==trt,]$s)*length(full.data[treat==0,]$s)*(length(full.data[treat==trt,]$s)*props[trt+1]+length(full.data[treat==0,]$s)*props[1])*(length(full.data[treat==trt,]$s)*(1-props[trt+1])+length(full.data[treat==0,]$s)*(1-props[1]))/n^2
    #
    pooleds = sqrt(sd(full.data[treat==trt,]$s)^2 + sd(full.data[treat==0,]$s)^2)/2
    theta.t[trt] = (mean(full.data[treat==trt,]$s) - mean(full.data[treat==0,]$s))/pooleds
    Q = sum(full.data[treat==trt,]$s^2) + sum(full.data[treat==0,]$s^2)
    D = sqrt((Q/n) - ((n.e*x.bar.e+n.c*x.bar.c)/n)^2)
    z.t[trt] = n.e*n.c*(x.bar.e-x.bar.c)/(n*D)
    v.t[trt] = n.e*n.c/n - (z.t[trt]^2 /(2*n))

  }
  # effect = z.t/v.t
  matrix(c(z.t,v.t),ncol=2)
  }

#' get.z.v.bin
#'
#' This function calculates the Z and V test statistics for a binomially distributed endpoint
#' @param full.data Data containing the binomally distributed endpoint and the treatment assignments
get.z.v.bin = function(full.data) {
  n.trt = length(unique(full.data$treat)) - 1
  trts = unique(full.data$treat)
  trts = trts[order(trts)]
  props = c()
  for ( trt in trts) {
    prop = sum(full.data[treat == trt,]$s)/length(full.data[treat == trt,]$s)
    props = c(props,prop)
  }
  theta.bin = rep(0,n.trt)
  for ( i in 1:max(n.trt)) theta.bin[i] = log((props[i+1]*(1-props[1]))/((1-props[i+1])*props[1]))

  z.bin = rep(0,n.trt)
  v.bin = rep(0,n.trt)
  for ( trt in trts[trts>0]){
    n=(length(full.data[treat==trt,]$s)+length(full.data[treat==0,]$s))
    n.e = length(full.data[treat==trt,]$s)
    n.c = length(full.data[treat==0,]$s)
    s.e = sum(full.data[treat==trt,]$s)
    s.c = sum(full.data[treat==0,]$s)
    s = s.e + s.c
    f = n - s
    # p.e = mean.bin[trt+1]
    # p.c = mean.bin[1]
    # # p.e = props[trt+1] # Not sure if this is correct
    # p.c = props[1] # Not sure if this is correct
    # x.bar.e = mean(full.data[treat==trt,]$s)
    # x.bar.c = mean(full.data[treat==0,]$s)
    z.bin[which(trts==trt)-1] = ((n.c*s.e) - (n.e*s.c))/n
    v.bin[which(trts==trt)-1] = (n.e*n.c*s*f)/n^3
    # z.bin[trt] = length(full.data[treat==trt,]$s)*length(full.data[treat==0,]$s)*(props[trt+1]-props[1])/n
    # v.bin[trt] = length(full.data[treat==trt,]$s)*length(full.data[treat==0,]$s)*(length(full.data[treat==trt,]$s)*props[trt+1]+length(full.data[treat==0,]$s)*props[1])*(length(full.data[treat==trt,]$s)*(1-props[trt+1])+length(full.data[treat==0,]$s)*(1-props[1]))/n^2


  }
  matrix(c(z.bin,v.bin),ncol=2)

  }
