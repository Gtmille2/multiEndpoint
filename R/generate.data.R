#'Generate Data Function
#'
#'This function is used to generate data for use in a simulation of the seamless trails with different endpoints
#'@param n1 n1 is the number of people recruited overall in the study
#'@param mean.bin Mean.bin is the mean of the binary endpoint
#'@param mean.t Mean.t is the mean of the normally distributed endpoint
#'@param sigma Sigma is the standard deviation of the endpoints, which is assumed to be equal
#'@param rho Rho is the correlation between the endpoints
generate.data = function(n1,mean.bin,mean.t,sigma, rho) {
var.type = list("BinomDist", "NormalDist")
# Outcome distribution parameters
placebo.par = parameters(parameters(prop = mean.bin[1]),
                         parameters(mean = mean.t[1], sd = sigma))

dosel.par1 = parameters(parameters(prop = mean.bin[2]),
                        parameters(mean = mean.t[2], sd = sigma))

doseh.par1 = parameters(parameters(prop = mean.bin[3]),
                        parameters(mean = mean.t[3], sd = sigma))

# Correlation between two endpoints
corr.matrix = matrix(c(sigma, rho*sigma*sigma,
                       rho*sigma*sigma, sigma), 2, 2)

# Outcome parameter set 1
outcome1.placebo = parameters(type = var.type,
                              par = placebo.par,
                              corr = corr.matrix)
outcome1.dosel = parameters(type = var.type,
                            par = dosel.par1,
                            corr = corr.matrix)
outcome1.doseh = parameters(type = var.type,
                            par = doseh.par1,
                            corr = corr.matrix)

# Data model
case.study5.data.model = DataModel() +
  OutcomeDist(outcome.dist = "MVMixedDist") +
  SampleSize(c(n1)) +
  Sample(id = list("Placebo Outcome 1", "Placebo Outcome 2"),
         outcome.par = parameters(outcome1.placebo)) +
  Sample(id = list("Treatment 1 Outcome 1", "Treatment 1 Outcome 2"),
         outcome.par = parameters(outcome1.dosel)) +
  Sample(id = list("Treatment 2 Outcome 1", "Treatment 2 Outcome 1"),
         outcome.par = parameters(outcome1.doseh))

# Simulation Parameters
case.study5.sim.parameters =  SimParameters(n.sims = numberInspectionsMax,
                                            proc.load = "full",
                                            seed = 42938001)

case.study5.data.stack = DataStack(data.model = case.study5.data.model,
                                   sim.parameters = case.study5.sim.parameters)

full.data = list()
for (j in 1:numberInspectionsMax) {
  # full.data.bin = NULL
  # full.data.t = NULL
  full.data.both = NULL
  for (i in 1:6) {
    full.data.both = cbind(full.data.both, case.study5.data.stack$data.set[[j]]$data.scenario[[1]]$sample[[i]]$data$outcome)
    # full.data.t = cbind(full.data.t,case.study5.data.stack$data.set[[100]]$data.scenario[[1]]$sample[[i+1]]$data$outcome)
    #
  }
  full.data[[j]] = full.data.both
}
return(full.data)
}
