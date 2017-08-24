# function bpv: calculates bayesian p-value and creates plot
# requires 'fit.new' and 'fit' objects within sims.list object from model fit with jagsUI
# created by WRF, 24 Aug 2017

bpv <- function(mcmc, title){
  # calculate Bayesian p-value
  bpv <- mean(mcmc$sims.list$fit.new>mcmc$sims.list$fit)
  # determine min/max for plot
  f.lim <- c(min(mcmc$sims.list$fit.new), max(mcmc$sims.list$fit.new), min(mcmc$sims.list$fit), max(mcmc$sims.list$fit))
  p.lim <- c(min(f.lim), max(f.lim))
  # create plot & add 1:1 line
  tex <- paste(title, "\n", "Bayesian p-value:", " ", round(bpv,4))
  plot(mcmc$sims.list$fit.new~mcmc$sims.list$fit,xlim=p.lim,ylim=p.lim, main=tex)
  abline(0,1,col="red")
  # return Bayesian p-value
  return(bpv)

}
