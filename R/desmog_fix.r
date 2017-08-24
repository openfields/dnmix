# desmog_fix.r
# models with fixed detection probability, .55
require(jagsUI)
# poisson
sink("pconst.jags")
cat("
    model {
    #omega ~ dunif(0, 1) # zero inflation parameter
    for(i in 1:nsites){
      p[i] ~ dunif(0.5, 0.6)
    }
    beta0 ~ dnorm(0,.1) # intercept
    beta1 ~ dnorm(0,.1) # pc1
    beta2 ~ dnorm(0,.1) # pc2
    beta3 ~ dnorm(0,.1) # pc3
    beta4 ~ dnorm(0,.1) # RT
    beta5 ~ dnorm(0,.1) # T1
    beta6 ~ dnorm(0,.1) # centrality parameter 1
    #p0 ~ dunif(0.5,0.6)
#     b4 ~ dnorm(0,.1) # pc1 - detection prob
#     b5 ~ dnorm(0,.1) # pc2 - detection prob
#     p0~dnorm(0,.1) # int - detection prob

    for(i in 1:nsites){
      p[i]<- p0 #p0 + b4*cov1[i] + b5*cov2[i]# could have covariates here
      mu[i,1]<- p[i]
      mu[i,2]<- p[i]*(1-p[i])
      mu[i,3]<- p[i]*(1-p[i])*(1-p[i])
      pi0[i]<- 1 - mu[i,1]-mu[i,2]-mu[i,3]
      pcap[i]<-1-pi0[i]

      for(j in 1:3){
        muc[i,j]<-mu[i,j]/pcap[i]
      }

      # 1. model part 1: the conditional multinomial
      y[i,1:3] ~ dmulti(muc[i,1:3],ncap[i])

      # 2. model for the observed count of uniques
      ncap[i] ~ dbin(pcap[i],N[i])

      # 3. abundance model
      #z[i] ~ dbern(omega)
      N[i] ~ dpois(lambda[i])
      #lam.eff[i] <- z[i]*lambda[i]
      log(lambda[i])<- beta0 + beta1*cov1[i] + beta2*cov2[i] + beta3*cov3[i] + beta4*cov4[i] + beta5*cov5[i]+ beta6*cov6[i]

      # fit stats
      for (j in 1:3){
        eval[i,j] <- p[i]*N[i]
        E[i,j] <- pow((y[i,j] - eval[i,j]),2)/(eval[i,j]+0.5)

        y.new[i,j] ~ dbin(p[i], N[i])
        E.new[i,j] <- pow((y.new[i,j] - eval[i,j]),2)/(eval[i,j]+0.5)
        }

    }

    # sum model fit stats
    fit <- sum(E[,])
    fit.new <- sum(E.new[,])

    }
    ",fill=TRUE)
sink()

dfa.pc.bc <- autojags(data.dfa.bc, inits.df, parameters, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)

