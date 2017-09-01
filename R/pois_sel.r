# poisson
sink("pois_sel.jags")
cat("
    model {
    for(i in 1:nsites){
      p[i] ~ dunif(0.2,0.5)
    }
    beta0 ~ dnorm(0,.1) # intercept
    beta1 ~ dnorm(0,.1) # pc1
    beta2 ~ dnorm(0,.1) # pc2
    beta3 ~ dnorm(0,.1) # pc3
    beta4 ~ dnorm(0,.1) # RT
    beta5 ~ dnorm(0,.1) # T1
    beta6 ~ dnorm(0,.1) # closeness
    beta7 ~ dnorm(0,.1) # betweenness    

    w1 ~ dbern(0.5)
    w2 ~ dbern(0.5)
    w3 ~ dbern(0.5)
    w4 ~ dbern(0.5)
    w5 ~ dbern(0.5)
    w6 ~ dbern(0.5)
    w7 ~ dbern(0.5)

    for(i in 1:nsites){
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
    N[i] ~ dpois(lambda[i])
    log(lambda[i])<- beta0 + w1*beta1*cov1[i] + w2*beta2*cov2[i] + w3*beta3*cov3[i] + w4*beta4*cov4[i] + w5*beta5*cov5[i]+ w6*beta6*cov6[i] + w7*beta7*cov7[i]
    
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
