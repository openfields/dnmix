library(jagsUI)

# poisson
sink("pois1.jags")
cat("
    model {
    #omega ~ dunif(0, 1) # zero inflation parameter
    beta0 ~ dnorm(0,.1) # intercept
    beta1 ~ dnorm(0,.1) # pc1
    beta2 ~ dnorm(0,.1) # pc2
    beta3 ~ dnorm(0,.1) # pc3
    beta4 ~ dnorm(0,.1) # RT
    beta5 ~ dnorm(0,.1) # T1
    beta6 ~ dnorm(0,.1) # centrality parameter 1

    b4 ~ dnorm(0,.1) # pc1 - detection prob
    b5 ~ dnorm(0,.1) # pc2 - detection prob
    p0~dnorm(0,.1) # int - detection prob

    for(i in 1:nsites){
      p[i]<- p0 + b4*cov1[i] + b5*cov2[i]# could have covariates here
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

dfa.pois.bc <- autojags(data.dfa.bc, inits.df, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dfa.pois.c1 <- autojags(data.dfa.c1, inits.df, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dfa.pois.c2 <- autojags(data.dfa.c2, inits.df, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dfa.pois.c3 <- autojags(data.dfa.c3, inits.df, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dfa.pois.c4 <- autojags(data.dfa.c4, inits.df, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dfa.pois.c5 <- autojags(data.dfa.c5, inits.df, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dfa.pois.c6 <- autojags(data.dfa.c6, inits.df, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dfa.pois.c7 <- autojags(data.dfa.c7, inits.df, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dfa.pois.c8 <- autojags(data.dfa.c8, inits.df, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dfa.pois.c9 <- autojags(data.dfa.c9, inits.df, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dfa.pois.c10 <- autojags(data.dfa.c10, inits.df, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)

dma.pois.bc <- autojags(data.dma.bc, inits.dm, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dma.pois.c1 <- autojags(data.dma.c1, inits.dm, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dma.pois.c2 <- autojags(data.dma.c2, inits.dm, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dma.pois.c3 <- autojags(data.dma.c3, inits.dm, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dma.pois.c4 <- autojags(data.dma.c4, inits.dm, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dma.pois.c5 <- autojags(data.dma.c5, inits.dm, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dma.pois.c6 <- autojags(data.dma.c6, inits.dm, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dma.pois.c7 <- autojags(data.dma.c7, inits.dm, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dma.pois.c8 <- autojags(data.dma.c8, inits.dm, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dma.pois.c9 <- autojags(data.dma.c9, inits.dm, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dma.pois.c10 <- autojags(data.dma.c10, inits.dm, parameters, "pois1.jags", n.chains=nc, n.thin=2, parallel=TRUE)

mean(dfa.pois.bc$sims.list$fit.new>dfa.pois.bc$sims.list$fit) ->bp.f.b
mean(dfa.pois.c1$sims.list$fit.new>dfa.pois.c1$sims.list$fit) -> bp.f.c1
mean(dfa.pois.c2$sims.list$fit.new>dfa.pois.c2$sims.list$fit) -> bp.f.c2
mean(dfa.pois.c3$sims.list$fit.new>dfa.pois.c3$sims.list$fit) -> bp.f.c3
mean(dfa.pois.c4$sims.list$fit.new>dfa.pois.c4$sims.list$fit) -> bp.f.c4
mean(dfa.pois.c5$sims.list$fit.new>dfa.pois.c5$sims.list$fit) -> bp.f.c5
mean(dfa.pois.c6$sims.list$fit.new>dfa.pois.c6$sims.list$fit) -> bp.f.c6
mean(dfa.pois.c7$sims.list$fit.new>dfa.pois.c7$sims.list$fit) -> bp.f.c7
mean(dfa.pois.c8$sims.list$fit.new>dfa.pois.c8$sims.list$fit) -> bp.f.c8
mean(dfa.pois.c9$sims.list$fit.new>dfa.pois.c9$sims.list$fit) -> bp.f.c9
mean(dfa.pois.c10$sims.list$fit.new>dfa.pois.c10$sims.list$fit) -> bp.f.c10
bp.f <- c(bp.f.b, bp.f.c1, bp.f.c2, bp.f.c3, bp.f.c4, bp.f.c5, bp.f.c6, bp.f.c7, bp.f.c8, bp.f.c9, bp.f.c10)

mean(dma.pois.bc$sims.list$fit.new>dma.pois.bc$sims.list$fit) ->bp.m.b
mean(dma.pois.c1$sims.list$fit.new>dma.pois.c1$sims.list$fit) -> bp.m.c1
mean(dma.pois.c2$sims.list$fit.new>dma.pois.c2$sims.list$fit) -> bp.m.c2
mean(dma.pois.c3$sims.list$fit.new>dma.pois.c3$sims.list$fit) -> bp.m.c3
mean(dma.pois.c4$sims.list$fit.new>dma.pois.c4$sims.list$fit) -> bp.m.c4
mean(dma.pois.c5$sims.list$fit.new>dma.pois.c5$sims.list$fit) -> bp.m.c5
mean(dma.pois.c6$sims.list$fit.new>dma.pois.c6$sims.list$fit) -> bp.m.c6
mean(dma.pois.c7$sims.list$fit.new>dma.pois.c7$sims.list$fit) -> bp.m.c7
mean(dma.pois.c8$sims.list$fit.new>dma.pois.c8$sims.list$fit) -> bp.m.c8
mean(dma.pois.c9$sims.list$fit.new>dma.pois.c9$sims.list$fit) -> bp.m.c9
mean(dma.pois.c10$sims.list$fit.new>dma.pois.c10$sims.list$fit) -> bp.m.c10
bp.m <- c(bp.m.b, bp.m.c1, bp.m.c2, bp.m.c3, bp.m.c4, bp.m.c5, bp.m.c6, bp.m.c7, bp.m.c8, bp.m.c9, bp.m.c10)
