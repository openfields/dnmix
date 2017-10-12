sink("ptune_w.jags")
cat("
    model {
    #omega ~ dunif(0, 1) # zero inflation parameter
    for(i in 1:nsites){
      p[i] <- .3 # ~ dunif(0,.4)
    }
    beta0 ~ dnorm(0,.1) # intercept
    beta1 ~ dnorm(0,.1) # pc1
    beta2 ~ dnorm(0,.1) # pc2
    beta3 ~ dnorm(0,.1) # pc3
    beta4 ~ dnorm(0,.1) # RT
    beta5 ~ dnorm(0,.1) # T1
    beta6 ~ dnorm(0,.1) # closeness centrality
    beta7 ~ dnorm(0,.1) # betweenness centrality
    w1 ~ dbern(0.5)
    w2 ~ dbern(0.5)
    w3 ~ dbern(0.5)
    w4 ~ dbern(0.5)
    w5 ~ dbern(0.5)
    w6 ~ dbern(0.5)
    w7 ~ dbern(0.5)

    for(i in 1:nsites){
      eps[i] ~ dnorm(0, tau.lam)
    }
    tau.lam <- 1/(sd.lam*sd.lam)
    sd.lam ~ dunif(0, 3)



#p0 ~ dunif(0.5,0.6)
    #     b4 ~ dnorm(0,.1) # pc1 - detection prob
    #     b5 ~ dnorm(0,.1) # pc2 - detection prob
    #     p0~dnorm(0,.1) # int - detection prob

    for(i in 1:nsites){
    #p[i]<- p0 #p0 + b4*cov1[i] + b5*cov2[i]# could have covariates here
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
    log(lambda[i])<- eps[i] + beta0 + beta1*cov1[i]*w1 + beta2*cov2[i]*w2 + beta3*cov3[i]*w3 + beta4*cov4[i]*w4 + beta5*cov5[i]*w5 + beta6*cov6[i]*w6 + beta7*cov7[i]*w7

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

params.w <- c(params.pt, "w1", "w2", "w3", "w4", "w5", "w6", "w7")

dma.pt.c2.w <- autojags(data.dma.c2, inits.dm, params.w, "ptune_w.jags", n.chains=nc, n.thin=1, parallel=TRUE)

w1 <- dma.pt.c2.w$sims.list$w1
w2 <- dma.pt.c2.w$sims.list$w2
w3 <- dma.pt.c2.w$sims.list$w3
w4 <- dma.pt.c2.w$sims.list$w4
w5 <- dma.pt.c2.w$sims.list$w5
w6 <- dma.pt.c2.w$sims.list$w6
w7 <- dma.pt.c2.w$sims.list$w7
mod <- paste(w1,w2,w3,w4,w5,w6,w7, sep="")

table(mod)
table(mod)/length(dma.pt.c2.w$sims.list$w1)
