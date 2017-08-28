require(jagsUI)
# poisson
sink("ptune.jags")
cat("
    model {
    #omega ~ dunif(0, 1) # zero inflation parameter
    for(i in 1:nsites){
    p[i] <- PT
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

params.pc
params.pt <- params.pc[-2]

data.dma.c1.pt <- data.dma.c1

# .3
dma.pt.c1 <- autojags(data.dma.c1.pt, inits.dm, params.pt, "ptune.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dma.pt.c1)
mc1_bpv <- bpv(dma.pt.c1, "monticola model: betweenness")

# set to .24
data.dma.c1.pt$PT <- .2
dma.pt.c1.low <- autojags(data.dma.c1.pt, inits.dm, params.pt, "ptune.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dma.pt.c1.low)
mc1_bpv.low <- bpv(dma.pt.c1.low, "monticola model: betweenness, low p")

# set to .36
data.dma.c1.pt$PT <- .4
dma.pt.c1.hi <- autojags(data.dma.c1.pt, inits.dm, params.pt, "ptune.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dma.pt.c1.hi)
mc1_bpv.hi <- bpv(dma.pt.c1.hi, "monticola model: betweenness, high p")

# boxplots
par(mfrow=c(1,3))
boxplot(dma.pt.c1.hi$sims.list$beta6, main="high p", ylim=c(-.8,.8))
boxplot(dma.pt.c1$sims.list$beta6, main="normal p", ylim=c(-.8,.8))
boxplot(dma.pt.c1.low$sims.list$beta6, main="low p", ylim=c(-.8,.8))


# 7:1, fuscus
dfa.pt.c7 <- autojags(data.dfa.c7, inits.df, params.pt, "ptune.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dfa.pt.c7)
mf7c_bpv <- bpv(dfa.pt.c7, "fuscus model: closeness (1:7)")

dfa.pt.c2 <- autojags(data.dfa.c2, inits.df, params.pt, "ptune.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dfa.pt.c2)
mf2c_bpv <- bpv(dfa.pt.c2, "fuscus model: closeness (1:2)")


# for WAIC calculation, see:
# https://jfiksel.github.io/2017-05-24-waic_aft_models_jags/
