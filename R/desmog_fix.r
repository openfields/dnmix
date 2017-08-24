# desmog_fix.r
# models with fixed detection probability, .55
require(jagsUI)
# poisson
sink("pconst.jags")
cat("
    model {
    #omega ~ dunif(0, 1) # zero inflation parameter
    for(i in 1:nsites){
      p[i] ~ dunif(0.25, 0.35)
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

params.pc <- c("N", "p", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "fit", "fit.new")

dfa.pc.bc.3 <- autojags(data.dfa.bc, inits.df, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
mean(dfa.pc.bc.3$sims.list$fit.new>dfa.pc.bc.3$sims.list$fit)
fbc_bpv <- bpv(dfa.pc.bc.3, "Model: betweenness")
#[1] 0.4026667

plot(dfa.pc.bc.3$sims.list$fit.new, dfa.pc.bc.3$sims.list$fit, xlim=c(20,60), ylim=c(20,60))
abline(0,1)

dfa.pc.c1 <- autojags(data.dfa.c1, inits.df, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dfa.pc.c1)
fc1_bpv <- bpv(dfa.pc.c1, "Model: closeness 1:1")

dfa.pc.c2 <- autojags(data.dfa.c2, inits.df, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dfa.pc.c2)
fc2_bpv <- bpv(dfa.pc.c2, "Model: closeness 1:2")

dfa.pc.c3 <- autojags(data.dfa.c3, inits.df, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dfa.pc.c3)
fc3_bpv <- bpv(dfa.pc.c3, "Model: closeness 1:3")

dfa.pc.c4 <- autojags(data.dfa.c4, inits.df, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dfa.pc.c4)
fc4_bpv <- bpv(dfa.pc.c4, "Model: closeness 1:4")

dfa.pc.c5 <- autojags(data.dfa.c5, inits.df, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dfa.pc.c5)
fc5_bpv <- bpv(dfa.pc.c5, "Model: closeness 1:5")

dfa.pc.c6 <- autojags(data.dfa.c6, inits.df, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dfa.pc.c6)
fc6_bpv <- bpv(dfa.pc.c2, "Model: closeness 1:6")

dfa.pc.c7 <- autojags(data.dfa.c7, inits.df, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dfa.pc.c7)
fc7_bpv <- bpv(dfa.pc.c7, "Model: closeness 1:7")

dfa.pc.c8 <- autojags(data.dfa.c8, inits.df, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dfa.pc.c8)
fc8_bpv <- bpv(dfa.pc.c8, "Model: closeness 1:8")

dfa.pc.c9 <- autojags(data.dfa.c9, inits.df, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dfa.pc.c9)
fc9_bpv <- bpv(dfa.pc.c9, "Model: closeness 1:9")

dfa.pc.c10 <- autojags(data.dfa.c10, inits.df, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dfa.pc.c10)
fc10_bpv <- bpv(dfa.pc.c10, "Model: closeness 1:10")


###------------------------------------
### monticola
###------------------------------------
dma.pc.c1 <- autojags(data.dma.c1, inits.dm, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dma.pc.c1)
mc1_bpv <- bpv(dma.pc.c1, "monticola model: betweenness")

dma.pc.c1 <- autojags(data.dma.c1, inits.dm, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dma.pc.c1)
mc1_bpv <- bpv(dma.pc.c1, "monticola model: closeness 1:1")

dma.pc.c2 <- autojags(data.dma.c2, inits.dm, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dma.pc.c2)
mc2_bpv <- bpv(dma.pc.c2, "monticola model: closeness 1:2")

dma.pc.c3 <- autojags(data.dma.c3, inits.dm, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dma.pc.c3)
mc3_bpv <- bpv(dma.pc.c3, "monticola model: closeness 1:3")

dma.pc.c4 <- autojags(data.dma.c4, inits.dm, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dma.pc.c4)
mc4_bpv <- bpv(dma.pc.c4, "monticola model: closeness 1:4")

dma.pc.c5 <- autojags(data.dma.c5, inits.dm, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dma.pc.c5)
mc5_bpv <- bpv(dma.pc.c5, "monticola model: closeness 1:5")

dma.pc.c6 <- autojags(data.dma.c6, inits.dm, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dma.pc.c6)
mc6_bpv <- bpv(dma.pc.c6, "monticola model: closeness 1:6")

dma.pc.c7 <- autojags(data.dma.c7, inits.dm, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dma.pc.c7)
mc7_bpv <- bpv(dma.pc.c7, "monticola model: closeness 1:7")

dma.pc.c8 <- autojags(data.dma.c8, inits.dm, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dma.pc.c8)
mc8_bpv <- bpv(dma.pc.c8, "monticola model: closeness 1:8")

dma.pc.c9 <- autojags(data.dma.c9, inits.dm, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dma.pc.c9)
mc9_bpv <- bpv(dma.pc.c9, "monticola model: closeness 1:9")

dma.pc.c10 <- autojags(data.dma.c10, inits.dm, params.pc, "pconst.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dma.pc.c10)
mc10_bpv <- bpv(dma.pc.c10, "monticola model: closeness 1:10")
