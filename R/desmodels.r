library(jagsUI)

# zero-inflated
sink("zi1.jags")
cat("
    model {
    omega ~ dunif(0, 1) # zero inflation parameter
    beta0 ~ dnorm(0,.1) # intercept
    beta1 ~ dnorm(0,.1) # pc1
    beta2 ~ dnorm(0,.1) # pc2
    beta3 ~ dnorm(0,.1) # pc3
    beta4 ~ dnorm(0,.1) # RT
    beta5 ~ dnorm(0,.1) # T1
    beta6 ~ dnorm(0,.1) # centrality parameter
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
    z[i] ~ dbern(omega)
    N[i] ~ dpois(lam.eff[i])
    lam.eff[i] <- z[i]*lambda[i]
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

# random site effects
sink("re1.jags")
cat("
    model {
    for(j in 1:nsites){
    alpha[j] ~ dnorm(0, tau.alpha)
    }

    tau.alpha <- 1/(sd.alpha*sd.alpha)
    sd.alpha ~ dunif(0,5)

    #omega ~ dunif(0, 1) # zero inflation parameter
    beta0 ~ dnorm(0,.1) # intercept
    beta4 ~ dnorm(0,.1) # RT
    beta5 ~ dnorm(0,.1) # T1
    beta6 ~ dnorm(0,.1) # centrality parameter
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
    log(lambda[i])<- beta0 + beta4*cov4[i] + beta5*cov5[i]+ beta6*cov6[i] + alpha[i]

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


# load data
read.csv('./data/dfadma_merge.csv', header=TRUE) -> desdata
desdata[,10:12]-> y.dfus
desdata[,26:28]-> y.dmon
y.dmon[is.na(y.dmon)] <- 0

# -------------------------------------------------------------------------
# fit models for D. fuscus: local habitat + centrality + network structure
# -------------------------------------------------------------------------
nsites <- dim(y.dfus)[1]
ncap.df<-apply(y.dfus,1,sum)
ymax.df<-ncap.df


# initial values
inits.df <- function(){
  list (p0=runif(1),beta0=runif(1,-1,1),N=ymax.df+1,z=rep(1,53))
}

# parameters to monitor
parameters <- c("N","p0","beta0","beta1","beta2","beta3","beta4","beta5","beta6","omega","b4","b5","fit","fit.new")
parameters.re <- c("N","alpha","p0","beta0","beta4","beta5","beta6","b4","b5","fit","fit.new")

# mcmc settings
nthin<-3
nc<-3
nb<-10000
ni<-30000

t1 <- Sys.time()

# betweenness centrality:
data.dfa.bc <- list(y=y.dfus,nsites=nsites,ncap=ncap.df,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$bcns))
dfa.zi.bc <- autojags(data.dfa.bc, inits.df, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE)
dfa.re.bc <- autojags(data.dfa.bc, inits.df, parameters.re, "re1.jags", n.chains=nc, n.thin=2, parallel=TRUE)

# closeness centrality - 1:1 (upstream/downstream)
data.dfa.c1 <- list(y=y.dfus,nsites=nsites,ncap=ncap.df,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c1n))
dfa.zi.c1 <- autojags(data.dfa.c1, inits.df, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)
dfa.re.c1 <- autojags(data.dfa.c1, inits.df, parameters.re, "re1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)

# closeness centrality - 1:2 (upstream/downstream)
data.dfa.c2 <- list(y=y.dfus,nsites=nsites,ncap=ncap.df,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c2n))
dfa.zi.c2 <- autojags(data.dfa.c2, inits.df, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)
dfa.re.c2 <- autojags(data.dfa.c2, inits.df, parameters.re, "re1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)

# closeness centrality - 1:3 (upstream/downstream)
data.dfa.c3 <- list(y=y.dfus,nsites=nsites,ncap=ncap.df,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c3n))
dfa.zi.c3 <- autojags(data.dfa.c3, inits.df, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)
dfa.re.c3 <- autojags(data.dfa.c3, inits.df, parameters.re, "re1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)

# closeness centrality - 1:4 (upstream/downstream)
data.dfa.c4 <- list(y=y.dfus,nsites=nsites,ncap=ncap.df,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c4n))
dfa.zi.c4 <- autojags(data.dfa.c4, inits.df, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)
dfa.re.c4 <- autojags(data.dfa.c4, inits.df, parameters.re, "re1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)

# closeness centrality - 1:5 (upstream/downstream)
data.dfa.c5 <- list(y=y.dfus,nsites=nsites,ncap=ncap.df,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c5n))
dfa.zi.c5 <- autojags(data.dfa.c5, inits.df, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)
dfa.re.c5 <- autojags(data.dfa.c5, inits.df, parameters.re, "re1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)

# closeness centrality - 1:6 (upstream/downstream)
data.dfa.c6 <- list(y=y.dfus,nsites=nsites,ncap=ncap.df,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c6n))
dfa.zi.c6 <- autojags(data.dfa.c6, inits.df, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)
dfa.re.c6 <- autojags(data.dfa.c6, inits.df, parameters.re, "re1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)

# closeness centrality - 1:7 (upstream/downstream)
data.dfa.c7 <- list(y=y.dfus,nsites=nsites,ncap=ncap.df,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c7n))
dfa.zi.c7 <- autojags(data.dfa.c7, inits.df, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)
dfa.re.c7 <- autojags(data.dfa.c7, inits.df, parameters.re, "re1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)

# closeness centrality - 1:8 (upstream/downstream)
data.dfa.c8 <- list(y=y.dfus,nsites=nsites,ncap=ncap.df,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c8n))
dfa.zi.c8 <- autojags(data.dfa.c8, inits.df, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)
dfa.re.c8 <- autojags(data.dfa.c8, inits.df, parameters.re, "re1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)

# closeness centrality - 1:9 (upstream/downstream)
data.dfa.c9 <- list(y=y.dfus,nsites=nsites,ncap=ncap.df,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c9n))
dfa.zi.c9 <- autojags(data.dfa.c9, inits.df, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)
dfa.re.c9 <- autojags(data.dfa.c9, inits.df, parameters.re, "re1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)

# closeness centrality - 1:10 (upstream/downstream)
data.dfa.c10 <- list(y=y.dfus,nsites=nsites,ncap=ncap.df,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                     cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c10n))
dfa.zi.c10 <- autojags(data.dfa.c10, inits.df, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)
dfa.re.c10 <- autojags(data.dfa.c10, inits.df, parameters.re, "re1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1)


# ----------------------------------------------------------------------------
# fit models for D. monticola: local habitat + centrality + network structure
# ----------------------------------------------------------------------------
# initial values
ncap.dm<-apply(y.dmon,1,sum)
ymax.dm<-ncap.dm
inits.dmzi <- function(){
  list (p0=runif(1),beta0=runif(1,-1,1),N=ymax.dm+1,z=rep(1,53))
}
inits.dm <- function(){
  list(p0=runif(1),beta0=runif(1,-1,1),N=ymax.dm+1)
}


# betweenness centrality:
data.dma.bc <- list(y=y.dmon,nsites=nsites,ncap=ncap.dm,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$bcns))
dma.zi.bc <- autojags(data.dma.bc, inits.dmzi, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1, max.iter = 300000)
dma.re.bc <- autojags(data.dma.bc, inits.dmzi, parameters, "re1.jags", n.chains = 3, n.thin = 2, parallel = TRUE, Rhat.limit = 1.1, max.iter = 300000)

# closeness centrality - 1:1 (upstream/downstream)
data.dma.c1 <- list(y=y.dmon,nsites=nsites,ncap=ncap.dm,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c1n))
dma.zi.c1 <- autojags(data.dma.c1, inits.dmzi, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1, max.iter = 300000)
dma.re.c1 <- autojags(data.dma.c1, inits.dmzi, parameters, "re1.jags", n.chains = 3, n.thin = 2, parallel = TRUE, Rhat.limit = 1.1, max.iter = 300000)

# closeness centrality - 1:2 (upstream/downstream)
data.dma.c2 <- list(y=y.dmon,nsites=nsites,ncap=ncap.dm,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c2n))
dma.zi.c2 <- autojags(data.dma.c2, inits.dmzi, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1, max.iter = 300000)
dma.re.c2 <- autojags(data.dma.c2, inits.dmzi, parameters, "re1.jags", n.chains = 3, n.thin = 2, parallel = TRUE, Rhat.limit = 1.1, max.iter = 300000)

# closeness centrality - 1:3 (upstream/downstream)
data.dma.c3 <- list(y=y.dmon,nsites=nsites,ncap=ncap.dm,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c3n))
dma.zi.c3 <- autojags(data.dma.c3, inits.dmzi, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1, max.iter = 300000)
dma.re.c3 <- autojags(data.dma.c3, inits.dmzi, parameters, "re1.jags", n.chains = 3, n.thin = 2, parallel = TRUE, Rhat.limit = 1.1, max.iter = 300000)

# closeness centrality - 1:4 (upstream/downstream)
data.dma.c4 <- list(y=y.dmon,nsites=nsites,ncap=ncap.dm,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c4n))
dma.zi.c4 <- autojags(data.dma.c4, inits.dmzi, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1, max.iter = 300000)
dma.re.c4 <- autojags(data.dma.c4, inits.dmzi, parameters, "re1.jags", n.chains = 3, n.thin = 2, parallel = TRUE, Rhat.limit = 1.1, max.iter = 300000)

# closeness centrality - 1:5 (upstream/downstream)
data.dma.c5 <- list(y=y.dmon,nsites=nsites,ncap=ncap.dm,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c5n))
dma.zi.c5 <- autojags(data.dma.c5, inits.dmzi, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1, max.iter = 300000)
dma.re.c5 <- autojags(data.dma.c5, inits.dmzi, parameters, "re1.jags", n.chains = 3, n.thin = 2, parallel = TRUE, Rhat.limit = 1.1, max.iter = 300000)

# closeness centrality - 1:6 (upstream/downstream)
data.dma.c6 <- list(y=y.dmon,nsites=nsites,ncap=ncap.dm,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c6n))
dma.zi.c6 <- autojags(data.dma.c6, inits.dmzi, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1, max.iter = 300000)
dma.re.c6 <- autojags(data.dma.c6, inits.dmzi, parameters, "re1.jags", n.chains = 3, n.thin = 2, parallel = TRUE, Rhat.limit = 1.1, max.iter = 300000)

# closeness centrality - 1:7 (upstream/downstream)
data.dma.c7 <- list(y=y.dmon,nsites=nsites,ncap=ncap.dm,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c7n))
dma.zi.c7 <- autojags(data.dma.c7, inits.dmzi, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1, max.iter = 300000)
dma.re.c7 <- autojags(data.dma.c7, inits.dmzi, parameters, "re1.jags", n.chains = 3, n.thin = 2, parallel = TRUE, Rhat.limit = 1.1, max.iter = 300000)

# closeness centrality - 1:8 (upstream/downstream)
data.dma.c8 <- list(y=y.dmon,nsites=nsites,ncap=ncap.dm,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c8n))
dma.zi.c8 <- autojags(data.dma.c8, inits.dmzi, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1, max.iter = 300000)
dma.re.c8 <- autojags(data.dma.c8, inits.dmzi, parameters, "re1.jags", n.chains = 3, n.thin = 2, parallel = TRUE, Rhat.limit = 1.1, max.iter = 300000)

# closeness centrality - 1:9 (upstream/downstream)
data.dma.c9 <- list(y=y.dmon,nsites=nsites,ncap=ncap.dm,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                    cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c9n))
dma.zi.c9 <- autojags(data.dma.c9, inits.dmzi, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1, max.iter = 300000)
dma.re.c9 <- autojags(data.dma.c9, inits.dmzi, parameters, "re1.jags", n.chains = 3, n.thin = 2, parallel = TRUE, Rhat.limit = 1.1, max.iter = 300000)

# closeness centrality - 1:10 (upstream/downstream)
data.dma.c10 <- list(y=y.dmon,nsites=nsites,ncap=ncap.dm,cov1=scale(desdata$PC1),cov2=scale(desdata$PC2),cov3=scale(desdata$PC3),
                     cov4=scale(desdata$rt),cov5=scale(desdata$t1),cov6=scale(desdata$c10n))
dma.zi.c10 <- autojags(data.dma.c10, inits.dmzi, parameters, "zi1.jags", n.chains=nc, n.thin=2, parallel=TRUE, Rhat.limit = 1.1, max.iter = 300000)
dma.re.c10 <- autojags(data.dma.c10, inits.dmzi, parameters, "re1.jags", n.chains = 3, n.thin = 2, parallel = TRUE, Rhat.limit = 1.1, max.iter = 300000)

t2 <- Sys.time()
