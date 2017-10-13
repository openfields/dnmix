# model file
sink("pt.jags")
cat("
    model {
    for(i in 1:nsites){
      p[i] ~ dunif(0,.4)
    }
    beta0 ~ dnorm(0,.1) # intercept
    beta1 ~ dnorm(0,.1) # pc1
    beta2 ~ dnorm(0,.1) # pc2
    beta3 ~ dnorm(0,.1) # pc3
    beta4 ~ dnorm(0,.1) # RT
    beta5 ~ dnorm(0,.1) # T1
    beta6 ~ dnorm(0,.1) # centrality parameter 1

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

# read data
read.csv('./data/dfadma_merge.csv', header=TRUE) -> desdata
desdata[,10:12]-> y.dfus
desdata[,26:28]-> y.dmon
y.dmon[is.na(y.dmon)] <- 0

# get sum of captures
apply(y.dmon,1,sum) -> ymax.dm
apply(y.dfus,1,sum) -> ymax.df

# compile data objects for fuscus and monticola
data.dfa.c7 <- list(y=y.dfus, nsites=dim(y.dfus)[1], ncap=ymax.df, cov1=scale(desdata$PC1), cov2=scale(desdata$PC2), cov3=scale(desdata$PC3), cov4=scale(desdata$rt), cov5=scale(desdata$t1), cov6=scale(desdata$c7n))
data.dma.c2 <- list(y=y.dmon, nsites=dim(y.dmon)[1], ncap=ymax.dm, cov1=scale(desdata$PC1), cov2=scale(desdata$PC2), cov3=scale(desdata$PC3), cov4=scale(desdata$rt), cov5=scale(desdata$t1), cov6=scale(desdata$c2n))

# generate starting values for models
inits.dm <- function(){list(N=ymax.dm+1)}
inits.df <- function(){list(N=ymax.df+1)}

# fit model, check fit
# closeness 7:1, fuscus
dfa.pt.c7 <- autojags(data.dfa.c7, inits.df, params.pt, "pt.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dfa.pt.c7)
mf7c_bpv <- bpv(dfa.pt.c7, "fuscus model: closeness (1:7)")

# closeness 2:1, monticola
dma.pt.c2 <- autojags(data.dma.c2, inits.dm, params.pt, "pt.jags", n.chains=nc, n.thin=1, parallel=TRUE)
print(dma.pt.c2)
mf7c_bpv <- bpv(dfa.pt.c7, "monticola model: closeness (1:2)")
