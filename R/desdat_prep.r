# data prep file
# compile data for use with models
# read file
read.csv("./data/dfadma_merge.csv", header=TRUE) -> desdat
str(desdat)

# create y matrices for observations
y.dm <- cbind(desdat$ma1,desdat$ma2,desdat$ma3)
y.dm[is.na(y.dm)]<-0

y.df <- cbind(desdat$fa1,desdat$fa2,desdat$fa3)

# create data object
data.df.c7 <- list(y=y.df, cov1=scale(desdat$PC1), cov2=scale(desdat$PC2), cov3=scale(desdat$PC3), cov4=scale(desdat$rt), cov5=scale(desdat$t1), cov6=scale(desdat$c7n), cov7=scale(desdat$bcns))

data.dm.c2 <- list(y=y.dm, cov1=scale(desdat$PC1), cov2=scale(desdat$PC2), cov3=scale(desdat$PC3), cov4=scale(desdat$rt), cov5=scale(desdat$t1), cov6=scale(desdat$c2n), cov7=scale(desdat$bcns))
