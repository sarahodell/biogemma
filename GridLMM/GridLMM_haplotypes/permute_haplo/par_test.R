#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
cores=args[[1]]

library('doParallel')


data(iris)

c1 <- makeCluster(spec=rep("c9-76",cores),cores=cores)
registerDoParallel(c1,cores=cores)
library('doParallel')
library('foreach')
print(clusterEvalQ(c1, library(doParallel)))

x <- iris[which(iris[,5] != "setosa"), c(1,5)]
trials <- 10000
ptime <- system.time({
  r <- foreach(icount(trials), combine=cbind) %dopar% {
    ind <- sample(100,100, replace=T)
    result1 <- glm(x[ind,2]~x[ind,1],family=binomial(logit))
    coeffficients(result1)
    }
})[3]

print("In parallel")
print(ptime)



print("Sequentially")
stime <- system.time({
  r <- foreach(icount(trials), .combine=cbind) %do% {
    ind<-sample(100,100,replace=T)
    result1<-glm(x[ind,2]~x[ind,1],family=binomial(logit))
    coefficients(result1)
  }
})[3]
print(stime)



