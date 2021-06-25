#!/usr/bin/env Rscript

library('data.table')

allreps=c()

for(r in 1:100){
  allres=c()
  for(c in 1:10){
    res=fread(sprintf('selection/bg%.0f_founder_chisq_rep%.0f_results.txt',c,r),data.table=F)
    allres=rbind(allres,res)
  }
  allreps=rbind(allreps,c(r,min(allres$p_chi)))
}
allreps=as.data.frame(allreps,stringsAsFactors=F)
names(allreps)=c('rep','minp')
thresh=0.05
threshold=quantile(allreps$minp,thresh,lower.tail=T)
#      5%
#8.828254
