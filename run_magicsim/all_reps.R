#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')



sig=c()
for(r in 1:100){
  total_tests=0
  tmpsig=c()
  for(c in 1:10){
    pchi=fread(sprintf('selection/bg%.0f_founder_chisq_rep%.0f_results.txt',c,r),data.table=F)
    tmpsig=rbind(tmpsig,pchi)
    total_tests=total_tests+nrow(pchi)
  }
  total_tests=total_tests*16
  bonf=-log10(0.05 / total_tests)
  s=tmpsig$p_chi[-log10(tmpsig$p_chi)>=bonf]
  sig=rbind(sig,c(r,length(s)))
}
