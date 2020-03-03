#!/usr/bin/env Rscript

library('data.table')

for(i in 1:10){
    pmap=fread(sprintf('Biogemma_pmap_c%.0f.csv',i),data.table=F)
    pmap$pos=pmap$pos*1e6
    fwrite(pmap,sprintf('Biogemma_pmap_c%.0f.csv',i),row.names=F,quote=F,sep=',')
}