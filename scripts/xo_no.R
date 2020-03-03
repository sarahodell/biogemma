#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')

for(c in 1:10){
   chr=fread(sprintf('Biogemma_xo_number_chr%.0f.txt',c),data.table=F)
   png(sprintf('crossover_chr%0.0f.png',c))
   print(ggplot(chr,aes(x=xo_no)) + geom_histogram(fill="green") + ggtitle(sprintf("Chromosome %.0f Crossover Number",c)))
   dev.off()
}