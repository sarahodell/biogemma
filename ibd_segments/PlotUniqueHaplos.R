#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])

library('data.table')
library('ggplot2')

ibd_segments=fread(sprintf('bg%s_ibd_blocks.txt',chr),data.table=F)

png(sprintf('chr%s_unique_haps.png',chr),width=960,height=480)
print(ggplot(ibd_segments, aes(start/1e6, n_haps)) + geom_segment(aes(xstart=start/1e6,xend = end/1e6,color="red", ystart=n_haps,yend=n_haps),lineend="butt",size=10) + ggtitle(sprintf("Number of Unique Chromosome %s Haplotypes",chr)) + xlab("Position (Mb)") + ylab("Haplotype Number")+guides(color=F))
dev.off()
