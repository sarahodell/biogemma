#!/usr/bin/env Rscript

library('data.table')
library('ggplot2')
library('dplyr')

perm=fread('rep1000_max_pvalues.txt',data.table=F)

pdf('perm_pval_x_chr.pdf')
print(ggplot(perm,aes(x=-log10(pval))) + geom_histogram(bins=20,fill="darkgreen",color="black") + facet_wrap(~chr) + theme_classic() + ggtitle('1000 Permutation P-value Distribution by Chromosome'))
dev.off()


pdf('perm_pval_x_hapgrp.pdf')
print(ggplot(perm,aes(x=-log10(pval))) + geom_histogram(bins=20,fill="darkgreen",color="black") + geom_vline(aes(xintercept=mean(-log10(pval))),col='red',size=1)+ facet_wrap(~hapgrp) + theme_classic() + ggtitle('1000 Permutation P-value Distribution by Haplotype Group'))
dev.off()


