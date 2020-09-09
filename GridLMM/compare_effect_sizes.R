#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)


library('data.table')
library('ggplot2')



qtls=fread('Biogemma_QTL.csv',data.table=F)

# Compare founder and QTL
