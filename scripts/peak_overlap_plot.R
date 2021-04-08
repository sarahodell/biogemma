#!/usr/bin/env Rscript

library('ggplot2')
library('data.table')

qtl_info=fread('GridLMM/effect_sizes/all_qtl_info.txt',data.table=F)
qtl=fread('GridLMM/Biogemma_QTL.csv',data.table=F)

name="female_flowering_d6_NERAC_2016_WD_qDTS3_1"

a=qtl[qtl$pheno_env_id==name,]
b=qtl_info[qtl_info$pheno_env_id==name,]

data=data.frame(variable=c('spos','fpos','hpos'),value=c(unique(b$spos),unique(b$fpos),unique(b$hpos)),height=c(1,1,1),stringsAsFactors=F)

png(sprintf('GridLMM/effect_sizes/%s_method_peak_positions.png',name),width=800,height=800)
print(ggplot(data,aes(x=value,y=height,color=variable)) + geom_point())
dev.off()
