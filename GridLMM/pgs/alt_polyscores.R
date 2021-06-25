#!/usr/bin/env Rscript

library('data.table')
library('dplyr')
library("xlsx")

data <- read.xlsx('Castelletti2020_effectsizes.xlsx', 1)
es=data[data$Trait=="Anthesis",]
es=data[data$Trait=="Anthesis" & data$Experiment=="Field 2015",]

#markers=c()
#for(c in 1:10){
#  map=fread(sprintf('../../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',c),data.table=F)
#  markers=c(markers,map$marker)
#}



genos=fread('../../genotypes/qtl2/Biogemma_DHgenos/rrBLUP_genotype.txt',data.table=F)
m=intersect(es$SNP.Name, genos$marker)

genos=genos[genos$marker %in% m,]
rownames(genos)=seq(1,nrow(genos))

es=es[es$SNP.Name %in% m,]
# effect sizes are for reference allele, not alternate (make for alternate allele)
es$effect_size=es$Allelic.effect * -1

us=fread('../../run_magicsim/pgs/marker_effect_sizes.txt',data.table=F)
us=us[us$marker %in% es$SNP.Name,]

X=genos[,4:347]
rownames(X)=genos$marker
X=t(X)
u=as.matrix(es$effect_size)
res=X %*% u

res=as.data.frame(res,stringsAsFactors=F)
names(res)=c('dta')
res$ind=rownames(res)
fwrite(res,'Castelletti_DTA_FT_PGS.txt',row.names=F,quote=F,sep='\t')

i=intersect(colnames(X),es$SNP.Name)
u2=us[match(i,us$marker),]$effect
u2=as.matrix(u2)
  #Zs=Zs[,i]
#X=X[,i]
res2=X %*% u2

res2=as.data.frame(res2,stringsAsFactors=F)
names(res2)=c('dta')
res2$ind=rownames(res2)
fwrite(res2,'Castelletti_markers_BG_DTA_FT_PGS.txt',row.names=F,quote=F,sep='\t')


#us=as.data.frame(us,stringsAsFactors=F)
#names(us)=c('marker','effect','chr')
fwrite(es,'../../run_magicsim/pgs/Castelletti_DTA_marker_effect_sizes.txt',row.names=F,quote=F,sep='\t')
