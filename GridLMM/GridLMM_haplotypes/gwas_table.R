#!/usr/bin/env Rscript
library('data.table')
library('dplyr')

args=commandArgs(trailingOnly=T)
env=as.character(args[[1]])
phenotypes=c("male_flowering_d6","female_flowering_d6","total_plant_height","grain_yield_15","tkw_15","harvest_grain_moisture")
base=c(6,10,6,7,9,9,9,9,8,7)

full_gwas=c()
for(i in seq(1,10)){
      pmap=fread(sprintf('../qtl2_startfiles/Biogemma_pmap_c%.0f.csv',i),data.table=F)
      for(h in seq(base[i],16)){
      	    tmp_gwas=c()
	    for(p in phenotypes){
      	    	  g=readRDS(sprintf('models/Biogemma_chr%.0f_haplogrp%.0f_%s_x_%s.rds',i,h,p,env))
		  gwas=g[,c('X_ID','p_value_ML')]
		  gwas$chr=i
		  gwas$hapgrp=h
		  m=match(gwas$X_ID,pmap$marker)
		  gwas$pos=pmap$pos[m]
		  pval=sprintf("%s_P",p)
		  names(gwas)=c("SNP",pval,"CHR","HAPGRP","BP")
		  gwas=gwas[complete.cases(gwas),]
	    	  if(is.null(dim(tmp_gwas))){
			tmp_gwas=gwas
	    	  }else{
	    	  tmp_gwas=left_join(tmp_gwas,gwas[,c('SNP',pval)],by="SNP",keep=F)
		  }
	    }
	    full_gwas=rbind(full_gwas,tmp_gwas)
     }
}

fwrite(full_gwas,sprintf('result_tables/Haplotype_GWAS_%s_results.txt',env),row.names=F,quote=F,sep='\t')