#!/usr/bin/env Rscript
library('data.table')
library('dplyr')

args=commandArgs(trailingOnly=T)
env=as.character(args[[1]])
f_envs=c("BLOIS_2014_OPT","GRANEROS_2015_OPT","NERAC_2016_WD","ALL")
if(env %in% f_envs){
  phenotypes=c("male_flowering_d6","male_flowering_days","female_flowering_d6","female_flowering_days","total_plant_height","grain_yield_15","tkw_15","harvest_grain_moisture","asi")

}else{
  phenotypes=c("male_flowering_d6","female_flowering_d6","total_plant_height","grain_yield_15","tkw_15","harvest_grain_moisture","asi")
}
#base=c(10,9,9,10,7,10,6,9,9,10)
base=c(8,7,7,8,6,7,7,8,7,7)

full_gwas=c()
for(i in seq(1,10)){
      pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',i),data.table=F)
      gmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_gmap_c%.0f.csv',i),data.table=F)
      for(h in seq(base[i],16)){
      	    tmp_gwas=c()
	    for(p in phenotypes){
      	    	  g=readRDS(sprintf('GridLMM_haplotypes/models/Biogemma_chr%.0f_haplogrp%.0f_%s_x_%s.rds',i,h,p,env))
		  gwas=g[,c('X_ID','p_value_ML')]
		  gwas$chr=i
		  gwas$hapgrp=h
		  m=match(gwas$X_ID,pmap$marker)
      m2=match(gwas$X_ID,gmap$marker)
		  gwas$pos=pmap$pos[m]
      gwas$cM=gmap$pos[m2]
		  pval=sprintf("%s_P",p)
		  names(gwas)=c("SNP",pval,"CHR","HAPGRP","BP","cM")
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
