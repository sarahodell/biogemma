#!/usr/bin/env Rscript
library('data.table')
library('dplyr')

args=commandArgs(trailingOnly=T)
env=as.character(args[[1]])

#env="BLOIS_2014_OPT"
phenotypes=c("male_flowering_d6","female_flowering_d6","total_plant_height","grain_yield_15","tkw_15","harvest_grain_moisture")


full_gwas=c()
for(i in seq(1,10)){
      tmp_gwas=c()
      for(p in phenotypes){
      	    g=readRDS(sprintf('GridLMM_founderprobs/models/Biogemma_chr%.0f_%s_x_%s_founderprobs.rds',i,p,env))
  	    pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',i),data.table=F)
	    gwas=g[,c('X_ID','p_value_ML')]
  	    gwas$chr=i
  	    m=match(gwas$X_ID,pmap$marker)
 	    gwas$pos=pmap$pos[m]
	    pval=sprintf("%s_P",p)
  	    names(gwas)=c("SNP",pval,"CHR","BP")
  	    gwas=gwas[complete.cases(gwas),]
	    if(is.null(dim(tmp_gwas))){
		tmp_gwas=gwas
	    }else{
	    	tmp_gwas=left_join(tmp_gwas,gwas[,c('SNP',pval)],by="SNP",keep=F)
	    }
     }
     full_gwas=rbind(full_gwas,tmp_gwas)
}

fwrite(full_gwas,sprintf('result_tables/Founder_GWAS_%s_results.txt',env),row.names=F,quote=F,sep='\t')
