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


full_gwas=c()
for(i in seq(1,10)){
      tmp_gwas=c()
      for(p in phenotypes){
      	    g=readRDS(sprintf('GridLMM_600KSNP/models/chr%.0f_%s_x_%s_600KSNP_ML.rds',i,p,env))$results
  	    pmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_pmap_c%.0f.csv',i),data.table=F)
        gmap=fread(sprintf('../genotypes/qtl2/startfiles/Biogemma_gmap_c%.0f.csv',i),data.table=F)
	    gwas=g[,c('X_ID','p_value_ML')]
  	    gwas$chr=i
  	    m=match(gwas$X_ID,pmap$marker)
        m2=match(gwas$X_ID,gmap$marker)
 	    gwas$pos=pmap$pos[m]
      gwas$cM=gmap$pos[m2]
	    pval=sprintf("%s_P",p)
  	    names(gwas)=c("SNP",pval,"CHR","BP","cM")
  	    gwas=gwas[complete.cases(gwas),]
	    if(is.null(dim(tmp_gwas))){
		tmp_gwas=gwas
	    }else{
	    	tmp_gwas=left_join(tmp_gwas,gwas[,c('SNP',pval)],by="SNP",keep=F)
	    }
     }
     full_gwas=rbind(full_gwas,tmp_gwas)
}

fwrite(full_gwas,sprintf('result_tables/600K_GWAS_%s_results.txt',env),row.names=F,quote=F,sep='\t')
