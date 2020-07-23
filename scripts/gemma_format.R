#!/usr/bin/env Rscript

library('data.table')
library('lme4')

phenotypes=fread('GridLMM/phenotypes_asi.csv',data.table=F)
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
envs=c('BLOIS_2014_OPT','BLOIS_2017_OPT','NERAC_2016_WD',
'GRANEROS_2015_OPT','STPAUL_2017_WD','SZEGED_2017_OPT')
phenos=c('male_flowering_d6','female_flowering_d6','male_flowering_days','female_flowering_days',
'total_plant_height','harvest_grain_moisture','grain_yield_15','tkw_15','asi')
for(p in phenos){
  for(e in envs){
    psub=phenotypes[phenotypes$Loc.Year.Treat==e,c('Genotype_code','Loc.Year.Treat',p)]
    bimbam=fread(cmd=sprintf('zcat genotypes/probabilities/allele_probs/bg10_wgs_alleleprobs.txt | head -n 4'),data.table=F)
    #sample_order=names(bimbam)[4:dim(bimbam)[2]]
    K=fread('GridLMM/K_matrices/K_matrix_chr1.txt',data.table=F)
    new_psub=data.frame(Genotype_code=K[,1],pheno=psub[match(K[,1],psub$Genotype_code),p],stringsAsFactors=F)
    fwrite(new_psub[,'pheno',drop=F],sprintf('gemma/phenotypes/%s_%s_phenotypes.txt',p,e),row.names=F,quote=F,col.names=F,na="NA")
    #for(c in 1:10){
    #  bimbam=fread(cmd=sprintf('zcat genotypes/probabilities/allele_probs/bg%.0f_wgs_alleleprobs.txt | head -n 4',c),data.table=F)
    #  sample_order=names(bimbam)[4:dim(bimbam)[2]]
    #  K=fread(sprintf('GridLMM/K_matrices/K_matrix_chr%.0f.txt',c),data.table=F)
    #  rownames(K)=K[,1]
    #  rownames(K)=gsub("-",".",rownames(K))
    #  K=as.matrix(K[,-1])
    #  colnames(K)=rownames(K)
    #  K=K[sample_order,sample_order]
    #  fwrite(K,sprintf('gemma/K_matrix_chr%.0f_%s_%s.txt',c,p,e),row.names=F,quote=F,col.names=F,sep='\t')
    #}
  }
  #BLUPS
  data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,p],stringsAsFactors=F)
  m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
  data_blup = as.data.frame(ranef(m1)$ID2)
  data_blup$ID = rownames(data_blup)
  data_blup$y=data_blup$`(Intercept)`
  data_blup=data_blup[,c('ID','y')]
  data_blup=data_blup[K[,1],]
  fwrite(data_blup[,'y',drop=F],sprintf('gemma/phenotypes/%s_BLUP_phenotypes.txt',p),row.names=F,quote=F,col.names=F,na="NA")
}
