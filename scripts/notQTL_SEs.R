library('data.table')
library('lme4')
library('lme4qtl')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

qtl=fread('GridLMM/Biogemma_QTL.csv',data.table=F)
qtl$pheno_env_id=paste0(qtl$pheno_env,'_',qtl$ID)

qtl_overlap=fread('GridLMM/Biogemma_Method_Overlap.csv',data.table=F)
qtl_overlap$pxe=paste0(qtl_overlap$Phenotype,'_',qtl_overlap$Environment)


f_nots=fread('GridLMM/effect_sizes/Founder_notSNP_highest_peaks.txt',data.table=F)
s_notf=fread('GridLMM/effect_sizes/SNP_notFounder_highest_peaks.txt',data.table=F)

f_noth=fread('GridLMM/effect_sizes/Founder_notHaplotype_highest_peaks.txt',data.table=F)
h_notf=fread('GridLMM/effect_sizes/Haplotype_notFounder_highest_peaks.txt',data.table=F)

s_notf_ids=unique(s_notf$pheno_env_id)
#m1=c('Founder_probs','600K_SNP')
ses=list()
count=1
for(q in s_notf_ids){
  line=qtl[qtl$pheno_env_id == q & qtl$Method=="600K_SNP",]
  #rownames(sub)=seq(1,nrow(sub))
  pheno=line$Phenotype
  env=line$Environment
  chr=as.character(line$Chromosome)
  snp=s_notf[s_notf$pheno_env_id==q,]$highest_SNP
  name=line$pheno_env_id
  founder_probs = readRDS(sprintf('genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
  K = fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  phenotype=fread('GridLMM/phenotypes_asi.csv',data.table=F)
  phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
  phenotype = phenotype[,c('Loc.Year.Treat','Genotype_code',pheno)]
  names(phenotype)=c('Loc.Year.Treat','Genotype_code','y')

  phenotype = subset(phenotype,Genotype_code %in% rownames(K))
  X = do.call(cbind,lapply(founder_probs,function(x) x[,snp]))
  colnames(X) = founders
  rownames(X) = dimnames(founder_probs[[1]])[[1]]
  lowrep=which(apply(X,MARGIN=2,function(x) sum(x>0.8)<5))
  if(length(lowrep)>0){
    X=X[,-lowrep]
  }


  if(env != "ALL"){
    phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
    phenotype = phenotype[!is.na(phenotype$y),]
    phenotype$y=phenotype$y-mean(phenotype$y)
    rownames(phenotype)=seq(1,nrow(phenotype))
    X = X[rownames(X) %in% phenotype$Genotype_code,]
    m0 = relmatLmer(y ~ 0 + X + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=K))
    se=as.data.frame(summary(m0)$coef)
  }
  else{
    m1=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
    data_blup = as.data.frame(ranef(m1)$Genotype_code)
    data_blup$ID = rownames(data_blup)
    data_blup$y=data_blup$`(Intercept)`
    data_blup=data_blup[,c('ID','y')]
    X = X[rownames(X) %in% data_blup$ID,]
    m0 = relmatLmer(y ~ 0 + X + (1|ID),data=data_blup,relmat = list(ID=K))
    se=as.data.frame(summary(m0)$coef)
  }

  if(dim(se)[1]!=16){
    rows=rownames(se)
    for(l in lowrep){
      se=rbind(se,rep(NA,3))
    }
    rownames(se)=c(rows,paste0('X',founders[lowrep]))
  }
  se=se[paste0('X',founders),]
  ses[[count]]=list(SE=se,qtl_id=name)
  count=count+1
}

saveRDS(ses,'GridLMM/effect_sizes/SNP_notFounder_prob_QTL_SEs.rds')

f_nots_ids=unique(f_nots$pheno_env_id)
#m1=c('Founder_probs','600K_SNP')
ses=list()
count=1
for(q in f_nots_ids){
  line=qtl[qtl$pheno_env_id == q & qtl$Method=="Founder_probs",]
  #rownames(sub)=seq(1,nrow(sub))
  pheno=line$Phenotype
  env=line$Environment
  chr=as.character(line$Chromosome)
  snp=f_nots[f_nots$pheno_env_id==q,]$highest_SNP
  name=line$pheno_env_id
  K = fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  phenotype=fread('GridLMM/phenotypes_asi.csv',data.table=F)
  phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
  phenotype = phenotype[,c('Loc.Year.Treat','Genotype_code',pheno)]
  names(phenotype)=c('Loc.Year.Treat','Genotype_code','y')
  X = fread(sprintf('genotypes/qtl2/Biogemma_DHgenos/DH_geno_chr%s_binary.csv',chr),data.table=F)
  rownames(X)=X$ind
  phenotype = subset(phenotype,Genotype_code %in% rownames(K))
  X = X[,snp,drop=F]

  if(env != "ALL"){
    phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
    phenotype = phenotype[!is.na(phenotype$y),]
    phenotype$y=phenotype$y-mean(phenotype$y)
    rownames(phenotype)=seq(1,nrow(phenotype))
    X = X[rownames(X) %in% phenotype$Genotype_code,,]
    m0 = relmatLmer(y ~ 1 + X + (1|Genotype_code),data=phenotype,
    relmat = list(Genotype_code=K))
    se=as.data.frame(summary(m0)$coef)
  }
  else{
    m1=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
    data_blup = as.data.frame(ranef(m1)$Genotype_code)
    data_blup$ID = rownames(data_blup)
    data_blup$y=data_blup$`(Intercept)`
    data_blup=data_blup[,c('ID','y')]
    X = X[rownames(X) %in% data_blup$ID,,]
    #X=unlist(unname(X[,1]))
    m0 = relmatLmer(y ~ 1 + X + (1|ID),data=data_blup,relmat = list(ID=K))
    se=as.data.frame(summary(m0)$coef)
  }
  ses[[count]]=list(SE=se,qtl_id=name,snp=snp)
  count=count+1
}

saveRDS(ses,'GridLMM/effect_sizes/Founder_notSNP_prob_QTL_SEs.rds')



h_notf_ids=unique(h_notf$pheno_env_id)
#m2=c('Founder_probs','Haplotype_probs')
ses=list()
count=1
for(q in h_notf_ids){
  line=qtl[qtl$pheno_env_id == q & qtl$Method=="Haplotype_probs",]
  #rownames(sub)=seq(1,nrow(sub))
  pheno=line$Phenotype
  env=line$Environment
  chr=as.character(line$Chromosome)
  snp=h_notf[h_notf$pheno_env_id==q,]$highest_SNP
  name=line$pheno_env_id
  founder_probs = readRDS(sprintf('genotypes/probabilities/geno_probs/bg%s_filtered_genotype_probs.rds',chr))
  K = fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  phenotype=fread('GridLMM/phenotypes_asi.csv',data.table=F)
  phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
  phenotype = phenotype[,c('Loc.Year.Treat','Genotype_code',pheno)]
  names(phenotype)=c('Loc.Year.Treat','Genotype_code','y')

  phenotype = subset(phenotype,Genotype_code %in% rownames(K))
  X = do.call(cbind,lapply(founder_probs,function(x) x[,snp]))
  colnames(X) = founders
  rownames(X) = dimnames(founder_probs[[1]])[[1]]
  lowrep=which(apply(X,MARGIN=2,function(x) sum(x>0.8)<5))
  if(length(lowrep)>0){
    X=X[,-lowrep]
  }
  if(env != "ALL"){
    phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
    phenotype = phenotype[!is.na(phenotype$y),]
    phenotype$y=phenotype$y-mean(phenotype$y)
    rownames(phenotype)=seq(1,nrow(phenotype))
    X = X[rownames(X) %in% phenotype$Genotype_code,]
    m0 = relmatLmer(y ~ 0 + X + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=K))
    se=as.data.frame(summary(m0)$coef)

  }
  else{
    m1=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
    data_blup = as.data.frame(ranef(m1)$Genotype_code)
    data_blup$ID = rownames(data_blup)
    data_blup$y=data_blup$`(Intercept)`
    data_blup=data_blup[,c('ID','y')]
    X = X[rownames(X) %in% data_blup$ID,]
    m0 = relmatLmer(y ~ 0 + X + (1|ID),data=data_blup,relmat = list(ID=K))
    se=as.data.frame(summary(m0)$coef)
  }
  if(dim(se)[1]!=16){
    rows=rownames(se)
    for(l in lowrep){
      se=rbind(se,rep(NA,3))
    }
    rownames(se)=c(rows,paste0('X',founders[lowrep]))
  }
  se=se[paste0('X',founders),]
  ses[[count]]=list(SE=se,qtl_id=name)
  count=count+1
}

saveRDS(ses,'GridLMM/effect_sizes/Haplotype_notFounder_prob_QTL_SEs.rds')

f_noth_ids=unique(f_noth$pheno_env_id)
#m2=c('Founder_probs','Haplotype_probs')
ses=list()
count=1
for(q in f_noth_ids){
  line=qtl[qtl$pheno_env_id == q & qtl$Method=="Founder_probs",]
  #rownames(sub)=seq(1,nrow(sub))
  pheno=line$Phenotype
  env=line$Environment
  chr=as.character(line$Chromosome)
  snp=f_noth[f_noth$pheno_env_id==q,]$highest_SNP
  name=line$pheno_env_id
  results=fread(sprintf('GridLMM/result_tables/Haplotype_GWAS_%s_results.txt',env),data.table=F)
  h=results[results$SNP==snp,]$HAPGRP
  haplo_probs = readRDS(sprintf('genotypes/probabilities/haplotype_probs/RefinedIBD_600K/bg%s_filtered_haplogroup%.0f_probs.rds',chr,h))
  K = fread(sprintf('GridLMM/K_matrices/K_matrix_chr%s.txt',chr),data.table = F,h=T)
  rownames(K) = K[,1]
  K = as.matrix(K[,-1])
  phenotype=fread('GridLMM/phenotypes_asi.csv',data.table=F)
  phenotype$Genotype_code = gsub('-','.',phenotype$Genotype_code,fixed=T)
  phenotype = phenotype[,c('Loc.Year.Treat','Genotype_code',pheno)]
  names(phenotype)=c('Loc.Year.Treat','Genotype_code','y')

  phenotype = subset(phenotype,Genotype_code %in% rownames(K))
  X = do.call(cbind,lapply(haplo_probs,function(x) x[,snp]))
  colnames(X) = paste0("HAPGRP_",seq(1,h))
  rownames(X) = dimnames(haplo_probs[[1]])[[1]]
  if(env != "ALL"){
    phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
    phenotype = phenotype[!is.na(phenotype$y),]
    phenotype$y=phenotype$y-mean(phenotype$y)
    rownames(phenotype)=seq(1,nrow(phenotype))
    X = X[rownames(X) %in% phenotype$Genotype_code,]
    m0 = relmatLmer(y ~ 0 + X + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=K))
    ses[[count]]=list(SE=summary(m0)$coef,qtl_id=name)
    count=count+1
  }
  else{
    m1=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
    data_blup = as.data.frame(ranef(m1)$Genotype_code)
    data_blup$ID = rownames(data_blup)
    data_blup$y=data_blup$`(Intercept)`
    data_blup=data_blup[,c('ID','y')]
    X = X[rownames(X) %in% data_blup$ID,]
    m0 = relmatLmer(y ~ 0 + X + (1|ID),data=data_blup,relmat = list(ID=K))
    ses[[count]]=list(SE=summary(m0)$coef,qtl_id=name)
    count=count+1
  }
}

saveRDS(ses,'GridLMM/effect_sizes/Founder_notHaplotype_prob_QTL_SEs.rds')
