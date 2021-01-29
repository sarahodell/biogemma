#!/usr/bin/env Rscript


library('data.table')
library('lme4')
library('lme4qtl')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

qtl=fread('GridLMM/Biogemma_QTL.csv',data.table=F)
qtl$pheno_env_id=paste0(qtl$pheno_env,'_',qtl$ID)
qtl=qtl[qtl$Method=="600K_SNP",]
#qtl=qtl[qtl$Phenotype==pheno & qtl$Environment==env & qtl$Chromosome==as.numeric(chr),]
rownames(qtl)=seq(1,nrow(qtl))
ses=list()
count=1
for(i in 1:nrow(qtl)){
  pheno=qtl[i,]$Phenotype
  env=qtl[i,]$Environment
  chr=as.character(qtl[i,]$Chromosome)
  snp=qtl[i,]$highest_SNP
  name=qtl[i,]$pheno_env_id

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
  #colnames(X) = founders
  #rownames(X) = dimnames(founder_probs[[1]])[[1]]
  #lowrep=which(apply(X,MARGIN=2,function(x) sum(x>0.8)<5))
  #if(length(lowrep)>0){
  #  X=X[,-lowrep]
  #}

  if(env != "ALL"){
    phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
    phenotype = phenotype[!is.na(phenotype$y),]
    phenotype$y=phenotype$y-mean(phenotype$y)
    rownames(phenotype)=seq(1,nrow(phenotype))
    #X = X[rownames(X) %in% phenotype$Genotype_code,,]
    phenotype$X=X[phenotype$Genotype_code,]
    m0 = relmatLmer(y ~ 0 + X + (1|Genotype_code),data=phenotype,
    relmat = list(Genotype_code=K))
    se=as.data.frame(summary(m0)$coef)
  }
  else{
    m1=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
    data_blup = as.data.frame(ranef(m1)$Genotype_code)
    data_blup$ID = rownames(data_blup)
    data_blup$y=data_blup$`(Intercept)`
    data_blup=data_blup[,c('ID','y')]
    data_blup$X=X[data_blup$ID,]
    #X = X[rownames(X) %in% data_blup$ID,,]
    #X=unlist(unname(X[,1]))
    m0 = relmatLmer(y ~ 0 + X + (1|ID),data=data_blup,relmat = list(ID=K))
    se=as.data.frame(summary(m0)$coef)
  }
  #if(dim(se)[1]!=16){
  #  rows=rownames(se)
  #  for(l in lowrep){
  #    se=rbind(se,rep(NA,3))
  #  }
  #  rownames(se)=c(rows,paste0('X',founders[lowrep]))
  #}
  #se=se[paste0('X',founders),]
  ses[[count]]=list(SE=se,qtl_id=name,snp=snp)
  count=count+1
    # pheno = cbind(pheno[match(rownames(X),pheno$Genotype_code),1:3],X_diff_mean)
    # model with no intercept
}

saveRDS(ses,'GridLMM/effect_sizes/600K_SNP_QTL_SEs.rds')

# model with intercept, dropping founder 1
#X_diff_mean = X %*% contr.sum(ncol(X))
#colnames(X_diff_mean) = paste0('F',1:ncol(X_diff_mean))
#m1 = relmatLmer(y ~ X_diff_mean + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=K))
#coefs1 = data.frame(summary(m1)$coef)[-1,]

# estimating SE for F16
#coefs1$Precision = 1/coefs1[,2]^2
#coefs1$nobs = colSums(X)[-16]
#m1_slope = lm(Precision~nobs,coefs1[-1,])
#with(coefs1,plot(Precision~nobs))
#abline(m1_slope)
#F16 = -sum(coefs1[,1])
#F16_SE = 1/sqrt(predict(m1_slope,newdata = list(nobs = colSums(X)[16])))

# model with intercept, dropping founders 1 and 16
#X_drop1_diff_mean = X[,-1] %*% contr.sum(ncol(X[,-1]))
#colnames(X_drop1_diff_mean) = paste0('F',1+1:ncol(X_drop1_diff_mean))
#m2 = relmatLmer(male_flowering_d6 ~ X_drop1_diff_mean + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=K))
#summary(m2)$coef
#coefs2 = data.frame(summary(m2)$coef)[-1,]

# estimating SE for F16
#coefs2$Precision = 1/coefs2[,2]^2
#coefs2$nobs = colSums(X)[-c(1,16)]
#m2_slope = lm(Precision~nobs,coefs2)
#with(coefs2,plot(Precision~nobs))
#abline(m2_slope)
#F16 = -sum(coefs2[,1])
#F16_SE = 1/sqrt(predict(m2_slope,newdata = list(nobs = colSums(X)[2])))

#fwrite(coefs,sprintf('GridLMM/effect_sizes/'))
