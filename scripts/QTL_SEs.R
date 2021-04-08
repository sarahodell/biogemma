#!/usr/bin/env Rscript


library('data.table')
library('lme4')
library('lme4qtl')
library('emmeans')
library('ggplot2')
library('GridLMM')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
has_mite=c(F,T,T,T,F,T,T,T,T,T,T,F,T,T,T,F)

qtl=fread('GridLMM/Biogemma_QTL.csv',data.table=F)
qtl$pheno_env_id=paste0(qtl$pheno_env,'_',qtl$ID)
qtl=qtl[qtl$Method=="Founder_probs",]

colorcodes=fread('GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder
colorcodes=colorcodes[founders,]
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
  print(name)
  pmap=fread(sprintf('genotypes/qtl2/startfiles/Biogemma_pmap_c%s.csv',chr),data.table=F)
  pos=pmap[pmap$marker==snp,]$pos
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
  #lowrep=which(apply(X,MARGIN=2,function(x) sum(x>0.8)<5))
  #if(length(lowrep)>0){
  #  X=X[,-lowrep]
  #}

  if(env != "ALL"){
    phenotype=phenotype[phenotype$Loc.Year.Treat==env,]
    phenotype = phenotype[!is.na(phenotype$y),]
    phenotype$y=phenotype$y-mean(phenotype$y)
    #rownames(phenotype)=seq(1,nrow(phenotype))
    #X = X[rownames(X) %in% phenotype$Genotype_code,]
    rownames(phenotype)=phenotype$Genotype_code
    i=intersect(phenotype$Genotype_code,rownames(X))
    X = X[i,]
    phenotype=phenotype[i,]
    subK=K[i,i]
    #drop2=which(colSums(X)<=5)
    #if(length(drop2)!=0){
    #  X_filt=X[,-drop2]
    #  f2=founders[-drop2]
    #}
    #else{
    #  X_filt=X
    #  f2=founders
    #}
    #phenotype=subset(phenotype,Genotype_code %in% rownames(X))
    #K=K[rownames(K) %in% rownames(X_filt),colnames(K) %in% rownames(X_filt)]
    #fmax=apply(X_filt,MARGIN=1,which.max)
    #phenotype$founder=f2[fmax]
    #phenotype$founder_f=factor(phenotype$founder,levels=f2)
    #dont drop anything
    #fmax=apply(X,MARGIN=1,which.max)
    #phenotype$founder=founders[fmax]
    #phenotype$founder_f=factor(phenotype$founder,levels=founders)
    m4 = relmatLmer(y ~ 0 + X + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=subK))
    se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
    names(se4)=c('value','se','tvalue')
    rownames(se4)=founders
    se4$founder=rownames(se4)
    se4$variable_f=factor(se4$founder,levels=se4$founder)
    #se4[2,]$value=se4[-1,]$value + se4[1,]$value
    #GridLMM
    null_model = GridLMM_ML(y~1 + (1|Genotype_code),phenotype,relmat=list(Genotype_code=subK),ML=T,REML=F,verbose=T)

    h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
    names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
    h2_start
    V_setup=null_model$setup
    vars = as.data.frame(VarCorr(m4))
    h2_start = vars$vcov[1]/sum(vars$vcov)
    names(h2_start) = 'Genotype_code'

    Y=as.matrix(phenotype$y)
    X_cov=null_model$lmod$X
    #X_cov is n x 0 matrix and can use full X_list_ordered
    #X_list_ordered=lapply(X_list,function(x) x[data_blup$ID,,drop=F])
    X_list_ordered = lapply(1:ncol(X),function(i) X[,i,drop=F])
    X_list_null=NULL
    #X_list_null
    cores=1
    gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=T,h2_step=1  )
    betas = unlist(gwas[1,grep('beta',colnames(gwas))])
    betas[-1] = betas[-1]+betas[1]
    corr=cor(betas,se4$value)

    png(sprintf('GridLMM/effect_sizes/founder_ES/%s_founder_effect_sizes_lme4qtl.png',name),width=1000,height=800)
    print(ggplot(se4,aes(x=variable_f,y=value,color=variable_f)) +
    geom_point() +
     geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
     scale_color_manual(values=colorcodes[levels(se4$variable_f),]$hex_color,labels=levels(se4$variable_f))+
    theme(axis.text.x=element_text(size=10)) +
    xlab("Founder") + ylab("Effect Size") +
    labs(title=sprintf("%s Effect Size Estimates",name)))
    dev.off()
  }
  else{
    phenotype = phenotype[!is.na(phenotype$y),]
    #phenotype$y=phenotype$y-mean(phenotype$y)
    m0=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
    data_blup = as.data.frame(ranef(m0)$Genotype_code)
    data_blup$ID = rownames(data_blup)
    data_blup$y=data_blup$`(Intercept)`
    data_blup=data_blup[,c('ID','y')]
    i=intersect(data_blup$ID,rownames(X))
    X = X[i,]
    data_blup=data_blup[i,]
    subK=K[i,i]
    #X = X[match(data_blup$ID, rownames(X)),]
    #drop2=which(colSums(X)<=5)
    #if(length(drop2)!=0){
    #  X_filt=X[,-drop2]
    #  f2=founders[-drop2]
    #}else{
    #  X_filt=X
    #  f2=founders
    #}
    #data_blup=subset(data_blup,ID %in% rownames(X_filt))
    #K=K[rownames(K) %in% rownames(X_filt),colnames(K) %in% rownames(X_filt)]

    m4 = relmatLmer(y ~ 0 + X + (1|ID),data=data_blup,relmat = list(ID=subK),REML=T)
    se4=as.data.frame(summary(m4)$coef,stringsAsFactors=F)
    rownames(se4)=founders
    names(se4)=c('value','se','tvalue')
    se4$founder=rownames(se4)
    se4$variable_f=factor(se4$founder,levels=se4$founder)
    #se4[2,]$value=se4[-1,]$value + se4[1,]$value
    #se4=se4[order(se4$value),]
    #GRIDLMM
    null_model = GridLMM_ML(y~1 + (1|ID),data_blup,relmat=list(ID=subK),ML=T,REML=F,verbose=F)

    h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
    names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
    h2_start
    V_setup=null_model$setup
    vars = as.data.frame(VarCorr(m4))
    h2_start = vars$vcov[1]/sum(vars$vcov)
    names(h2_start) = 'ID'


    Y=as.matrix(data_blup$y)
    X_cov=null_model$lmod$X
    #X_cov is n x 0 matrix and can use full X_list_ordered
    #X_list_ordered=lapply(X_list,function(x) x[data_blup$ID,,drop=F])
    X_list_ordered = lapply(1:ncol(X),function(i) X[,i,drop=F])
    X_list_null=NULL
    #X_list_null
    cores=1
    gwas=run_GridLMM_GWAS(Y,X_cov,X_list_ordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=T,h2_step=1  )
    betas = unlist(gwas[1,grep('beta',colnames(gwas))])
    betas[-1] = betas[-1]+betas[1]
    corr=cor(betas,se4$value)
    #

    #fmax=apply(X_filt,MARGIN=1,which.max)
    #data_blup$founder=f2[fmax]
    #data_blup$founder_f=factor(data_blup$founder,levels=f2)

    #m4 = relmatLmer(y ~ 0 + founder_f + (1|ID),data=data_blup,relmat = list(ID=K))
    #se4=as.data.frame(emmeans(m4,'founder_f'),stringsAsFactors=F)
    #se4=se4[order(se4$emmean),]
    #rownames(se4)=seq(1,nrow(se4))
    #se4$variable=sapply(seq(1,nrow(se4)),function(x) strsplit(rownames(se4)[x],'founder')[[1]][2])
    #se4$variable_f=factor(se4$founder,levels=se4$founder)
    #se4$mite=has_mite[match(se4$founder,founders)]

    png(sprintf('GridLMM/effect_sizes/founder_ES/%s_BLUP_founder_effect_sizes_lme4qtl.png',name),width=1000,height=800)
    print(ggplot(se4,aes(x=variable_f,y=value,color=variable_f)) +
    geom_point() +
    geom_errorbar(aes(ymin=value - 2*se,ymax=value + 2*se)) +
    scale_color_manual(values=colorcodes[levels(se4$variable_f),]$hex_color,labels=levels(se4$variable_f))+
    theme(axis.text.x=element_text(size=10)) +
    xlab("Founder") + ylab("Effect Size") +
    labs(title=sprintf("%s BLUP Effect Size Estimates",name)))
    dev.off()
  }

  #if(dim(se4)[1]!=16){
  #  rows=rownames(se4)
  #  for(l in lowrep){
  #    se4=rbind(se4,rep(NA,3))
  #  }
  #  rownames(se4)=c(rows,paste0('X',founders[lowrep]))
  #}
  #se4=se4[paste0('X',founders),]
  ses[[count]]=list(SE=se4,qtl_id=name,snp=snp,pos=pos,corr=corr)
  count=count+1
    # pheno = cbind(pheno[match(rownames(X),pheno$Genotype_code),1:3],X_diff_mean)
    # model with no intercept
}

saveRDS(ses,'GridLMM/effect_sizes/founder_ES/Founder_prob_QTL_SEs.rds')

# model with intercept, dropping founder 1
#X_diff_mean = X %*% contr.sum(ncol(X))
#colnames(X_diff_mean) = founders[2:16]
#X_diff_mean=X_diff_mean[rownames(X_diff_mean) %in% data_blup$ID,]
#m1 = relmatLmer(y ~ X_diff_mean + (1|ID),data=data_blup,relmat = list(ID=K))
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
#snp="AX-91102858"
# Using only very certain lines




#For single environments
