#!/usr/bin/env Rscript


library('data.table')
library('lme4')
library('lme4qtl')
library('emmeans')
library('ggplot2')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
has_mite=c(F,T,T,T,F,T,T,T,T,T,T,F,T,T,T,F)

qtl=fread('GridLMM/Biogemma_QTL.csv',data.table=F)
qtl$pheno_env_id=paste0(qtl$pheno_env,'_',qtl$ID)
qtl=qtl[qtl$Method=="Founder_probs",]
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
    phenotype = phenotype[!is.na(phenotype$y),]
    phenotype$y=phenotype$y-mean(phenotype$y)
    m0=lmer(y~Loc.Year.Treat + (1|Genotype_code),phenotype)
    data_blup = as.data.frame(ranef(m0)$Genotype_code)
    data_blup$ID = rownames(data_blup)
    data_blup$y=data_blup$`(Intercept)`
    data_blup=data_blup[,c('ID','y')]
    X = X[rownames(X) %in% data_blup$ID,]
    m1 = relmatLmer(y ~ 0 + X + (1|ID),data=data_blup,relmat = list(ID=K))
    se=as.data.frame(summary(m1)$coef)
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
    # pheno = cbind(pheno[match(rownames(X),pheno$Genotype_code),1:3],X_diff_mean)
    # model with no intercept
}

saveRDS(ses,'GridLMM/effect_sizes/Founder_prob_QTL_SEs.rds')

# model with intercept, dropping founder 1
X_diff_mean = X %*% contr.sum(ncol(X))
colnames(X_diff_mean) = founders[2:16]
X_diff_mean=X_diff_mean[rownames(X_diff_mean) %in% data_blup$ID,]
m1 = relmatLmer(y ~ X_diff_mean + (1|ID),data=data_blup,relmat = list(ID=K))
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

X_clean=t(sapply(seq(1,nrow(X)),function(x) ifelse(X[x,]>0.85,1,0)))
rownames(X_clean)=rownames(X)
drop=which(rowSums(X_clean)<1)
X_filt=X_clean[-drop,]
X_filt=X_filt[,-5]
phenotype=phenotype[rownames(X_filt),]
#data_blup=data_blup[data_blup$ID %in% rownames(X_filt),]
K=K[rownames(K) %in% rownames(X_filt),colnames(K) %in% rownames(X_filt)]
m2 = relmatLmer(y ~ 0 + X_filt + (1|ID),data=data_blup,relmat = list(ID=K))
se2=as.data.frame(summary(m2)$coef)
names(se2)=c('value','se','tvalue')
se2=se2[order(se2$value),]

f2=founders[founders!="OH43_inra"]
#Using discrete founder variables
fmax=apply(X_filt,MARGIN=1,which.max)
data_blup$founder=f2[fmax]
data_blup$founder_f=factor(data_blup$founder,levels=f2)
m3 = relmatLmer(y ~ 0 + founder + (1|ID),data=data_blup,relmat = list(ID=K))
se3=as.data.frame(summary(m3)$coef)
names(se3)=c('value','se','tvalue')
se3=se3[order(se3$value),]

#Drop A654 as well
d2=data_blup[data_blup$founder !="A654_inra",]
rownames(d2)=seq(1,nrow(d2))
K2=K[rownames(K) %in% d2$ID, colnames(K) %in% d2$ID]

m4 = relmatLmer(y ~ 0 + founder + (1|ID),data=d2,relmat = list(ID=K2))
se4=as.data.frame(summary(m4)$coef)
names(se4)=c('value','se','tvalue')
se4=se4[order(se4$value),]
se4$variable=sapply(seq(1,nrow(se4)),function(x) strsplit(rownames(se4)[x],'founder')[[1]][2])
se4$variable_f=factor(se4$variable,levels=se4$variable)
se4$mite=has_mite[match(se4$variable,founders)]

se4=as.data.frame(emmeans(m4,'founder'),stringsAsFactors=F)
se4=se4[order(se4$emmean),]
rownames(se4)=seq(1,nrow(se4))
#se4$variable=sapply(seq(1,nrow(se4)),function(x) strsplit(rownames(se4)[x],'founder')[[1]][2])
se4$variable_f=factor(se4$founder,levels=se4$founder)
se4$mite=has_mite[match(se4$founder,founders)]


png('GridLMM/effect_sizes/vgt1_BLUP_effect_sizes_lme4qtl.png',width=1000,height=800)
print(ggplot(se4,aes(x=variable_f,y=emmean,color=mite)) +
geom_point() +
 geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL)) +
theme(axis.text.x=element_text(size=10)) +
xlab("Founder") + ylab("Effect Size (ggd)") +
labs(title="Days to Anthesis BLUP Effect Size Estimates",subtitle="GWAS_SNP Alt. Allele SE=-4.548011, StdEr=1.32739",color="MITE Present"))
dev.off()

# Using 1 intercept
fmax=apply(X_filt,MARGIN=1,which.max)
data_blup$founder=f2[fmax]
data_blup$founder_f=factor(data_blup$founder,levels=f2)
m5 = relmatLmer(y ~ 1 + founder_f + (1|ID),data=data_blup,relmat = list(ID=K))
se5=as.data.frame(summary(m5)$coef)
names(se5)=c('value','se','tvalue')
se5=se5[order(se5$value),]

#fmax=apply(X_filt,MARGIN=1,which.max)
#data_blup$founder=f2[fmax]
#data_blup$founder_f=factor(data_blup$founder,levels=f2)
m6 = relmatLmer(y ~ founder_f + (1|ID),data=data_blup,relmat = list(ID=K))
se6=as.data.frame(summary(m6)$coef)
names(se6)=c('value','se','tvalue')
se6=se6[order(se6$value),]

#For single environments
rownames(phenotype)=phenotype$Genotype_code
#X_clean=t(sapply(seq(1,nrow(X)),function(x) ifelse(X[x,]>0.80,1,0)))
#rownames(X_clean)=rownames(X)
#drop=which(rowSums(X)<.90)
#X_filt=X_clean[-drop,]
drop2=which(colSums(X)<=5)
if(length(drop2)!=0){
  X_filt=X[,-drop2]
}
phenotype=subset(phenotype,Genotype_code %in% rownames(X_filt))
#phenotype=phenotype[rownames(phenotype) %in% rownames(X_filt),]
#phenotype=phenotype[rownames(X_filt),]
#data_blup=data_blup[data_blup$ID %in% rownames(X_filt),]
K=K[rownames(K) %in% rownames(X_filt),colnames(K) %in% rownames(X_filt)]
m2 = relmatLmer(y ~ 0 + X_filt + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=K))
se2=as.data.frame(summary(m2)$coef)
names(se2)=c('value','se','tvalue')
se2=se2[order(se2$value),]

f2=founders[-drop2]
#f2=founders[founders!="OH43_inra"]
#Using discrete founder variables
fmax=apply(X_filt,MARGIN=1,which.max)
phenotype$founder=f2[fmax]
phenotype$founder_f=factor(phenotype$founder,levels=f2)
m3 = relmatLmer(y ~ 0 + founder + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=K))
se3=as.data.frame(summary(m3)$coef)
names(se3)=c('value','se','tvalue')
se3=se3[order(se3$value),]

#Drop A654 as well
d2=data_blup[data_blup$founder !="A654_inra",]
rownames(d2)=seq(1,nrow(d2))
K2=K[rownames(K) %in% d2$ID, colnames(K) %in% d2$ID]

se3=as.data.frame(emmeans(m3,'founder'),stringsAsFactors=F)
#se3$founder=sort(f2)
se3=se3[order(se3$emmean),]
rownames(se3)=seq(1,nrow(se3))
#se4$variable=sapply(seq(1,nrow(se4)),function(x) strsplit(rownames(se4)[x],'founder')[[1]][2])
se3$variable_f=factor(se3$founder,levels=se3$founder)
#se3$mite=has_mite[match(se3$founder,founders)]



#dont drop anything
fmax=apply(X,MARGIN=1,which.max)
phenotype$founder=founders[fmax]
phenotype$founder_f=factor(phenotype$founder,levels=founders)
m4 = relmatLmer(y ~ 0 + founder + (1|Genotype_code),data=phenotype,relmat = list(Genotype_code=K))
se4=as.data.frame(summary(m4)$coef)
names(se4)=c('value','se','tvalue')
se4=se4[order(se4$value),]
se4=as.data.frame(emmeans(m4,'founder'),stringsAsFactors=F)
#se3$founder=sort(f2)
se4=se4[order(se4$emmean),]
rownames(se4)=seq(1,nrow(se4))
#se4$variable=sapply(seq(1,nrow(se4)),function(x) strsplit(rownames(se4)[x],'founder')[[1]][2])
se4$founder=as.character(se4$founder)
se4$variable_f=factor(se4$founder,levels=se4$founder)

colorcodes=fread('GridLMM/effect_sizes/founder_color_codes.txt',data.table=F)
rownames(colorcodes)=colorcodes$founder
colorcodes=colorcodes[founders,]

png('GridLMM/effect_sizes/qDTA3_1_BLOIS_2017_OPT_effect_sizes_lme4qtl.png',width=1000,height=800)
print(ggplot(se4,aes(x=variable_f,y=emmean,color=variable_f)) +
geom_point() +
 geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL)) +
 scale_color_manual(values=colorcodes[levels(se4$variable_f),]$hex_color,labels=levels(se4$variable_f))+
theme(axis.text.x=element_text(size=10)) +
xlab("Founder") + ylab("Effect Size (ggd)") +
labs(title="Days to Anthesis, BLOIS_2017_OPT Effect Size Estimates",subtitle="qDTA3_1"))
dev.off()
