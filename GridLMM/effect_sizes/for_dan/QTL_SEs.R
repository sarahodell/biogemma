
library('data.table')
library('lme4qtl')

founder_probs = readRDS('male_flowering_d6_BLOIS_2017_OPT_qDTA3_1_founder_probs.rds')
K = fread('K_matrix_chr3.txt',data.table = F,h=T)
rownames(K) = K[,1]
K = as.matrix(K[,-1])
pheno = fread('male_flowering_d6_BLOIS_2017_OPT_phenotypes.csv',data.table = F)
pheno$Genotype_code = gsub('-','.',pheno$Genotype_code,fixed=T)

X = do.call(cbind,founder_probs)
colnames(X) = paste0('F',1:ncol(X))
rownames(X) = names(founder_probs[[1]])
pheno = subset(pheno,Genotype_code %in% rownames(K))
X = X[rownames(X) %in% pheno$Genotype_code,]
# pheno = cbind(pheno[match(rownames(X),pheno$Genotype_code),1:3],X_diff_mean)
# model with no intercept
m0 = relmatLmer(male_flowering_d6 ~ 0 + X + (1|Genotype_code),data=pheno,relmat = list(Genotype_code=K))
summary(m0)$coef

# model with intercept, dropping founder 16
X_diff_mean = X %*% contr.sum(ncol(X))
colnames(X_diff_mean) = paste0('F',1:ncol(X_diff_mean))
m1 = relmatLmer(male_flowering_d6 ~ X_diff_mean + (1|Genotype_code),data=pheno,relmat = list(Genotype_code=K))
coefs1 = data.frame(summary(m1)$coef)[-1,]

# estimating SE for F16
coefs1$Precision = 1/coefs1[,2]^2
coefs1$nobs = colSums(X)[-16]
m1_slope = lm(Precision~nobs,coefs1[-1,])
with(coefs1,plot(Precision~nobs))
abline(m1_slope)
F16 = -sum(coefs1[,1])
F16_SE = 1/sqrt(predict(m1_slope,newdata = list(nobs = colSums(X)[16])))

# model with intercept, dropping founders 1 and 16
X_drop1_diff_mean = X[,-1] %*% contr.sum(ncol(X[,-1]))
colnames(X_drop1_diff_mean) = paste0('F',1+1:ncol(X_drop1_diff_mean))
m2 = relmatLmer(male_flowering_d6 ~ X_drop1_diff_mean + (1|Genotype_code),data=pheno,relmat = list(Genotype_code=K))
summary(m2)$coef
coefs2 = data.frame(summary(m2)$coef)[-1,]

# estimating SE for F16
coefs2$Precision = 1/coefs2[,2]^2
coefs2$nobs = colSums(X)[-c(1,16)]
m2_slope = lm(Precision~nobs,coefs2)
with(coefs2,plot(Precision~nobs))
abline(m2_slope)
F16 = -sum(coefs2[,1])
F16_SE = 1/sqrt(predict(m2_slope,newdata = list(nobs = colSums(X)[2])))
