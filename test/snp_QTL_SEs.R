#!/usr/bin/env Rscript


library('data.table')
library('lme4')
library('lme4qtl')
library('GridLMM')

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra",
"C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

K = fread('test_K_matrix.txt',data.table = F,h=T)
rownames(K) = names(K)[-1]
K = as.matrix(K[,-1])

data_blup = fread('test_data_blups.txt',data.table=F)
rownames(data_blup)=data_blup$ID

X=fread('test_X_table.txt',data.table=F)
rownames(X)=X[,1]
X=as.matrix(X[,-1,drop=F])

new_X=fread('test_GridLMM_X.txt',data.table=F)
rownames(new_X)=new_X[,1]
new_X=as.matrix(new_X[,-1,drop=F])

m0 = relmatLmer(y ~ 0 + X + (1|ID),data=data_blup,relmat = list(ID=K))
se=as.data.frame(summary(m0)$coef)
names(se)=c('value','se','tvalue')
se$value[-1]=se$value[-1]+se$value[1]

cores=1
vars = as.data.frame(VarCorr(m0))
h2_start = vars$vcov[1]/sum(vars$vcov)
names(h2_start) = 'ID'
gwas = GridLMM_GWAS(
  formula = y~1 + (1|ID),
  test_formula = ~1,
  reduced_formula = ~0,
  data = data_blup,
  weights = NULL,
  X = new_X,
  X_ID = 'ID',
  h2_start = h2_start,
  h2_step=1,
  max_steps = 100,
  relmat = list(ID=K),
  centerX = FALSE,
  scaleX = FALSE,
  fillNAX = FALSE,
  method = 'REML',
  mc.cores = cores,
  verbose = TRUE
)

betas = unlist(gwas$results[1,grep('beta',colnames(gwas$results))])
betas[-1] = betas[-1]+betas[1]
corr=(betas['beta.2']-se$value[2])
