#!/usr/bin/env Rscript

library('data.table')
library('qtl2')
library('ggplot2')

founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
pr=readRDS('MAGICSim_121718_chr10_genoprobs.rds')

bg=read_cross2('MAGICSim_121718_c10.json')
pmap=fread('Physical_map_c10.csv',data.table=F)
dimnames(pr[[1]])[[2]]=founders

fgeno=t(bg$founder_geno[[1]])
sdp=calc_sdp(fgeno)
pos=pmap[pmap$marker %in% names(sdp),]$pos
snpinfo=data.frame(chr=c("10"),pos=pos,sdp=sdp,snp=names(sdp),stringsAsFactors=F)
snpinfo=index_snps(bg$pmap,snpinfo)
snp_prbs=genoprob_to_snpprob(pr,snpinfo)
saveRDS(snp_prbs,'bg10_fulldata_SNPprobs_010819.rds')

rsq<-function(x,y) cor(x,y,method=c("pearson")) ^2

markers=dimnames(snp_prbs[[1]])[[3]]
included=bg$geno[[1]][,dimnames(bg$geno[[1]])[[2]] %in% markers]

r1=0
w1=0
r3=0
w3=0

r_squared=c()

for(i in 1:400){
   p_ind=unlist(unname(snp_prbs[[1]][i,,]))
   p_snps=sapply(seq(1,dim(p_ind)[2]),function(x) ifelse(max(p_ind[,x])==p_ind[1,x],1,3))
   a_snps=unlist(unname(included[i,]))
   for(x in 1:dim(p_ind)[2]){
      if(p_snps[x]==1 & a_snps[x]==1){
         r1=r1+1
      }
      else if(p_snps[x]==1 & a_snps[x]==3){
        w3=w3+1
      }
      else if(p_snps[x]==3 & a_snps[x]==1){
        w1=w1+1
      }
      else{
        r3=r3+1
      }
   }
   r_squared=c(r_squared,rsq(a_snps,p_snps))
}

size=c(r1,w1,w3,r3)
result=data.frame(actual=as.factor(c(1,1,3,3)),predicted=as.factor(c(1,3,1,3)),size=size)
fwrite(result,'imputation_accuracy_010819.txt',row.names=F,quote=F,sep='\t')

r2=mean(r_squared)

png('MAGICSim_Accuracy_plot_010819.png',width=460,height=480)
print(ggplot(result,aes(x=actual,y=predicted,size=size)) + geom_point() + scale_x_discrete("Actual", labels=c(sprintf('Reference, n=%.0f',r1+w3),sprintf('Alternate, n=%.0f',r3+w1))) + scale_y_discrete("Predicted",labels=c('Reference','Alternate')) + labs(title="Imputation Accuracy of 400 Simulated MAGIC lines",subtitle=sprintf("Mean R-Squared: %.3f , n=%.0f , %.3f %% Correct",r2,sum(size),(r1+r3)/sum(size))) + theme_bw() + guides(size=F) + theme(text=element_text(family="Georgia")))
dev.off()

#actual_mat=c()

#for(i in seq(1,dim(pr[[1]][3]))){
#   p=array(pr[[1]][,,i],dim=c(400,16))
#   marker=dimnames(pr[[1]])[[3]][i]
#   pos=pmap[pmap$marker==marker,]$pos
#   block=actual[(actual$start<pos) & (actual$end > pos),]$donor1
#   a=sapply(seq(1,10),function(x) founders==block[x]
#   r_squared=abind(r_squared,rsq(t(p),a))
#   a=array(t(a),dim=c(10,16,1))
#   actual_mat=abind(actual_mat,a)
#}

#fwrite(actual_mat,'Actual_MAGICSim_121718_geno.csv',row.names=F,quote=F,sep=',')
#saveRDS(r_squared,'MAGICSim_121718_rsquared.rds')
