# Break up haplotype probability files so that each chromosome has h files where
# h is the number of unique haplotype groups and hap[[i]] is an nxp matrix of probabilities
# that individual n has haplotype i at marker p
# Markers that are more than 0.95 correlated are dropped and an R object with dropped markers is saved

args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])

hap=readRDS(sprintf('bg%s_haplotype_probs_010319.rds',chr))
for(h in 2:length(hap)){
   hsub=hap[[h]]
   if(is.null(dim(hsub))==F){
	regrp=c()
	for(i in 1:h){regrp[[i]]=hsub[,,i]}
	saveRDS(regrp,sprintf('bg%s_haplogroup%.0f_reformat.rds',chr,h))
   }
}



