# Break up haplotype probability files so that each chromosome has h files where
# h is the number of unique haplotype groups and hap[[i]] is an nxp matrix of probabilities
# that individual n has haplotype i at marker p
# Markers that are more than 0.95 correlated are dropped and an R object with dropped markers is saved

args=commandArgs(trailingOnly=T)
chr=as.character(args[[1]])
#date=format(Sys.time(), "%m%d%y")

cutoff=0.95

hap=readRDS(sprintf('bg%s_haplotype_probs_020320.rds',chr))
for(h in 2:length(hap)){
   hsub=hap[[h]]
   if(is.null(dim(hsub))==F){
	regrp=c()
	for(i in 1:h){regrp[[i]]=hsub[,,i]}
	# Filter Correlated Markers in Haplotype Probability
	m_names=names(regrp[[1]][1,])
	dropped=list()
	count=1
	dropped[[count]]=list(marker=c(m_names[1]),linked=c())
	
	keep=c(1)
	start=c()
	for(x in 1:h){start=c(start,as.vector(regrp[[x]][,1]))}
	for(i in 2:dim(regrp[[1]])[2]){
  	      ind_max = c()
	      for(x in 1:h){ind_max=c(ind_max,as.vector(regrp[[x]][,i]))}
  	      if(cor(start,ind_max)<cutoff){
		keep=c(keep,i)
		start=ind_max
		start_ind=i
		count=count+1
		dropped[[count]]=list(marker=c(m_names[start_ind]),linked=c())
  	      }
	      else{
	        dropped[[count]]$linked=c(dropped[[count]]$linked,m_names[i])
	      }
	}
	saveRDS(dropped,sprintf('bg%s_haplogroup%0.f_dropped_markers.rds',chr,h))
	filtered_hap=c()
	for(k in 1:h){filtered_hap[[k]]=regrp[[k]][,keep]}
	saveRDS(filtered_hap,sprintf('bg%s_filtered_haplogroup%.0f_probs_2.rds',chr,h))
	sprintf("Keeping %0.f markers in haplogroup %.0f on chromosome %s with correlation filter of %.2f",length(keep),h,chr,cutoff)
   }
}