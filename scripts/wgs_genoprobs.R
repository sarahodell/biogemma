#!/usr/bin/env Rscript
### Creating Matrices of Genotype Probabilities for WGS data
library('data.table')
library('dplyr')
library('tibble')

founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

args=commandArgs(trailingOnly=TRUE)
c = as.character(args[1])
start = args[2]
end = args[3]
wgsfull=fread('biogemma/WGS_positions.txt',data.table=F,header=F)
names(wgsfull)=c('chr','pos')
#for each chromosome
#read in files
pr=readRDS(sprintf('Biogemma_082318/600K_DHgenoprobs/bg%s_0823_genoprobs.rds',c))
wgs=wgsfull[wgsfull$chr==c,]
wgs=wgs %>% mutate(marker=paste0('S',chr,'_',pos))
wgs=wgs[,c('marker','chr','pos')]
pmap=fread(sprintf('Biogemma_082318/Biogemma_pmap_c%s.csv',c),data.table=F)
pmap$pos=pmap$pos*1e6
positions=rbind(wgs,pmap)
positions=positions[order(positions$pos),]
rownames(positions)=seq(1,nrow(positions))
#create list
prob_list=list(chr=c,names=names(pr[[1]][,1,1]),positions=positions)
#for each DH line
for(i in seq(start,end)){
    print(i)
    new_ind=c()
    ind=t(pr[[1]][i,,])
    name=names(pr[[1]][start:end,1,1])[i]
    ind=as.data.frame(ind)
    ind<-rownames_to_column(ind,"marker")
    names(ind)=c("marker",founders)
    ind <- merge(ind,pmap,by.x='marker',by.y='marker')
    ind=ind[order(ind$pos),]
    rownames(ind)=seq(1,nrow(ind))
    ind=ind[,c('marker','chr','pos',founders)]
    start_pos=0
    max_pos=max(ind$pos)
    for(p in 1:nrow(ind)){
	current_line=ind[p,]
	current_pos=ind$pos[p]
	below=wgs[wgs$pos<current_pos & wgs$pos>start_pos,]
	if(start_pos==0){
            t=matrix(unlist(unname(current_line[4:19])),nrow=nrow(below),ncol=16,byrow=T)
	    new_ind=rbind(new_ind,t)
        }
	else if(start_pos==max_pos){
	    below=wgs[wgs$pos>start_pos,]
	    t=matrix(unlist(unname(start_line[4:19])),nrow=nrow(below),ncol=16,byrow=T)
	    new_ind=rbind(new_ind,t)
	}
	else if(dim(below)[1]==0){
            t = unlist(unname(start_line[4:19]))
	    new_ind=rbind(new_ind,t)
	}
	else{
            interp=c()
	    for(f in founders){
	        get_gp=approxfun(x=c(start_pos,current_pos),y=c(start_line[f],current_line[f]))
	        f_probs=get_gp(below$pos)
	        interp=cbind(interp,f_probs)
	    }
            t=unlist(unname(start_line[4:19]))
	    interp=unlist(interp)
	    new_ind=rbind(new_ind,t)
	    new_ind=rbind(new_ind,interp)
        }
        start_pos=current_pos
        start_line=current_line
    }
    colnames(new_ind)=c(founders)
    #rownames(new_ind)=positions$marker
    new_ind=t(new_ind)
    prob_list[[name]]=new_ind
    saveRDS(prob_list,file=sprintf('bg%s_WGS_geno_probs.rds',c))
}				      