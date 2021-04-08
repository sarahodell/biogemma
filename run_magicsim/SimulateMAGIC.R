#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
rep=as.numeric(args[[1]])

### Script to simulate 400 Double Haploid MAGIC lines
#require("devtools")
#require("roxygen2")

#install_path='/home/sodell/R/x86_64-pc-linux-gnu-library/3.6/'
#devtools::install_github('sarahodell/magicsim',lib=install_path,force=T)
#date=format(Sys.time(), "%m%d%y")
set.seed(rep)
library('magicsim')
library("data.table")
library("tidyverse")

mapfile='../genotypes/misc/ogutmap_v4_ordered.csv'
gmap=fread(mapfile,data.table=F)

#founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra",
#           "OH43_inra","A654_inra","FV2_inra","C103_inra",
#           "EP1_inra","D105_inra","W117_inra","B96","DK63",
#           "F492","ND245","VA85")
#cross_order=c(4,2,14,15,3,8,10,16,1,9,5,6,12,7,11,13)
#cross_info=founders[cross_order]

### Simulate chromsome 10
#c=10


#founder_pop=magicsim::pop_init(n=16,c=10,donors=cross_info,gmap)

#Make f1s
#f1s=new("Pop",nIndv=8,indvlist=vector("list",length=8))

#count=1
#for(i in seq(1,16,2)){
#  f1s@indvlist[[count]]=make_f1(founder_pop[i],founder_pop[i+1],gmap,chroms=10)
#  count=count+1
#}

#magic_pop=make_magic_pop(start_pop=f1s,c=10,g_map=gmap,popsize=2000)

#synth_pop=outcross(magic_pop,3,10,gmap)

#dh_pop=make_dh_pop(synth_pop,n=400,c=10,gmap)
#saveRDS(dh_pop,sprintf('breaktables/MAGIC_DH_400_Sim_rep%.0f.rds',rep))

dh_pop=readRDS(sprintf('breaktables/MAGIC_DH_400_Sim_rep%.0f.rds',rep))

pop_breaks_files=c()
#Convert the dh_pop object into breakpoints fille format
for(chr in 1:10){
  pop_breaks_file=make_pop_breaktable(dh_pop,400,chr,het=F)
  pop_breaks_files=rbind(pop_breaks_files,pop_breaks_file)
}

fwrite(pop_breaks_files,sprintf('breaktables/MAGIC_DH_Sim_rep%.0f_breaktable.txt',rep),row.names=F,quote=F,sep='\t')
