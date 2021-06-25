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

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra",
           "OH43_inra","A654_inra","FV2_inra","C103_inra",
           "EP1_inra","D105_inra","W117_inra","B96","DK63",
           "F492","ND245","VA85")

#if(rep==1){
#  cross_order=c(1,4,14,15,3,8,10,16,2,9,5,6,12,7,11,13)
#}else{
#  cross_order=sample(seq(1,16))
#}
cross_order=c(1,4,14,15,3,8,10,16,2,9,5,6,12,7,11,13)
cross_info=founders[cross_order]
#line=c(rep,cross_info)
#names(line)=c('rep',seq(1,16))
#line=data.frame(t(line),stringsAsFactors=F)
#fwrite(line,'founder_cross_info.txt',quote=F,row.names=F,sep='\t',append=T)
### Simulate chromsome 10
c=10


founder_pop=magicsim::pop_init(n=16,c=10,donors=cross_info,gmap)

#Make f1s
f1s=new("Pop",nIndv=8,indvlist=vector("list",length=8))

count=1
for(i in seq(1,16,2)){
  f1s@indvlist[[count]]=make_f1(founder_pop[i],founder_pop[i+1],gmap,chroms=10)
  count=count+1
}

f2_list=vector("list",length=4)
#4-way hybrid
# X crosses, 138 meiotic events
f2_count=70
count2=1
for(p in seq(1,8,2)){
  f2s=new("Pop",nIndv=f2_count,indvlist=vector("list",length=f2_count))
  count=1
  for(i in seq(1,f2s@nIndv)){
    f2s@indvlist[[count]]=offspring(f1s[p],f1s[p+1],c,gmap)
    count=count+1
  }
  f2_list[[count2]]=f2s
  count2=count2+1
}

#8-way hybrid
# harvest at least 60 ears per 8-way hybrids
# 138 crosses, 720 meiotic events
f3_list=vector("list",length=2)

f3_count=360
count=1
for(p in seq(1,4,2)){
  draw1=sample(seq(1,f2_count),replace=T,f3_count)
  draw2=sample(seq(1,f2_count),replace=T,f3_count)
  f3s=new("Pop",nIndv=f3_count,indvlist=vector("list",length=f3_count))
  for(i in seq(1,f3s@nIndv)){
    f3s@indvlist[[i]]=offspring(f2_list[[p]][draw1[i]],f2_list[[p+1]][draw2[i]],c,gmap)
  }
  f3_list[[count]]=f3s
  count=count+1
}

#16-way hybrid
#288 crosses, 2200 meoitic events
f4_count=2200

draw1=sample(seq(1,f3_count),replace=T,f4_count)
draw2=sample(seq(1,f3_count),replace=T,f4_count)
f4s=new("Pop",nIndv=f4_count,indvlist=vector("list",length=f4_count))
for(i in seq(1,f4s@nIndv)){
  f4s@indvlist[[i]]=offspring(f3_list[[1]][draw1[i]],f3_list[[2]][draw2[i]],c,gmap)
}


synth_pop=outcross(f4s,3,10,gmap)

dh_pop=make_dh_pop(synth_pop,n=344,c=10,gmap)
saveRDS(dh_pop,sprintf('breaktables/MAGIC_DH_344_Sim_rep%.0f_v2.rds',rep))

#dh_pop=readRDS(sprintf('breaktables/MAGIC_DH_344_Sim_rep%.0f.rds',rep))

pop_breaks_files=c()
#Convert the dh_pop object into breakpoints fille format
for(chr in 1:10){
  pop_breaks_file=make_pop_breaktable(dh_pop,344,chr,het=F)
  pop_breaks_files=rbind(pop_breaks_files,pop_breaks_file)
}

fwrite(pop_breaks_files,sprintf('breaktables/MAGIC_DH_Sim_rep%.0f_breaktable_v2.txt',rep),row.names=F,quote=F,sep='\t')


#total_xo=c()
#for(i in seq(1,344)){total_xo=c(total_xo,dh_pop[i][1]@h1@xo_no)}

#crossovers=data.frame(ind=seq(1,344),xo_no=total_xo)
#ggplot(crossovers,aes(x=xo_no)) + geom_histogram(bins=9,fill="darkgreen",color="black") + theme_classic() +
#  ggtitle("Distribution of Crossovers in Simulation Synth3 Population") + xlab("Crossover Number") +
#  ylab("Frequency") + geom_vline(xintercept=mean(crossovers$xo_no))
