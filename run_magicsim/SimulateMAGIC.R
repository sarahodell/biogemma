### Script to simulate 400 Double Haploid MAGIC lines
require("devtools")
require("roxygen2")

devtools::install_github('sarahodell/magicsim')
date=format(Sys.time(), "%m%d%y")

library('magicsim')
library("data.table")
library("tidyverse")

ogutmap=fread('ogutmap_v4_ordered.csv',data.table=F)

founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra",
           "OH43_inra","A654_inra","FV2_inra","C103_inra",
           "EP1_inra","D105_inra","W117_inra","B96","DK63",
           "F492","ND245","VA85")
cross_order=c(2,4,14,15,3,8,10,16,1,9,5,6,12,7,11,13)
cross_info=founders[cross_order]

### Simulate chromsome 10
c=10
### Initialize founder chromosomes for the 16 founders
founder_list=list()
for(p in seq(1,length(cross_info))){
  founder_list[[p]]=list()
  founder_list[[p]][[c]]=magicsim::chrom_init(c,ogutmap,h1_donor=cross_info[p],h2_donor=cross_info[p])
}
### Create F1s from the first round of crossing
f1s=list()
count=1
for(i in seq(1,length(founder_list),2)){
  f1=list()
  for(x in seq(1,10)){f1[[x]]=make_f1(founder_list[[i]][[x]],founder_list[[i+1]][[x]],ogutmap,x)}
  f1s[[count]]=f1
  count=count+1
}

#create synthetic population of 2000 lines from MAGIC founder funnel
synth_pop=list()
for(x in seq(1,2000)){synth_pop[[x]]=make_magic(f1s,10)}

# Random outcrossing within the synthetic population for 3 generations
synth3=outcross(synth_pop,3,10)

### Simulate Double Haploids
#Randomly select 400 lines from the synthetic population and
# then choose one of the chromosomes at random to duplicate
draw=sample(seq(1,2000),400,replace=F)
dh_pop=list()
for(d in seq(1,400)){
  line=draw[d]
  if(runif(1)>0.5){
    dh_pop[[d]]=list(chr=10,h1=synth3[[line]][[10]]$h1,h2=synth3[[line]][[10]]$h1)
  }
  else{
    dh_pop[[d]]=list(chr=10,h1=synth3[[line]][[10]]$h2,h2=synth3[[line]][[10]]$h2)
  }
}

saveRDS(dh_pop,sprintf('MAGIC_DHSim_%s.rds',date))



no=sapply(seq(1,400),function(x) length(dh_pop[[x]]$h1$donors))
#add the crossover count as xo_no to dh_pop
for(i in seq(1,400)){dh_pop[[i]]$xo_no=count_breakpoints(dh_pop[[i]],hom=T)}
#what is the distribution of xo_no in the population?
total_xo=c()
for(i in seq(1,400)){total_xo=c(total_xo,dh_pop[[i]]$xo_no)}

crossovers=data.frame(ind=seq(1,400),xo_no=total_xo)
ggplot(crossovers,aes(x=total_xo)) + geom_histogram(bins=9,fill="darkgreen",color="black") + theme_classic() +
  ggtitle("Distribution of Crossovers in Simulation DH Population") + xlab("Crossover Number") + 
  ylab("Frequency")

#Convert the dh_pop object into breakpoints fille format

#This should be its own function
breaks_file=c()
for(i in seq(1,400)){
  sample=sprintf('Sim%.0f',i)
  chr=c
  ind=dh_pop[[i]]$h1
  first=ind$donors[1]
  start=0
  if(length(ind$breakpoints)==1){
    line=c(sample,chr,start,ind$breakpoints,first,first)
    breaks_file=rbind(breaks_file,line)
  }
  else if(length(unique(ind$donors))==1){
    end_index=length(ind$breakpoints)
    line=c(sample,chr,start,ind$breakpoints[end_index],first,first)
    breaks_file=rbind(breaks_file,line)
  }
  else{
    for(j in seq(1,length(ind$breakpoints))){
      if(ind$donors[j]!=first){
        line=c(sample,chr,start,ind$breakpoints[j],first,first)
        breaks_file=rbind(breaks_file,line)
        first=ind$donors[j]
        start=ind$breakpoints[j]+1
      }
    }
  }
}
breaks_file=as.data.frame(breaks_file)
names(breaks_file)=c('sample','chr','start','end','donor1','donor2')

fwrite(breaks_file,sprintf('MAGIC_DHSim_%s_breakpoints.txt',date),row.names =F,quote=F,sep='\t')

