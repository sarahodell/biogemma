#!/usr/bin/env Rscript

### Analyzing DH line genotype probabilities

library("dplyr")
library("data.table")
library("ggplot2")
library("reshape2")
library("tibble")
library('cowplot')

real_breaks=c()
actual=c()
for(c in 1:10){
  a=fread(sprintf('../stats/crossovers/Biogemma_xo_number_chr%.0f.txt',c),data.table=F)
  b=fread(sprintf('../stats/crossovers/Biogemma_breakpoints_chr%.0f.txt',c),data.table=F)
  real_breaks=rbind(real_breaks,b)
  actual=rbind(actual,a)
}
# remove unrealistically small break points
real_breaks$size=real_breaks$end - real_breaks$start
real_breaks=real_breaks[real_breaks$size>100,]

r=ggplot(real_breaks)

asum=real_breaks %>% group_by(sample) %>% count()
acsum=real_breaks %>% group_by(sample,chr) %>% count()
#asum=actual %>% group_by(sample) %>% summarize(xo_no=sum(xo_no))
founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra","OH43_inra","A654_inra",
"FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")
asum$rep=101
asum$real=T
asum=as.data.frame(asum,stringsAsFactors=F)
asum=asum[asum$sample!="EB.10H.H.00023",]
names(asum)=c('sample','n','rep','real')
ssums=c()
all_breaks=c()
for(r in 1:100){
  xo=fread(sprintf('breaktables/MAGIC_DH_Sim_rep%.0f_breaktable_v2.txt',r),data.table=F)
  xo$rep=r
  xo$size=xo$end-xo$start
  xo=xo[xo$size>100,]
  all_breaks=rbind(all_breaks,xo)
  ssum=xo %>% group_by(sample) %>% count()
  ssum$rep=r
  ssums=rbind(ssums,ssum)
}
ssums$real=F
ssums=as.data.frame(ssums,stringsAsFactors=F)
ssums=rbind(ssums,asum)


p=ggplot(ssums,aes(x=n,group=rep,color=real,fill=real)) + geom_density(alpha=0.3) +
 xlab("Crossover Counts") + ylab("Density") +
  scale_fill_discrete(name = "Type", labels = c("Simulated","Actual")) + guides(color=F)

png('xo_count_density.png')
print(p)
dev.off()

all_breaks$size=all_breaks$end - all_breaks$start

real_breaks$size=real_breaks$end - real_breaks$start
real_breaks=real_breaks[,c('sample','chr','start','end','donor1','size')]
names(real_breaks)=c('sample','chr','start','end','donor','size')
real_breaks$rep=101
real_breaks=real_breaks[,c('sample','chr','start','end','donor','rep','size')]
real_breaks$real=T
all_breaks=all_breaks[,c('sample','chr','start','end','donor','rep','size')]
all_breaks$real=F
all_breaks=rbind(all_breaks,real_breaks)

p=ggplot(all_breaks,aes(x=size,group=rep,color=real,fill=real)) + geom_density(alpha=0.3) +
 xlab("Recombination Block Size") + ylab("Density") +
  scale_fill_discrete(name = "Type", labels = c("Simulated","Actual")) + guides(color=F)

png('segment_size_density.png')
print(p)
dev.off()


p=ggplot(real_breaks,aes(x=size/1e6)) + geom_histogram(bins=100) +
 xlab("Recombination Block Size (Mb)") + ylab("Frequency") +
 geom_vline(xintercept=mean(real_breaks$size/1e6),color='red')

png('real_segment_size_hist.png')
print(p)
dev.off()

summary=ssums %>% group_by(rep) %>% summarize(mean=mean(n),var=var(n))
a=summary[summary$rep==101,]$mean

summary=summary[summary$rep!=101,]
upper.cl=quantile(summary$mean,0.95)
lower.cl=quantile(summary$mean,0.05)

dens <- density(summary$mean)
data <- tibble(x = dens$x, y = dens$y) %>%
  mutate(variable = case_when(
    (x >= lower.cl & x <= upper.cl) ~ "On",
    (x < lower.cl | x > upper.cl) ~ "Off",
    TRUE ~ NA_character_))


p2=ggplot(data,aes(x,y)) + geom_line() +
 geom_area(data=filter(data,variable=="On"),fill="grey") +
 geom_vline(xintercept=a,color='red') + xlab("Crossover Number Population Mean")
#  geom_segment(x=lower.cl,xend=upper.cl,y=1,yend=1,color='grey')

png('xo_mean_density.png')
print(p2)
dev.off()
