#!/usr/bin/env Rscript

library('ggplot2')
library('data.table')
library('dplyr')

#actual=fread('../GridLMM/pgs/Biogemma_FT_PGS.txt',data.table=F)
actual=fread('../GridLMM/pgs/Biogemma_FT_PGS.txt',data.table=F)
actual=fread('../GridLMM/pgs/Castelletti_DTA_FT_PGS.txt',data.table=F)
actual=fread('../GridLMM/pgs/Castelletti_markers_BG_DTA_FT_PGS.txt',data.table=F)
actual$rep=101
actual$real=T
#actual
pgs=c()
for(r in 1:100){
  p=fread(sprintf('pgs/Rep%.0f_Castelletti_markers_BG_DTA_FT_PGS.txt',r),data.table=F)
  p$rep=r
  pgs=rbind(pgs,p)
}

pgs$real=F
pgs=rbind(pgs,actual[,c('dta','ind','rep','real')])

p=ggplot(pgs,aes(x=dta,group=rep,color=real,fill=real)) + geom_density(alpha=0.3) +
 xlab("Male Flowering Time Polygenic Score") + ylab("Density") +
  scale_fill_discrete(name = "Type", labels = c("Simulated","Actual")) + guides(color=F)

png('castelletti_markers_bg_ft_pgs_density.png')
print(p)
dev.off()

summary=pgs %>% group_by(rep) %>% summarize(mean=mean(dta),var=var(dta))
a=summary[summary$rep==101,]$mean

ecdf(summary$mean)(a) #77% of simulated pop means were greater than actual
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
 geom_vline(xintercept=a,color='red') + xlab("DTS Population Mean")
#  geom_segment(x=lower.cl,xend=upper.cl,y=1,yend=1,color='grey')

png('castelletti_marker_bg_ft_pgs_mean_density.png')
print(p2)
dev.off()

summary=pgs %>% group_by(rep) %>% summarize(mean=mean(dta),var=var(dta))
a=summary[summary$rep==101,]$var

ecdf(summary$var)(a) #100% of simulated pop means were greater than actual
summary=summary[summary$rep!=101,]
upper.cl=quantile(summary$var,0.95)
lower.cl=quantile(summary$var,0.05)

dens <- density(summary$var)
data <- tibble(x = dens$x, y = dens$y) %>%
  mutate(variable = case_when(
    (x >= lower.cl & x <= upper.cl) ~ "On",
    (x < lower.cl | x > upper.cl) ~ "Off",
    TRUE ~ NA_character_))

p2=ggplot(data,aes(x,y)) + geom_line() +
 geom_area(data=filter(data,variable=="On"),fill="grey") +
 geom_vline(xintercept=a,color='red') + xlab("DTS Population Variance")
#  geom_segment(x=lower.cl,xend=upper.cl,y=1,yend=1,color='grey')

png('castelletti_marker_bg_ft_dta_pgs_var_density.png')
print(p2)
dev.off()
