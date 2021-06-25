#!/usr/bin/env Rscript

library('data.table')
library('dplyr')
library('igraph')
library('ggplot2')
library('reshape2')

total=0
incompletes=list()
count=1

for(c in 1:10){
  allgraph=readRDS(sprintf('refinedibd/600K/bg%.0f_ibd_graph.rds',c))
  nh=dim(allgraph)[3]
  total=total+nh
  for(i in 1:nh){
    adj=as.matrix(allgraph[,,i])
    graph=graph_from_adjacency_matrix(adj,mode=c('undirected'))
    #blocks=unname(components(graph)$membership)
    #n_grps=components(graph)$no
    csize=components(graph)$csize
    comps=which(csize>2)
    dg=decompose.graph(graph)
    for(d in comps){
      t=simplify(dg[[d]])
      n=vcount(t)
      e=ecount(t)
      complete_e=((n*(n-1))/2)
      if(complete_e!=e){
        incompletes[[count]]=list(adjacency=adj,chr=c,index=i,members=names(components(t)$membership),complete_e=complete_e,edges=e,vertices=n)
        count=count+1
      }
    }
  }
}
ntotal=6929
perc_incomplete=length(incompletes)/ntotal * 100
#17.34738

uniq_haps=length(unique(lapply(incompletes,function(x) paste0(x$chr,'_',x$index))))

perc_incomplete2=uniq_haps/ntotal * 100
# 16.62578

vertices=unlist(unname(lapply(incompletes,function(x) x$vertices)))
actual_edges=unlist(unname(lapply(incompletes,function(x) x$edges)))
complete_edges=unlist(unname(lapply(incompletes,function(x) x$complete_e)))
diffs=data.frame(vertices=vertices,actual_edges=actual_edges,complete_edges=complete_edges,stringsAsFactors=F)
diffm=melt(diffs,'vertices')
diffm$vertices_f=factor(diffm$vertices,levels=sort(unique(diffm$vertices)))
diffm$variable=as.character(diffm$variable)
#p=ggplot(diffm,aes(x=vertices_f,y=value,group=vertices_f)) + geom_boxplot(aes(group=variable,fill=variable))

t= diffm %>% group_by(vertices,variable) %>% summarize(mean=mean(value),sd=sd(value),n=length(value))
subt=t[t$variable=="actual_edges",]

p=ggplot(t,aes(x=vertices,y=mean,group=vertices)) +
 geom_point(aes(group=variable,color=variable),size=4) +
  geom_errorbar(data=subt,aes(ymin = mean-(2*sd), ymax = mean+(2*sd),color=variable), position=position_dodge(-0.9),width = 0.2) +
  xlab("Haplotype Members") + ylab("Number of Pairwise IBD Connections") +
  scale_y_continuous(limits=c(0,40),breaks=seq(0,40,5)) + scale_x_continuous(breaks=c(3,4,5,6,7,8,9)) +
  scale_color_discrete(name = "Group", labels = c("Actual","Complete")) +
  theme_classic() +
  theme(panel.grid.major.y=element_line(),
  panel.grid.major.x=element_line(),
  legend.position = c(0.1, 0.9))


png('incomplete_graph_deviation.png',width=600,height=400)
print(p)
dev.off()
