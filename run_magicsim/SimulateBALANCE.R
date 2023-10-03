#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
rep=as.numeric(args[[1]])

### Script to simulate 400 Double Haploid MAGIC lines
#require("devtools")
#require("roxygen2")

#install_path='/home/sodell/R/x86_64-pc-linux-gnu-library/4.2/'
#devtools::install_github('sarahodell/magicsim',lib=install_path,force=T)
#date=format(Sys.time(), "%m%d%y")
set.seed(rep)
library('magicsim',lib='/home/sodell/R/x86_64-pc-linux-gnu-library/4.2/')
library("data.table")
library("tidyverse")
library(stringr)
mapfile='../genotypes/misc/ogutmap_v4_ordered.csv'
gmap=fread(mapfile,data.table=F)

founders=c("B73_inra","A632_usa","CO255_inra","FV252_inra",
           "OH43_inra","A654_inra","FV2_inra","C103_inra",
           "EP1_inra","D105_inra","W117_inra","B96","DK63",
           "F492","ND245","VA85")

if(rep==1){
  cross_order=c(1,4,14,15,3,8,10,16,2,9,5,6,12,7,11,13)
}else{
  cross_order=sample(seq(1,16))
}
#cross_order=c(1,4,14,15,3,8,10,16,2,9,5,6,12,7,11,13)
cross_info=founders[cross_order]
#line=c(rep,cross_info)
#names(line)=c('rep',seq(1,16))
#line=data.frame(t(line),stringsAsFactors=F)
#fwrite(line,'founder_cross_info.txt',quote=F,row.names=F,sep='\t',append=T)
### Simulate chromsome 10
c=10


founder_pop=pop_init(n=16,c=10,donors=cross_info,gmap)

#Make f1s
f1s=new("Pop",nIndv=8,indvlist=vector("list",length=8))

count=1
for(i in seq(1,16,2)){
  f1s@indvlist[[count]]=make_f1(founder_pop[i],founder_pop[i+1],gmap,chroms=10)
  count=count+1
}

# 8 identical F1s


# 2nd generation
# F7 used as reference (very early)


#4-way hybrid
# X crosses, 138 meiotic events
#f2_list[[1]] HS1 x HS2
#f2_list[[2]] HS3 x HS4

f2_list=vector("list",length=4)
f2_count=200
count2=1
for(p in seq(1,8,2)){
	f2s=new("Pop",nIndv=f2_count,indvlist=vector("list",length=f2_count))
	count=1
	draw=sample(c("T","M","P"),f2_count,prob=c(0.16,0.68,0.16),replace=TRUE)
	for(i in seq(1,f2s@nIndv)){
    	f2=offspring(f1s[p],f1s[p+1],c,gmap)
    	# give the F2 a phenotyp
    	newpheno=new("Phenotype",label="ft",category=draw[i])
		f2s@indvlist[[count]]=new("IndvPheno",f2,pheno=newpheno)
    	count=count+1
  }
  f2_list[[count2]]=f2s
  count2=count2+1
}

#8-way hybrid
# harvest at least 60 ears per 8-way hybrids
# 138 crosses, 720 meiotic events
f3_list=vector("list",length=2)

#f3_count=360
count2=1
types=c("P-P","M-M","T-T","P-T","T-P")
ft_gen3_list=list(cross1=list("P-P"=list(ncrosses=12,nkpe=4),'M-M'=list(ncrosses=19,nkpe=8),"T-T"=list(ncrosses=13,nkpe=4),"P-T"=list(ncrosses=14,nkpe=4),"T-P"=list(ncrosses=10,nkpe=4),totalcrosses=68,totalkernels=348),cross2=list("P-P"=list(ncrosses=18,nkpe=4),'M-M'=list(ncrosses=23,nkpe=8),"T-T"=list(ncrosses=12,nkpe=4),"P-T"=list(ncrosses=9,nkpe=4),"T-P"=list(ncrosses=8,nkpe=4),totalcrosses=70,totalkernels=372))
for(p in seq(1,4,2)){
	if(p==1){
		crosscounts=ft_gen3_list[["cross1"]]
		#grab phenotypes of population1
	}else{
		crosscounts=ft_gen3_list[["cross2"]]
	}
	pop1=f2_list[[p]]
	pop2=f2_list[[p+1]]

	phenos1=c()
	for(n in 1:pop1@nIndv){
		phenos1=c(phenos1,pop1[n]@pheno@category)
	}
	phenos2=c()
	for(n in 1:pop2@nIndv){
		phenos2=c(phenos2,pop2[n]@pheno@category)
	}
	count=1
	f3_count=crosscounts[["totalkernels"]]
	f3s=new("Pop",nIndv=f3_count,indvlist=vector("list",length=f3_count))
	for(t in types){
		t1=strsplit(t,'-')[[1]][1]
		t2=strsplit(t,'-')[[1]][2]
		pos1=which(phenos1==t1)
		pos2=which(phenos2==t2)
		ccount=crosscounts[[t]][["ncrosses"]]
		kcount=crosscounts[[t]][["nkpe"]]
		draw1=sample(pos1,ccount,replace=T)
		draw2=sample(pos2,ccount,replace=T)
		for(c1 in 1:ccount){
			for(k1 in 1:kcount){
				newind=offspring(pop1[draw1[c1]],pop2[draw2[c1]],c,gmap)
				# offspring has 50/50 chance to be phenotype of either parent
				newdraw=sample(c(t1,t2),1)
				newpheno=new("Phenotype",label="ft",category=newdraw)
				f3s@indvlist[[count]]=new("IndvPheno",newind,pheno=newpheno)
    			count=count+1
			}	
		}
	}
 	f3_list[[count2]]=f3s
 	count2=count2+1
}

# Half of the kernels were planted in one and half were planted in the other - I don't think it makes a difference here

########16-way hybrid ##########
#288 crosses, 2200 meoitic events


# 300 genotypes samples from 1st pop and 300 genotypes sampled from 2nd pop
#f3_count=360
# assume two kernels per ear
types=c("M-M","M-P","M-T","P-M","P-P","P-T","T-M","T-P","T-T")
crosscounts=list("M-M"=172,"M-P"=17,"M-T"=24,"P-M"=22,"P-P"=15,"P-T"=1,"T-M"=25,"T-P"=1,"T-T"=11,totalcrosses=288,totalkernels=2304)

pop1=f3_list[[1]]
pop2=f3_list[[2]]

phenos1=c()
for(n in 1:pop1@nIndv){
	phenos1=c(phenos1,pop1[n]@pheno@category)
}
phenos2=c()
for(n in 1:pop2@nIndv){
	phenos2=c(phenos2,pop2[n]@pheno@category)
}
count=1
kcount=8

f4_count=crosscounts[["totalkernels"]]
f4s=new("Pop",nIndv=f4_count,indvlist=vector("list",length=f4_count))
for(type in types){
	#print(type)
	t1=strsplit(type,'-')[[1]][1]
	t2=strsplit(type,'-')[[1]][2]
	pos1=which(phenos1==t1)
	pos2=which(phenos2==t2)
	ccount=crosscounts[[type]]
	draw1=sample(pos1,ccount,replace=T)
	draw2=sample(pos2,ccount,replace=T)
	for(c1 in 1:ccount){
		#print(c1)
		for(k1 in 1:kcount){
			#print(count)
			newind=offspring(pop1[draw1[c1]],pop2[draw2[c1]],c,gmap)
			newdraw=sample(c(t1,t2),1)
			newpheno=new("Phenotype",label="ft",category=newdraw)
			f4s@indvlist[[count]]=new("IndvPheno",newind,pheno=newpheno)
			count=count+1
		}
	}
}


###### 5th generation - 1st intercross ########

# 4 blocks with 550 plants each. ~ 1600 genotypes, 801 ears harvested
types=c("M-M","M-P","M-T","P-M","P-P","T-M","T-T")
# No T-P or P-T here
crosscounts=list("M-M"=427,"M-P"=75,"M-T"=74,"P-M"=33,"P-P"=43,"T-M"=87,"T-T"=61,totalcrosses=800,totalkernels=1600)


phenos1=c()
for(n in 1:f4s@nIndv){
	phenos1=c(phenos1,f4s[n]@pheno@category)
}

kcount=2
count=1
f5_count=crosscounts[["totalkernels"]]
f5s=new("Pop",nIndv=f5_count,indvlist=vector("list",length=f5_count))
for(type in types){
	#print(type)
	t1=strsplit(type,'-')[[1]][1]
	t2=strsplit(type,'-')[[1]][2]
	pos1=which(phenos1==t1)
	pos2=which(phenos1==t2)
	ccount=crosscounts[[type]]
	if(type %in% c('M-M','T-T','P-P')){
		# If the same type sample without replacement take index i and i+2 to cross
		draw1=sample(pos1,ccount*2,replace=F)
		for(c1 in seq(1,ccount*2,2)){
			for(k1 in 1:kcount){
				#print(count)
				newind=offspring(f4s[draw1[c1]],f4s[draw1[c1+1]],c,gmap)
				newdraw=sample(c(t1,t2),1)
				newpheno=new("Phenotype",label="ft",category=newdraw)
				f5s@indvlist[[count]]=new("IndvPheno",newind,pheno=newpheno)
				count=count+1
			}
		}
	}else{
		# If different type, sample without replacement and pair i from each list
		draw1=sample(pos1,ccount,replace=F)
		draw2=sample(pos2,ccount,replace=F)
		for(c1 in 1:ccount){
			#print(c1)
			for(k1 in 1:kcount){
				#print(count)
				newind=offspring(f4s[draw1[c1]],f4s[draw2[c1]],c,gmap)
				newdraw=sample(c(t1,t2),1)
				newpheno=new("Phenotype",label="ft",category=newdraw)
				f5s@indvlist[[count]]=new("IndvPheno",newind,pheno=newpheno)
				count=count+1
			}
		}
	}
}


###### 6th generation #### 2nd intercross

types=c("M-M","M-P","M-T","P-M","P-P","T-M","T-T")
# No T-P or P-T here
crosscounts=list("M-M"=522,"M-P"=54,"M-T"=27,"P-M"=64,"P-P"=81,"T-M"=25,"T-T"=65,totalcrosses=838,totalkernels=1656)

phenos1=c()
for(n in 1:f5s@nIndv){
	phenos1=c(phenos1,f5s[n]@pheno@category)
}

kcount=2
count=1
f6_count=crosscounts[["totalkernels"]]
f6s=new("Pop",nIndv=f6_count,indvlist=vector("list",length=f6_count))
for(type in types){
	#print(type)
	t1=strsplit(type,'-')[[1]][1]
	t2=strsplit(type,'-')[[1]][2]
	pos1=which(phenos1==t1)
	pos2=which(phenos1==t2)
	ccount=crosscounts[[type]]
	if(type %in% c('M-M','T-T','P-P')){
		# If the same type sample without replacement take index i and i+2 to cross
		draw1=sample(pos1,ccount*2,replace=F)
		for(c1 in seq(1,ccount*2,2)){
			for(k1 in 1:kcount){
				#print(count)
				newind=offspring(f5s[draw1[c1]],f5s[draw1[c1+1]],c,gmap)
				newdraw=sample(c(t1,t2),1)
				newpheno=new("Phenotype",label="ft",category=newdraw)
				f6s@indvlist[[count]]=new("IndvPheno",newind,pheno=newpheno)
				count=count+1
			}
		}
	}else{
		# If different type, sample without replacement and pair i from each list
		draw1=sample(pos1,ccount,replace=F)
		draw2=sample(pos2,ccount,replace=F)
		for(c1 in 1:ccount){
			#print(c1)
			for(k1 in 1:kcount){
				#print(count)
				newind=offspring(f5s[draw1[c1]],f5s[draw2[c1]],c,gmap)
				newdraw=sample(c(t1,t2),1)
				newpheno=new("Phenotype",label="ft",category=newdraw)
				f6s@indvlist[[count]]=new("IndvPheno",newind,pheno=newpheno)
				count=count+1
			}
		}
	}
}

###### 7th generation ##### 3rd intercross

# 1 kernel kept to make double haploids

types=c("M-M","M-P","M-T","P-M","P-P","T-M","T-T")
# No T-P or P-T here
crosscounts=list("M-M"=522,"M-P"=54,"M-T"=27,"P-M"=64,"P-P"=81,"T-M"=25,"T-T"=65,totalcrosses=838,totalkernels=838)

phenos1=c()
for(n in 1:f6s@nIndv){
	phenos1=c(phenos1,f6s[n]@pheno@category)
}

kcount=1
count=1
f7_count=crosscounts[["totalkernels"]]
f7s=new("Pop",nIndv=f7_count,indvlist=vector("list",length=f7_count))
for(type in types){
	#print(type)
	t1=strsplit(type,'-')[[1]][1]
	t2=strsplit(type,'-')[[1]][2]
	pos1=which(phenos1==t1)
	pos2=which(phenos1==t2)
	ccount=crosscounts[[type]]
	if(type %in% c('M-M','T-T','P-P')){
		# If the same type sample without replacement take index i and i+2 to cross
		draw1=sample(pos1,ccount*2,replace=F)
		for(c1 in seq(1,ccount*2,2)){
			for(k1 in 1:kcount){
				#print(count)
				newind=offspring(f6s[draw1[c1]],f6s[draw1[c1+1]],c,gmap)
				newdraw=sample(c(t1,t2),1)
				newpheno=new("Phenotype",label="ft",category=newdraw)
				f7s@indvlist[[count]]=new("IndvPheno",newind,pheno=newpheno)
				count=count+1
			}
		}
	}else{
		# If different type, sample without replacement and pair i from each list
		draw1=sample(pos1,ccount,replace=F)
		draw2=sample(pos2,ccount,replace=F)
		for(c1 in 1:ccount){
			#print(c1)
			for(k1 in 1:kcount){
				#print(count)
				newind=offspring(f6s[draw1[c1]],f6s[draw2[c1]],c,gmap)
				newdraw=sample(c(t1,t2),1)
				newpheno=new("Phenotype",label="ft",category=newdraw)
				f7s@indvlist[[count]]=new("IndvPheno",newind,pheno=newpheno)
				count=count+1
			}
		}
	}
}


#synth_pop=outcross(f4s,3,10,gmap)

dh_pop=make_dh_pop(f7s,n=325,c=10,gmap)
saveRDS(dh_pop,sprintf('breaktables/MAGIC_DH_344_Sim_rep%.0f_v3.rds',rep))

#dh_pop=readRDS(sprintf('breaktables/MAGIC_DH_344_Sim_rep%.0f.rds',rep))

pop_breaks_files=c()
#Convert the dh_pop object into breakpoints fille format
for(chr in 1:10){
  pop_breaks_file=make_pop_breaktable(dh_pop,325,chr,het=F)
  pop_breaks_files=rbind(pop_breaks_files,pop_breaks_file)
}

fwrite(pop_breaks_files,sprintf('breaktables/MAGIC_DH_Sim_rep%.0f_breaktable_v3.txt',rep),row.names=F,quote=F,sep='\t')


#total_xo=c()
#for(i in seq(1,344)){total_xo=c(total_xo,dh_pop[i][1]@h1@xo_no)}

#crossovers=data.frame(ind=seq(1,344),xo_no=total_xo)
#ggplot(crossovers,aes(x=xo_no)) + geom_histogram(bins=9,fill="darkgreen",color="black") + theme_classic() +
#  ggtitle("Distribution of Crossovers in Simulation Synth3 Population") + xlab("Crossover Number") +
#  ylab("Frequency") + geom_vline(xintercept=mean(crossovers$xo_no))
