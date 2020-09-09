#!/usr/bin/env Rscript
###Run R/qtl2 on 344 16-way MAGIC DH lines

library('broman')
library('qtl2')
library('qtlcharts')

args=commandArgs(trailingOnly=TRUE)
c=as.character(args[1])
cores=as.character(args[[2]])

control_file=sprintf('genotypes/qtl2/startfiles/Biogemma_c%s.json',c)
outfile=sprintf("genotypes/probabilities/geno_probs/raw/bg%s_genoprobs_010319.rds",c)

bg<-read_cross2(control_file)
bg<-drop_nullmarkers(bg)

# Compare samples and find bad samples

percent_missing <- n_missing(bg, "ind", "prop")*100
print(round(sort(percent_missing, decreasing=TRUE)[1:19], 1))


cg <- compare_geno(bg, cores=0)
print(head(summary(cg)))


# Genotype frequencies

g <- bg$geno[[c]]
fg <- bg$founder_geno[[c]]
g <- g[,colSums(fg==0)==0]
fg <- fg[,colSums(fg==0)==0]
fgn <- colSums(fg==3)

gf_ind <- vector("list", 8)
for(i in 1:8) {
    gf_ind[[i]] <- t(apply(g[,fgn==i], 1, function(a) table(factor(a, 1:3))/sum(a != 0)))
}

pr <- calc_genoprob(bg,error_prob=0.002,cores=cores)
m <- maxmarg(pr, minprob=0.5,cores=cores)
nxo <- count_xo(m, cores=1)
totxo <- rowSums(nxo)

png(sprintf('qc/bg%s_xono.png',c))
print(iplot(seq_along(totxo)[percent_missing < 19.97],
      totxo[percent_missing < 19.97],
      chartOpts=list(xlab="Ind", ylab="Number of crossovers",
                     margin=list(left=80,top=40,right=40,bottom=40,inner=5),
                     axispos=list(xtitle=25,ytitle=50,xlabel=5,ylabel=5))))
dev.off()

png(sprintf('qc/bg%s_sample_matches.png',c))
par(mar=c(5.1,0.6,0.6, 0.6))
print(hist(cg[upper.tri(cg)], breaks=seq(0, 1, length=201),
     main="", yaxt="n", ylab="", xlab="Proportion matching genotypes"))
rug(cg[upper.tri(cg)])
dev.off()

# Genotyping error LOD scores
e <- calc_errorlod(bg, pr, cores=cores)

errors_ind <- rowSums(e[[c]]>2)/n_typed(bg)*100
png(sprintf('qc/bg%s_genotyping_errors.png',c))
lab <- paste0(names(errors_ind), " (", round(percent_missing,1), "%)")
print(iplot(seq_along(errors_ind), errors_ind, indID=lab,
      chartOpts=list(xlab="Ind", ylab="Percent genotyping errors", ylim=c(0, 10),
                     axispos=list(xtitle=25, ytitle=50, xlabel=5, ylabel=5))))
dev.off()

#Proportion missing markers
pmis_mar <- n_missing(bg, "marker", "proportion")*100
errors_mar <- colSums(e[[c]]>2)/n_typed(bg, "marker")*100
png(sprintf('qc/bg%s_missing_markers_hist.png',c))
par(mar=c(5.1,0.6,0.6, 0.6))
print(hist(pmis_mar, breaks=seq(0, 100, length=201),
     main="", yaxt="n", ylab="", xlab="Percent missing genotypes"))
rug(pmis_mar)
dev.off()

#Percent missing by error
png(sprintf('qc/bg%s_missing_markers_by_error.png',c))
print(grayplot(pmis_mar, errors_mar,
         xlab="Proportion missing", ylab="Proportion genotyping errors"))
dev.off()

#Clean data (drop markers above a error margin of 5)

gmap <- bg$gmap
pmap <- bg$pmap
bg <- drop_markers(bg, names(errors_mar)[errors_mar > 5])

prcl <- calc_genoprob(bg,error_prob=0.002,cores=cores)
pr_clean <- clean_genoprob(prcl)
#print(dim(pr[[1]]))


saveRDS(pr_clean,outfile)
