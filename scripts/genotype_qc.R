#!/usr/bin/env Rscript

library('broman')
library('qtl2')
library('qtlcharts')
library('data.table')


print("Reading in Crossfiles")
bg=c()
for(i in 1:10){
  if(i==1){
    bg=read_cross2(sprintf('genotypes/qtl2/startfiles/Biogemma_c%.0f.json',i))
  }else{
    tmp=read_cross2(sprintf('genotypes/qtl2/startfiles/Biogemma_c%.0f.json',i))
    bg$geno[[i]]=tmp$geno[[1]]
    bg$founder_geno[[i]]=tmp$founder_geno[[1]]
    bg$pmap[[i]]=tmp$pmap[[1]]
    bg$gmap[[i]]=tmp$gmap[[1]]
  }
}

print("Droppin Null markers")
bg=drop_nullmarkers(bg)

for(chr in seq_along(bg$founder_geno)) {
    fg <- bg$founder_geno[[chr]]
    g <- bg$geno[[chr]]
    f1 <- colSums(fg==1)/colSums(fg != 0)

    fg[fg==0] <- NA
    g[g==0] <- NA

    fg[,f1 < 0.5] <- 4 - fg[,f1 < 0.5]
    g[,f1 < 0.5]  <- 4 - g[,f1 < 0.5]

    fg[is.na(fg)] <- 0
    g[is.na(g)] <- 0

    bg$founder_geno[[chr]] <- fg
    bg$geno[[chr]] <- g
}

print("Percent Missing")
percent_missing <- n_missing(bg, "ind", "prop")*100

png('qc/percent_missing.png')
labels <- paste0(names(percent_missing), " (", round(percent_missing), "%)")
print(iplot(seq_along(percent_missing), percent_missing, indID=labels,
      chartOpts=list(xlab="Sample", ylab="Percent missing genotype data",
                     ylim=c(0, 100))))
dev.off()

print(summary(percent_missing))
#fwrite(summary(percent_missing),'qc/summary.txt')

print("Comparing genos")

cg <- compare_geno(bg, cores=1)
print(summary(cg))
#fwrite(summary(cg),'qc/summary.txt',append=T)

png('qc/matching_genos.png')
par(mar=c(5.1,0.6,0.6, 0.6))
print(hist(cg[upper.tri(cg)], breaks=seq(0, 1, length=201),
     main="", yaxt="n", ylab="", xlab="Proportion matching genotypes"))
rug(cg[upper.tri(cg)])
dev.off()

percent_missing <- n_missing(bg, "ind", "prop")*100
round(sort(percent_missing, decreasing=TRUE)[1:19], 1)

#Genotype frequencies
g <- do.call("cbind", bg$geno[1:10])
fg <- do.call("cbind", bg$founder_geno[1:10])
g <- g[,colSums(fg==0)==0]
fg <- fg[,colSums(fg==0)==0]
fgn <- colSums(fg==3)

print("Calc_genoprob")

pr <- calc_genoprob(bg, error_prob=0.002, cores=4)
m <- maxmarg(pr, minprob=0.5,cores=4)
nxo <- count_xo(m, cores=4)
totxo <- rowSums(nxo)

png('qc/no_xo.png')
print(iplot(seq_along(totxo)[percent_missing < 19.97],
      totxo[percent_missing < 19.97],
      chartOpts=list(xlab="Ind", ylab="Number of crossovers",
                     margin=list(left=80,top=40,right=40,bottom=40,inner=5),
                     axispos=list(xtitle=25,ytitle=50,xlabel=5,ylabel=5))))
dev.off()

tmp <- cbind(percent_missing=round(percent_missing), total_xo=totxo)[percent_missing >= 19.97,]

print(tmp[order(tmp[,1]),])
#fwrite(tmp[order(tmp[,1]),],'qc/summary.txt')

# Genotyping error LOD scores
e <- calc_errorlod(bg, pr, cores=4)

e <- do.call("cbind", e)
errors_ind <- rowSums(e>2)/n_typed(bg)*100

png('qc/genotype_errors.png')
lab <- paste0(names(errors_ind), " (", myround(percent_missing,1), "%)")
print(iplot(seq_along(errors_ind), errors_ind, indID=lab,
      chartOpts=list(xlab="Sample", ylab="Percent genotyping errors", ylim=c(0, 4.1),
                     axispos=list(xtitle=25, ytitle=50, xlabel=5, ylabel=5))))
dev.off()

# Bad markers
n <- n_ind(bg)
bg <- bg[percent_missing < 19.97,]
bg( n_ind(bg) == n-9 )

# update other stuff
e <- e[ind_ids(bg),]
g <- g[ind_ids(bg),]

pmis_mar <- n_missing(bg, "marker", "proportion")*100

png('qc/missing_data.png')
par(mar=c(5.1,0.6,0.6, 0.6))
print(hist(pmis_mar, breaks=seq(0, 100, length=201),
     main="", yaxt="n", ylab="", xlab="Percent missing genotypes"))
rug(pmis_mar)
dev.off()


# Genotyping erros
errors_mar <- colSums(e>2)/n_typed(bg, "marker")*100

png('qc/errors_by_missing.png')
print(grayplot(pmis_mar, errors_mar,
         xlab="Proportion missing", ylab="Proportion genotyping errors"))
dev.off()

gmap <- bg$gmap
pmap <- bg$pmap
bg <- drop_markers(bg, names(errors_mar)[errors_mar > 5])

prcl <- calc_genoprob(bg,gmap,error_prob=0.002, cores=4)
prcl <- prcl[ind_ids(bg),]
pr <- pr[ind_ids(bg),]

prdiff <- vector("list", length(pr))
for(i in seq_along(prdiff)) prdiff[[i]] <- apply(abs(pr[[i]] - prcl[[i]]), c(1,3), sum)
names(prdiff) <- names(pr)

png('qc/EB.09S.H.00002_c10.png')
print(plot_genoprobcomp(pr, prcl, pmap, ind="EB.09S.H.00002", chr="10", threshold=0.25))
dev.off()

png('qc/EB.09S.H.00002_c10_zoom.png')
par(mfrow=c(2,1), mar=c(3.1, 3.1, 2.6, 0.6))
grayplot(pmap[["10"]], pr[["10"]]["EB.09S.H.00002","A632_usaA632_usa",], type="l",
         ylim=c(0,1), col="slateblue",
         xlab="Chr 10 position (Mbp)", ylab="Genotype probability",
         main="Line EB.09S.H.00002\n(before cleaning)",
         mgp.x=c(1.4,0.3,0), mgp.y=c(1.9,0.3,0), lwd=2)
lines(pmap[["10"]], pr[["10"]]["EB.09S.H.00002","B73_inraB73_inra",], col="violetred", lwd=2)
legend("topright", lwd=2, col=c("slateblue", "violetred"),
       c("A632_usa", "B73_inra"), bg="gray92")

grayplot(pmap[["10"]], prcl[["10"]]["EB.09S.H.00002","A632_usaA632_usa",], type="l",
         ylim=c(0,1), col="slateblue",
         xlab="Chr 1- position (Mbp)", ylab="Genotype probability",
         main="Line EB.09S.H.00002\n(after cleaning)",
         mgp.x=c(1.4,0.3,0), mgp.y=c(1.9,0.3,0), lwd=2)
lines(pmap[["10"]], prcl[["10"]]["EB.09S.H.00002","B73_inraB73_inra",], col="violetred", lwd=2)
legend("topright", lwd=2, col=c("slateblue", "violetred"),
       c("A632_usa", "B73_inra"), bg="gray92")
dev.off()
