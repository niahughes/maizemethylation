###################################################################
### GENOMATION DIAGRAMS
### R script for genomation diagrams 
###################################################################

# set up environment
setwd("/scratch/nia/manuFinal")
load(file=".RData")
library(methylKit)
library(GenomicRanges)
library(genomation)

#get gene annotation object
gene.obj <- readTranscriptFeatures("/home/nia/project/nia/annotation/geneAnnotationb73.bed", unique.prom=FALSE,  up.flank=2000, down.flank=0)
# list DM tiles by comparison
DMregions <- list(CpG.3x.100.oChi.50p.all, aCpG.3x.100.50p.all, bCpG.3x.100.50p.all, cCpG.3x.100.50p.all,
                  CHG.3x.100.oChi.50p.all, aCHG.3x.100.50p.all, bCHG.3x.100.50p.all, cCHG.3x.100.50p.all,
                  CHH.3x.100.oChi.50p.all, aCHH.3x.100.50p.all, bCHH.3x.100.50p.all, cCHH.3x.100.50p.all)
DMnames <- list("CpG", "aCpG", "bCpG", "cCpG", "CHG", "aCHG", "bCHG", "cCHG", "CHH", "aCHH", "bCHH", "cCHH")


Map(function(x,y) {
        print(y)
        gr <- as(x, "GRanges")
        
        pdf(file=paste0("/home/nia/project/nia/results/3.2-genomation/DMByFeature/" ,y, ".pdf")) # plot target annotation across features
        plotTargetAnnotation(annotateWithGeneParts(gr, gene.obj), precedence=TRUE, main="All DM regions by their annotated features", sub=y)
        dev.off()
        
        pdf(file=paste0("/home/nia/project/nia/results/3.2-genomation/DMByFeature/hyper-" ,y, ".pdf")) # hyper
        plotTargetAnnotation(annotateWithGeneParts(gr[gr$meth.diff>0], gene.obj), precedence=TRUE, main="Hypermethylated regions by their annotated features", sub=y)
        dev.off()
        
        pdf(file=paste0("/home/nia/project/nia/results/3.2-genomation/DMByFeature/hypo-" ,y, ".pdf")) # hypo
        plotTargetAnnotation(annotateWithGeneParts(gr[gr$meth.diff<0], gene.obj), precedence=TRUE, main="Hypomethylated regions by their annotated features", sub=y)
        dev.off()
        
        pdf(file=paste0("/home/nia/project/nia/results/3.2-genomation/DistFromTSS/",y,".pdf")) # plot distance from TSS
        tss.assoc=getAssociationWithTSS(annotateWithGeneParts(gr, gene.obj))
        hist(tss.assoc$dist.to.feature[abs(tss.assoc$dist.to.feature)<=100000],main=paste("Distance of DM", y,"tiles to nearest TSS"), xlab="Distance (bp)",
             breaks=50, col="darkred")
        dev.off()
        
        pdf(file=paste0("/home/nia/project/nia/results/3.2-genomation/DistFromTSS/hyper-",y,".pdf")) # hyper
        tss.assoc=getAssociationWithTSS(annotateWithGeneParts(gr[gr$meth.diff>0], gene.obj))
        hist(tss.assoc$dist.to.feature[abs(tss.assoc$dist.to.feature)<=100000],main=paste("Distance of hypermethylated", y,"tiles to nearest TSS"), xlab="Distance (bp)",
             breaks=50, col="darkred")
        dev.off()
        
        pdf(file=paste0("/home/nia/project/nia/results/3.2-genomation/DistFromTSS/hypo-",y,".pdf")) # hypo
        tss.assoc=getAssociationWithTSS(annotateWithGeneParts(gr[gr$meth.diff<0], gene.obj))
        hist(tss.assoc$dist.to.feature[abs(tss.assoc$dist.to.feature)<=100000],main=paste("Distance of hypomethylated", y,"tiles to nearest TSS"), xlab="Distance (bp)",
             breaks=50, col="darkred")
        dev.off()

        return(NULL)
}, DMregions, DMnames)

q(save="no")