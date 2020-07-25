###################################################################
### HEATMAPS OF RELATIVE METHYLATION CHANGE AT EACH FEATURE TYPE
### R script 
###################################################################

# set up environment
setwd("/scratch/nia/manuFinal")
library(methylKit)
library(GenomicRanges)
library(genomation)
library(rtracklayer)
library(GenomicFeatures)
library(reshape2)
library(cowplot)
library(ggplot2)
library(ggsci)
library(circlize)
library(ComplexHeatmap)
library(genomation)

# Annotation objects ------------------------------------------------------
# Only needs to be run once
# Genes
geneObj=readBed("/home/nia/project/nia/annotation/geneAnnotationb73.bed", remove.unusual = TRUE)
geneBits <- readTranscriptFeatures("/home/nia/project/nia/annotation/geneAnnotationb73.bed", unique.prom=FALSE,  up.flank=2000, down.flank=0)
# TEs - most recent annotation
te.file=keepStandardChromosomes(import.gff("/home/nia/project/nia/annotation/B73.structuralTEv2.fulllength.2018-09-19.gff3", format="gff3"), pruning.mode="coarse")
te.df=as.data.frame(te.file)
te.df$code=substr(te.df$ID,1,3)
te.df$class=ifelse(substr(te.df$ID,1,1) == "R",1,2)
teObj.all=makeGRangesFromDataFrame(te.df, keep.extra.columns = TRUE)
# MITEs
legacy.file=keepStandardChromosomes(import.gff("/home/nia/project/nia/annotation/legacy.gff3", format="gff3"), pruning.mode="coarse")
mite = legacy.file[legacy.file$source=="detectMITE"]

# Collect into one object
regionObjs <- list(genes=geneObj,
                   promoters=geneBits$promoters,
                   exons=geneBits$exons,
                   introns=geneBits$introns,
                   allTEs=teObj.all,
                   classI=teObj.all[teObj.all$class == 1],
                   classII=teObj.all[teObj.all$class == 2],
                   LTR=teObj.all[teObj.all$type == "LTR_retrotransposon" | teObj.all$type == "solo_LTR"],
                   LINE=teObj.all[teObj.all$type == "LINE_element"],
                   SINE=teObj.all[teObj.all$type == "SINE_element"],
                   TIR=teObj.all[teObj.all$type == "terminal_inverted_repeat_element"],
                   helitron=teObj.all[teObj.all$type == "helitron"],
                   RLC=teObj.all[teObj.all$code == "RLC"],
                   RLG=teObj.all[teObj.all$code == "RLG"],
                   RLX=teObj.all[teObj.all$code == "RLX"],
                   RIL=teObj.all[teObj.all$code == "RIL"],
                   RIT=teObj.all[teObj.all$code == "RIT"],
                   DTA=teObj.all[teObj.all$code == "DTA"],
                   DTC=teObj.all[teObj.all$code == "DTC"],
                   DTH=teObj.all[teObj.all$code == "DTH"],
                   DTM=teObj.all[teObj.all$code == "DTM"],
                   DTT=teObj.all[teObj.all$code == "DTT"],
                   DTX=teObj.all[teObj.all$code == "DTX"],
                   MITE=mite)
saveRDS(regionObjs, "3.8-data/regionObjs.rds")


# Heatmaps ----------------------------------------------------------------
# Functions
averageAcrossRegion <- function(pMeth, context) { #get average methylation across each region
  print(sub('\\..*$', '', basename(pMeth)))
  
  methyl <- readRDS(pMeth) # read in percent methylation bedgraph
  
  out <- lapply(regionObjs, function(x) { # calculate for every annotation type
    print(x$code[1])
    y <- subsetByOverlaps(methyl, x)
    return(mean(y$perc.meth))
  })
  
  return(out)
  
}

averageAcrossAverages <- function(pMeth, context) { # get average methylation across each region and then average the regions
  print(sub('\\..*$', '', basename(pMeth)))
  
  methyl <- readRDS(pMeth) # read in percent methylation bedgraph
  
  out <- lapply(regionObjs, function(x) { # calculate for every annotation type
    print(x$code[1])
    bin=ScoreMatrixBin(target=methyl, windows=x, bin.op="mean", bin.num=1, is.noCovNA = TRUE, weight.col="perc.meth") #score
    avg <- plotMeta(bin) # average
    output <- list(matrix=bin@.Data, avg=avg)
    return(output)
  })
  
  return(out)
  
}



# Mean across genome ------------------------------------------------------
files <- list.files(path="3.5-data/", pattern="(CpG|CHG|CHH)", full.names=TRUE, recursive=FALSE)
names(files) <- sub('\\..*$', '', basename(list.files(path="3.5-data/", pattern="(CpG|CHG|CHH)", full.names=TRUE, recursive=FALSE)))

percMeth <- lapply(files, function(x) {
  y <- readRDS(x)
  z <- mean(y$perc.meth)
  return(z)
})


# Average across regions --------------------------------------------------
files <- list.files(path="3.5-data/", pattern="(CpG|CHG|CHH)", full.names=TRUE, recursive=FALSE)
names(files) <- sub('\\..*$', '', basename(list.files(path="3.5-data/", pattern="(CpG|CHG|CHH)", full.names=TRUE, recursive=FALSE)))
regionObjs <- readRDS("3.8-data/regionObjs.rds")

CpG.averages <- lapply(files[c("aMutCpG", "aWTCpG", "bMutCpG", "bWTCpG", "cMutCpG", "cWTCpG")], averageAcrossRegion, context="CpG")
saveRDS(CpG.averages, file="3.8-data/avg/CpGaverages.rds")
CHG.averages <- lapply(files[c("aMutCHG", "aWTCHG", "bMutCHG", "bWTCHG", "cMutCHG", "cWTCHG")], averageAcrossRegion, context="CHG")
saveRDS(CHG.averages, file="3.8-data/avg/CHGaverages.rds")
CHH.averages <- lapply(files[c("aMutCHH", "aWTCHH", "bMutCHH", "bWTCHH", "cMutCHH", "cWTCHH")], averageAcrossRegion, context="CHH")
saveRDS(CHH.averages, file="3.8-data/avg/CHHaverages.rds")

## Calculate relative methylation loss
## Average of methylation calls mapping to feature
contextList <- list(CpG = "3.8-data/avg/CpGaverages.rds", CHG="3.8-data/avg/CHGaverages.rds", CHH="3.8-data/avg/CHHaverages.rds")
heatMat <- lapply(contextList, function(a) {
  x <- readRDS(a)
  names(x) <- c("aMut", "aWT", "bMut", "bWT", "cMut", "cWT")
  A <- c((x$aMut$genes-x$aWT$genes)/x$aWT$genes, (x$aMut$promoters-x$aWT$promoters)/x$aWT$promoters, (x$aMut$exons-x$aWT$exons)/x$aWT$exons,
         (x$aMut$introns-x$aWT$introns)/x$aWT$introns, (x$aMut$allTEs-x$aWT$allTEs)/x$aWT$allTEs, (x$aMut$classI-x$aWT$classI)/x$aWT$classI,
         (x$aMut$LTR-x$aWT$LTR)/x$aWT$LTR, (x$aMut$RLC-x$aWT$RLC)/x$aWT$RLC, (x$aMut$RLG-x$aWT$RLG)/x$aWT$RLG, (x$aMut$RLX-x$aWT$RLX)/x$aWT$RLX,
         (x$aMut$LINE-x$aWT$LINE)/x$aWT$LINE, (x$aMut$RIL-x$aWT$RIL)/x$aWT$RIL, (x$aMut$RIT-x$aWT$RIT)/x$aWT$RIT,
         (x$aMut$SINE-x$aWT$SINE)/x$aWT$SINE, (x$aMut$classII-x$aWT$classII)/x$aWT$classII, (x$aMut$TIR-x$aWT$TIR)/x$aWT$TIR,  (x$aMut$DTA-x$aWT$DTA)/x$aWT$DTA,
         (x$aMut$DTC-x$aWT$DTC)/x$aWT$DTC, (x$aMut$DTH-x$aWT$DTH)/x$aWT$DTH, (x$aMut$DTM-x$aWT$DTM)/x$aWT$DTM,
         (x$aMut$DTT-x$aWT$DTT)/x$aWT$DTT, (x$aMut$DTX-x$aWT$DTX)/x$aWT$DTX, (x$aMut$MITE-x$aWT$MITE)/x$aWT$MITE, (x$aMut$helitron-x$aWT$helitron)/x$aWT$helitron)
  B <- c((x$bMut$genes-x$bWT$genes)/x$bWT$genes, (x$bMut$promoters-x$bWT$promoters)/x$bWT$promoters, (x$bMut$exons-x$bWT$exons)/x$bWT$exons,
         (x$bMut$introns-x$bWT$introns)/x$bWT$introns, (x$bMut$allTEs-x$bWT$allTEs)/x$bWT$allTEs, (x$bMut$classI-x$bWT$classI)/x$bWT$classI,
         (x$bMut$LTR-x$bWT$LTR)/x$bWT$LTR, (x$bMut$RLC-x$bWT$RLC)/x$bWT$RLC, (x$bMut$RLG-x$bWT$RLG)/x$bWT$RLG, (x$bMut$RLX-x$bWT$RLX)/x$bWT$RLX,
         (x$bMut$LINE-x$bWT$LINE)/x$bWT$LINE, (x$bMut$RIL-x$bWT$RIL)/x$bWT$RIL, (x$bMut$RIT-x$bWT$RIT)/x$bWT$RIT,
         (x$bMut$SINE-x$bWT$SINE)/x$bWT$SINE, (x$bMut$classII-x$bWT$classII)/x$bWT$classII, (x$bMut$TIR-x$bWT$TIR)/x$bWT$TIR,  (x$bMut$DTA-x$bWT$DTA)/x$bWT$DTA,
         (x$bMut$DTC-x$bWT$DTC)/x$bWT$DTC, (x$bMut$DTH-x$bWT$DTH)/x$bWT$DTH, (x$bMut$DTM-x$bWT$DTM)/x$bWT$DTM,
         (x$bMut$DTT-x$bWT$DTT)/x$bWT$DTT, (x$bMut$DTX-x$bWT$DTX)/x$bWT$DTX, (x$bMut$MITE-x$bWT$MITE)/x$bWT$MITE, (x$bMut$helitron-x$bWT$helitron)/x$bWT$helitron)
  C <- c((x$cMut$genes-x$cWT$genes)/x$cWT$genes, (x$cMut$promoters-x$cWT$promoters)/x$cWT$promoters, (x$cMut$exons-x$cWT$exons)/x$cWT$exons,
         (x$cMut$introns-x$cWT$introns)/x$cWT$introns, (x$cMut$allTEs-x$cWT$allTEs)/x$cWT$allTEs, (x$cMut$classI-x$cWT$classI)/x$cWT$classI,
         (x$cMut$LTR-x$cWT$LTR)/x$cWT$LTR, (x$cMut$RLC-x$cWT$RLC)/x$cWT$RLC, (x$cMut$RLG-x$cWT$RLG)/x$cWT$RLG, (x$cMut$RLX-x$cWT$RLX)/x$cWT$RLX,
         (x$cMut$LINE-x$cWT$LINE)/x$cWT$LINE, (x$cMut$RIL-x$cWT$RIL)/x$cWT$RIL, (x$cMut$RIT-x$cWT$RIT)/x$cWT$RIT,
         (x$cMut$SINE-x$cWT$SINE)/x$cWT$SINE, (x$cMut$classII-x$cWT$classII)/x$cWT$classII, (x$cMut$TIR-x$cWT$TIR)/x$cWT$TIR,  (x$cMut$DTA-x$cWT$DTA)/x$cWT$DTA,
         (x$cMut$DTC-x$cWT$DTC)/x$cWT$DTC, (x$cMut$DTH-x$cWT$DTH)/x$cWT$DTH, (x$cMut$DTM-x$cWT$DTM)/x$cWT$DTM,
         (x$cMut$DTT-x$cWT$DTT)/x$cWT$DTT, (x$cMut$DTX-x$cWT$DTX)/x$cWT$DTX, (x$cMut$MITE-x$cWT$MITE)/x$cWT$MITE, (x$cMut$helitron-x$cWT$helitron)/x$cWT$helitron)
  mat <- cbind(E, G, H)
  rownames(mat) <- c("Genes", "  Promoters", "  Exons", "  Introns", "All TEs", "  Class I", "    LTR", "      RLC", "      RLG", "      RLX",
                     "    LINE", "      RIL", "      RIT", "    SINE", "  Class II","    TIR",
                     "      DTA", "      DTC", "      DTH", "      DTM", "      DTT", "      DTX", "      MITE", "   Helitron")
  mat <- mat*100
  return(mat)
})

# Plot heatmap
col_fun = colorRamp2(c(-20, 0, 20), c("red", "white", "blue"))
pdf("/home/nia/project/nia/results/3.8-heatmaps/heatmapSymm.pdf")
heat <- Heatmap(heatMat$CpG, col=col_fun, cluster_rows = FALSE, cluster_columns = FALSE, column_title = "CpG", show_heatmap_legend = TRUE,
                heatmap_legend_param = list(title = "Relative methylation\ndifference", direction="horizontal")) +
  Heatmap(heatMat$CHG, col=col_fun, cluster_rows = FALSE, cluster_columns = FALSE, column_title = "CHG", show_heatmap_legend = FALSE)
draw(heat,  heatmap_legend_side = "bottom")
dev.off()

col_fun2 = colorRamp2(c(-100, 0, 100), c("red", "white", "blue"))
pdf("/home/nia/project/nia/results/3.8-heatmaps/heatmapCHH.pdf")
heat2 <- Heatmap(heatMat$CHH, col=col_fun2, cluster_rows = FALSE, cluster_columns = FALSE,column_title = "CHH", show_heatmap_legend = TRUE,
                 heatmap_legend_param = list(title = "Relative methylation\ndifference", direction="horizontal"))
draw(heat2,  heatmap_legend_side = "bottom")
dev.off()




# Meta-average ------------------------------------------------------------

regionObjs$exons[width(regionObjs$exons) == 1] <- NULL # Eliminate exons of length 1 to prevent scoreMatrixBin errors

CpG.averages.plotMeta <- lapply(files[c("aMutCpG", "aWTCpG", "bMutCpG", "bWTCpG", "cMutCpG", "cWTCpG")], averageAcrossAverages, context="CpG")
saveRDS(CpG.averages.plotMeta, file="3.8-data/avg/CpGaverages-plotMeta.rds")
CHG.averages.plotMeta <- lapply(files[c("aMutCHG", "aWTCHG", "bMutCHG", "bWTCHG", "cMutCHG", "cWTCHG")], averageAcrossAverages, context="CHG")
saveRDS(CHG.averages.plotMeta, file="3.8-data/avg/CHGaverages-plotMeta.rds")
CHH.averages.plotMeta <- lapply(files[c("aMutCHH", "aWTCHH", "bMutCHH", "bWTCHH", "cMutCHH", "cWTCHH")], averageAcrossAverages, context="CHH")
saveRDS(CHH.averages.plotMeta, file="3.8-data/avg/CHHaverages-plotMeta.rds")

## Weighted average of per-feature means
contextListMeta <- list(CpG = "3.8-data/avg/CpGaverages-plotMeta.rds", CHG="3.8-data/avg/CHGaverages-plotMeta.rds", CHH="3.8-data/avg/CHHaverages-plotMeta.rds")
heatMeta <- lapply(contextListMeta, function(a) {
  x <- readRDS(a)
  names(x) <- c("aMut", "aWT", "bMut", "bWT", "cMut", "cWT")
  A <- c((x$aMut$genes$avg-x$aWT$genes$avg)/x$aWT$genes$avg, (x$aMut$promoters$avg-x$aWT$promoters$avg)/x$aWT$promoters$avg, 
         (x$aMut$exons$avg-x$aWT$exons$avg)/x$aWT$exons$avg, (x$aMut$introns$avg-x$aWT$introns$avg)/x$aWT$introns$avg,
         (x$aMut$allTEs$avg-x$aWT$allTEs$avg)/x$aWT$allTEs$avg, (x$aMut$classI$avg-x$aWT$classI$avg)/x$aWT$classI$avg,
         (x$aMut$LTR$avg-x$aWT$LTR$avg)/x$aWT$LTR$avg, (x$aMut$RLC$avg-x$aWT$RLC$avg)/x$aWT$RLC$avg, 
         (x$aMut$RLG$avg-x$aWT$RLG$avg)/x$aWT$RLG$avg, (x$aMut$RLX$avg-x$aWT$RLX$avg)/x$aWT$RLX$avg,
         (x$aMut$LINE$avg-x$aWT$LINE$avg)/x$aWT$LINE$avg, (x$aMut$RIL$avg-x$aWT$RIL$avg)/x$aWT$RIL$avg, 
         (x$aMut$RIT$avg-x$aWT$RIT$avg)/x$aWT$RIT$avg, (x$aMut$SINE$avg-x$aWT$SINE$avg)/x$aWT$SINE$avg, 
         (x$aMut$classII$avg-x$aWT$classII$avg)/x$aWT$classII$avg, (x$aMut$TIR$avg-x$aWT$TIR$avg)/x$aWT$TIR$avg, 
         (x$aMut$DTA$avg-x$aWT$DTA$avg)/x$aWT$DTA$avg, (x$aMut$DTC$avg-x$aWT$DTC$avg)/x$aWT$DTC$avg, 
         (x$aMut$DTH$avg-x$aWT$DTH$avg)/x$aWT$DTH$avg, (x$aMut$DTM$avg-x$aWT$DTM$avg)/x$aWT$DTM$avg,
         (x$aMut$DTT$avg-x$aWT$DTT$avg)/x$aWT$DTT$avg, (x$aMut$DTX$avg-x$aWT$DTX$avg)/x$aWT$DTX$avg,
         (x$aMut$MITE$avg-x$aWT$MITE$avg)/x$aWT$MITE$avg, (x$aMut$helitron$avg-x$aWT$helitron$avg)/x$aWT$helitron$avg)
  B <- c((x$bMut$genes$avg-x$bWT$genes$avg)/x$bWT$genes$avg, (x$bMut$promoters$avg-x$bWT$promoters$avg)/x$bWT$promoters$avg, 
         (x$bMut$exons$avg-x$bWT$exons$avg)/x$bWT$exons$avg, (x$bMut$introns$avg-x$bWT$introns$avg)/x$bWT$introns$avg,
         (x$bMut$allTEs$avg-x$bWT$allTEs$avg)/x$bWT$allTEs$avg, (x$bMut$classI$avg-x$bWT$classI$avg)/x$bWT$classI$avg,
         (x$bMut$LTR$avg-x$bWT$LTR$avg)/x$bWT$LTR$avg, (x$bMut$RLC$avg-x$bWT$RLC$avg)/x$bWT$RLC$avg, 
         (x$bMut$RLG$avg-x$bWT$RLG$avg)/x$bWT$RLG$avg, (x$bMut$RLX$avg-x$bWT$RLX$avg)/x$bWT$RLX$avg,
         (x$bMut$LINE$avg-x$bWT$LINE$avg)/x$bWT$LINE$avg, (x$bMut$RIL$avg-x$bWT$RIL$avg)/x$bWT$RIL$avg, 
         (x$bMut$RIT$avg-x$bWT$RIT$avg)/x$bWT$RIT$avg, (x$bMut$SINE$avg-x$bWT$SINE$avg)/x$bWT$SINE$avg, 
         (x$bMut$classII$avg-x$bWT$classII$avg)/x$bWT$classII$avg, (x$bMut$TIR$avg-x$bWT$TIR$avg)/x$bWT$TIR$avg, 
         (x$bMut$DTA$avg-x$bWT$DTA$avg)/x$bWT$DTA$avg, (x$bMut$DTC$avg-x$bWT$DTC$avg)/x$bWT$DTC$avg, 
         (x$bMut$DTH$avg-x$bWT$DTH$avg)/x$bWT$DTH$avg, (x$bMut$DTM$avg-x$bWT$DTM$avg)/x$bWT$DTM$avg,
         (x$bMut$DTT$avg-x$bWT$DTT$avg)/x$bWT$DTT$avg, (x$bMut$DTX$avg-x$bWT$DTX$avg)/x$bWT$DTX$avg,
         (x$bMut$MITE$avg-x$bWT$MITE$avg)/x$bWT$MITE$avg, (x$bMut$helitron$avg-x$bWT$helitron$avg)/x$bWT$helitron$avg)
  C <- c((x$cMut$genes$avg-x$cWT$genes$avg)/x$cWT$genes$avg, (x$cMut$promoters$avg-x$cWT$promoters$avg)/x$cWT$promoters$avg, 
         (x$cMut$exons$avg-x$cWT$exons$avg)/x$cWT$exons$avg, (x$cMut$introns$avg-x$cWT$introns$avg)/x$cWT$introns$avg,
         (x$cMut$allTEs$avg-x$cWT$allTEs$avg)/x$cWT$allTEs$avg, (x$cMut$classI$avg-x$cWT$classI$avg)/x$cWT$classI$avg,
         (x$cMut$LTR$avg-x$cWT$LTR$avg)/x$cWT$LTR$avg, (x$cMut$RLC$avg-x$cWT$RLC$avg)/x$cWT$RLC$avg, 
         (x$cMut$RLG$avg-x$cWT$RLG$avg)/x$cWT$RLG$avg, (x$cMut$RLX$avg-x$cWT$RLX$avg)/x$cWT$RLX$avg,
         (x$cMut$LINE$avg-x$cWT$LINE$avg)/x$cWT$LINE$avg, (x$cMut$RIL$avg-x$cWT$RIL$avg)/x$cWT$RIL$avg, 
         (x$cMut$RIT$avg-x$cWT$RIT$avg)/x$cWT$RIT$avg, (x$cMut$SINE$avg-x$cWT$SINE$avg)/x$cWT$SINE$avg, 
         (x$cMut$classII$avg-x$cWT$classII$avg)/x$cWT$classII$avg, (x$cMut$TIR$avg-x$cWT$TIR$avg)/x$cWT$TIR$avg, 
         (x$cMut$DTA$avg-x$cWT$DTA$avg)/x$cWT$DTA$avg, (x$cMut$DTC$avg-x$cWT$DTC$avg)/x$cWT$DTC$avg, 
         (x$cMut$DTH$avg-x$cWT$DTH$avg)/x$cWT$DTH$avg, (x$cMut$DTM$avg-x$cWT$DTM$avg)/x$cWT$DTM$avg,
         (x$cMut$DTT$avg-x$cWT$DTT$avg)/x$cWT$DTT$avg, (x$cMut$DTX$avg-x$cWT$DTX$avg)/x$cWT$DTX$avg,
         (x$cMut$MITE$avg-x$cWT$MITE$avg)/x$cWT$MITE$avg, (x$cMut$helitron$avg-x$cWT$helitron$avg)/x$cWT$helitron$avg)
  mat <- cbind(E, G, H)
  rownames(mat) <- c("Genes", "  Promoters", "  Exons", "  Introns", "All TEs", "  Class I", "    LTR", "      RLC", "      RLG", "      RLX",
                     "    LINE", "      RIL", "      RIT", "    SINE", "  Class II","    TIR",
                     "      DTA", "      DTC", "      DTH", "      DTM", "      DTT", "      DTX", "      MITE", "   Helitron")
  mat <- mat*100
  return(mat)
})

# # Plot heatmap
col_fun = colorRamp2(c(-20, 0, 20), c("red", "white", "blue"))
pdf("/home/nia/project/nia/results/3.8-heatmaps/heatmapSymm-meta.pdf")
Heatmap(heatMeta$CpG, col=col_fun, cluster_rows = FALSE, cluster_columns = FALSE, column_title = "CpG", show_heatmap_legend = TRUE, heatmap_legend_param = list(title = "Relative\nmethylation\ndifference")) +
  Heatmap(heatMeta$CHG, col=col_fun, cluster_rows = FALSE, cluster_columns = FALSE, column_title = "CHG", show_heatmap_legend = FALSE)
dev.off()

col_fun2 = colorRamp2(c(-100, 0, 100), c("red", "white", "blue"))
pdf("/home/nia/project/nia/results/3.8-heatmaps/heatmapCHH-meta.pdf")
Heatmap(heatMeta$CHH, col=col_fun2, cluster_rows = FALSE, cluster_columns = FALSE,column_title = "CHH", show_heatmap_legend = TRUE, heatmap_legend_param = list(title = "Relative\nmethylation\ndifference"))
dev.off()



# Histograms --------------------------------------------------------------

cL2 <- c(readRDS("3.8-data/avg/CpGaverages-plotMeta.rds"), readRDS("3.8-data/avg/CHGaverages-plotMeta.rds"), readRDS("3.8-data/avg/CHHaverages-plotMeta.rds")) # convert to one list
cL3 <- lapply(cL2, function(x) { # remove average 
  a <- lapply(x, function(y) do.call(rbind, y))
  return(a)
})
cL <- melt(cL3)
colnames(cL) <- c("rowname", "filler", "value", "feature", "background")
cL$status <- substr(cL$background, 2, nchar(cL$background)-3)
cL$geno <- toupper(substr(cL$background, 1, 1))
cL$context <- substr(cL$background, nchar(cL$background)-2, nchar(cL$background))
cL <- na.omit(cL)
cL.a <- split(cL, f=cL$feature)


aCpGgenes <- ggplot(subset(cL.a$genes, geno == "E" & context=="CpG"), aes(x=value, y=stat(density*width), fill = status)) + 
  geom_histogram(alpha=0.2, bins=30, position="identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title="A CpG genes", x="Percent methylation", y="% features", color = "Legend") +
  scale_fill_manual(values=c("red", "blue")) +
  theme(legend.position = "none")
bCpGgenes <- ggplot(subset(cL.a$genes, geno == "G" & context=="CpG"), aes(x=value, y=stat(density*width), fill = status)) + 
  geom_histogram(alpha=0.2, bins=30, position="identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title="B CpG genes", x="Percent methylation", y="% features", color = "Legend") +
  scale_fill_manual(values=c("red", "blue")) +
  theme(legend.position = "none")
cCpGgenes <- ggplot(subset(cL.a$genes, geno == "H" & context=="CpG"), aes(x=value, y=stat(density*width), fill = status)) + 
  geom_histogram(alpha=0.2, bins=30, position="identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title="C CpG genes", x="Percent methylation", y="% features", color = "Legend") +
  scale_fill_manual(values=c("red", "blue")) +
  theme(legend.position = "none")

aCHGintrons <- ggplot(subset(cL.a$introns, geno == "E" & context=="CHG"), aes(x=value, y=stat(density*width), fill = status)) + 
  geom_histogram(alpha=0.2, bins=30, position="identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title="A CHG introns",x="Percent methylation", y="% features", color = "Legend") +
  scale_fill_manual(values=c("red", "blue")) +
  theme(legend.position = "none")
bCHGintrons <- ggplot(subset(cL.a$introns, geno == "G" & context=="CHG"), aes(x=value, y=stat(density*width), fill = status)) + 
  geom_histogram(alpha=0.2, bins=30, position="identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title="B CHG introns",x="Percent methylation", y="% features", color = "Legend") +
  scale_fill_manual(values=c("red", "blue")) +
  theme(legend.position = "none") 
cCHGintrons <- ggplot(subset(cL.a$introns, geno == "H" & context=="CHG"), aes(x=value, y=stat(density*width), fill = status)) + 
  geom_histogram(alpha=0.2, bins=30, position="identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title="C CHG introns", x="Percent methylation", y="% features", color = "Legend") +
  scale_fill_manual(values=c("red", "blue")) +
  theme(legend.position = "none")

aCHHMITE <- ggplot(subset(cL.a$MITE, geno == "E" & context=="CHH"), aes(x=value, y=stat(density*width), fill = status)) + 
  geom_histogram(alpha=0.2, bins=30, position="identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title="A CHH MITEs", x="Percent methylation", y="% features", color = "Legend") +
  scale_fill_manual(values=c("red", "blue")) +
  theme(legend.position = "none")
bCHHMITE <- ggplot(subset(cL.a$MITE, geno == "G" & context=="CHH"), aes(x=value, y=stat(density*width), fill = status)) + 
  geom_histogram(alpha=0.2, bins=30, position="identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title="B CHH MITEs",x="Percent methylation", y="% features", color = "Legend") +
  scale_fill_manual(values=c("red", "blue")) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center") 
cCHHMITE <- ggplot(subset(cL.a$MITE, geno == "H" & context=="CHH"), aes(x=value, y=stat(density*width), fill = status)) + 
  geom_histogram(alpha=0.2, bins=30, position="identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title="C CHH MITEs", x="Percent methylation", y="% features", color = "Legend") +
  scale_fill_manual(values=c("red", "blue")) +
  theme(legend.position = "none")

theme_set(theme_cowplot(font_size=8))

leg <- get_legend(bCHHMITE)
a1 <- plot_grid(aCpGgenes + coord_cartesian(ylim=c(0,0.12)),  bCpGgenes + coord_cartesian(ylim=c(0,0.12)),  cCpGgenes + coord_cartesian(ylim=c(0,0.12)), 
          aCHGintrons + coord_cartesian(ylim=c(0,0.03)), bCHGintrons + coord_cartesian(ylim=c(0,0.03)), cCHGintrons + coord_cartesian(ylim=c(0,0.03)), 
          aCHHMITE + coord_cartesian(ylim=c(0,0.45)), bCHHMITE + coord_cartesian(ylim=c(0,0.45)) + theme(legend.position="none"), cCHHMITE + coord_cartesian(ylim=c(0,0.45)), 
          labels="AUTO", ncol=3)
a2 <- plot_grid(a1, leg, rel_heights = c(1, 0.03), ncol=1)

a3 <- plot_grid(aCHGintrons + coord_cartesian(ylim=c(0,0.7)), bCHGintrons + coord_cartesian(ylim=c(0,0.7)), cCHGintrons + coord_cartesian(ylim=c(0,0.7)),
                labels="AUTO", ncol=3)

a4 <- plot_grid(a3, leg, rel_heights=c(1 ,0.07), ncol=1)