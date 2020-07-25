###################################################################
### RELATIVE LINE PLOTS OF ABSOLUTE METHYLATION LEVELS PER SAMPLE
### ACROSS FEATURE TYPES
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


# Methylation percentages -------------------------------------------------
## Functions
readPercMeth <-function(infile,RDSname) {
  print(RDSname)
  x <- readGeneric(file=paste0("/home/nia/project/nia/results/bedgraphs/mKpercMeth/", infile), chr=1, start=2, end=3, meta.cols=list(perc.meth=4),skip=1)
  x <- GRanges(seqnames(x), IRanges(start(x), width=1), perc.meth=mcols(x)$perc.meth)
  print("subsetting")
  y <- subsetByOverlaps(x, chr2region)
  print("saving")
  saveRDS(y, file=paste0("3.5-data/", RDSname,".rds"))
  return(NULL)
}

backupReadPercMeth <- function(infile, RDSname) {
  print(RDSname)
  x=import(con=paste0("/home/nia/project/nia/results/bedgraphs/mKpercMeth/",infile), format="bedGraph")
  y=as(x, "GRanges")
  y$perc.meth<-y$score
  y$score <- NULL
  print("subsetting")
  z <- subsetByOverlaps(y, chr2region)
  print("saving")
  saveRDS(z, file=paste0("3.5-data/", RDSname, ".rds"))
  return(NULL)
}

## Read in (only needs to be done once)
chr2region <- GRanges(seqnames=c("1", "2", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
                      ranges=IRanges(start=c(1,1,182000000,1,1,1,1,1,1,1,1),
                                     end=c(307041717,22000000,244442276,235667834,246994605,223902240,174033170,182381542,181122637,159769782,150982314)),
                      strand="*")
seqlengths(chr2region)=c(307041717,244442276,235667834,246994605,223902240,174033170,182381542,181122637,159769782,150982314)

CpGIn <- list("sorted.aMutmKCpGpercMeth3x.bedgraph", "sorted.aWTmKCpGpercMeth3x.bedgraph", "sorted.bMutmKCpGpercMeth3x.bedgraph",
              "sorted.bWTmKCpGpercMeth3x.bedgraph", "sorted.cMutmKCpGpercMeth3x.bedgraph", "sorted.cWTmKCpGpercMeth3x.bedgraph")
CpGNames <- list("aMutCpG", "aWTCpG", "bMutCpG", "bWTCpG", "cMutCpG", "cWTCpG")
Map(readPercMeth, infile=CpGIn, RDSname=CpGNames)

CHGIn <- list("sorted.aMutmKCHGpercMeth3x.bedgraph", "sorted.aWTmKCHGpercMeth3x.bedgraph", "sorted.bMutmKCHGpercMeth3x.bedgraph",
              "sorted.bWTmKCHGpercMeth3x.bedgraph", "sorted.cMutmKCHGpercMeth3x.bedgraph", "sorted.cWTmKCHGpercMeth3x.bedgraph")
CHGNames <- list("aMutCHG", "aWTCHG", "bMutCHG", "bWTCHG", "cMutCHG", "cWTCHG")
Map(readPercMeth, infile=CHGIn, RDSname=CHGNames)

CHHIn <- list("sorted.aMutmKCHHpercMeth3x.bedgraph", "sorted.aWTmKCHHpercMeth3x.bedgraph",
              "sorted.bWTmKCHHpercMeth3x.bedgraph", "sorted.cMutmKCHHpercMeth3x.bedgraph", "sorted.cWTmKCHHpercMeth3x.bedgraph")
CHHNames <- list("aMutCHH", "aWTCHH", "bWTCHH", "cMutCHH", "cWTCHH")
Map(readPercMeth, infile=CHHIn, RDSname=CHHNames)

# To import bMut CHH, which does not import properly using readGeneric (no trailing numbers for first 68 lines
# of metacol so readGeneric gets upset when it reaches the first), and all pooled bedgraphs (also don't play nice
# with readGeneric)
Map(backupReadPercMeth, infile="sorted.bMutmKCHHpercMeth3x.bedgraph", RDSname="bMutCHH")
poolIn <- list("poolCpGMut.bedgraph", "poolCpGWT.bedgraph", "poolCHGMut.bedgraph", "poolCHGWT.bedgraph", "poolCHHMut.bedgraph", "poolCHHWT.bedgraph")
poolNames <- list("pool/CpGMut", "pool/CpGWT", "pool/CHGMut", "pool/CHGWT", "pool/CHHMut", "pool/CHHWT")
Map(backupReadPercMeth, infile=poolIn, RDSname=poolNames)

# Annotation objects ------------------------------------------------------
## Genes
geneObj=readBed("/home/nia/project/nia/annotation/geneAnnotationb73.bed", remove.unusual = TRUE)
geneBits <- readTranscriptFeatures("/home/nia/project/nia/annotation/geneAnnotationb73.bed", unique.prom=FALSE,  up.flank=2000, down.flank=0)
geneBits$TSSes <-NULL
## TEs - most recent annotation
te.file=keepStandardChromosomes(import.gff("/home/nia/project/nia/annotation/B73.structuralTEv2.fulllength.2018-09-19.gff3", format="gff3"), pruning.mode="coarse")
te.df=as.data.frame(te.file)
te.df$code=substr(te.df$ID,1,3)
te.df$class=ifelse(substr(te.df$ID,1,1) == "R",1,2)
teObj.all=makeGRangesFromDataFrame(te.df, keep.extra.columns = TRUE)
## MITEs
legacy.file=keepStandardChromosomes(import.gff("/home/nia/project/nia/annotation/legacy.gff3", format="gff3"), pruning.mode="coarse")
mite = legacy.file[legacy.file$source=="detectMITE"]

## Collect into one object
regionObjs <- list(genes=geneObj,
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
saveRDS(regionObjs, "3.5-data/regionObjs.rds")



# Pooled ------------------------------------------------------------------
## Functions
pooledRegionMatrix <- function(percMeth, featRegion, outName) {
  name <- print(sub('\\..*$', '', basename(percMeth)))
  name
  
  methyl <- readRDS(percMeth) # read in percent methylation bedgraphs
  
  up <- flank(featRegion, width=2000, start=TRUE, both=FALSE, use.names=TRUE, ignore.strand=FALSE) # get upstream flank of feature region
  down <- flank(featRegion, width=2000, start=FALSE, both=FALSE, use.names=TRUE, ignore.strand=FALSE) # get downstream flank of feature region
  
  up=ScoreMatrixBin(target=methyl, windows=up, bin.op="mean", bin.num=20, is.noCovNA = TRUE, weight.col="perc.meth")
  bin=ScoreMatrixBin(target=methyl, windows=featRegion, bin.op="mean", bin.num=60, is.noCovNA = TRUE, weight.col="perc.meth")
  down=ScoreMatrixBin(target=methyl, windows=down, bin.op="mean", bin.num=20, is.noCovNA = TRUE, weight.col="perc.meth")
  
  upPlot <- plotMeta(up)
  binPlot <- plotMeta(bin)
  downPlot <- plotMeta(down)
  
  mat <- setNames(as.data.frame(t(cbind(upPlot, binPlot, downPlot))), name)
  mat$bin <- 1:100
  mat <- melt(mat, id.vars=c("bin"), variable.name="genotype", value.name="perc.meth")
  mat$context <- as.factor(substr(mat$genotype, 1, 3))
  mat$status <- as.factor(substr(mat$genotype, 4, nchar(as.character(mat$genotype))))
  
  return(mat)
  
}

pooledPlotFeat <- function(matrices, names) {
  matrix <- readRDS(matrices)
  
  df <- do.call(rbind, matrix)
  rownames(df) <- NULL
  
  df$status <- factor(df$status, levels = c("WT", "Mut"))
  df$context <- factor(df$context, levels = c("CpG", "CHG", "CHH"))
  
  plot <- ggplot(df) +
    geom_line(aes(bin, perc.meth, linetype=status, colour=context)) +
    labs(y="% mC") +
    ylim(0,100) +
    theme_classic() +
    geom_vline(xintercept=20, color="gray") +
    geom_vline(xintercept=80, color="gray") +
    scale_color_lancet(name="Context") +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank())  +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.box = "horizontal"
    ) +
    scale_x_continuous(breaks = c(0,20, 80, 100), labels=c("-2kb", "Start", "End", "+2kb")) +
    scale_linetype_discrete(name="Methylation Status")
  
  title <- ggdraw() +
    draw_label(paste(names), fontface = 'bold', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  a1 <- plot_grid(title, plot, labels=NULL, ncol=1, rel_heights=c(.1,1))
  save_plot(paste0("/home/nia/project/nia/results/3.5-featureBodyPlots/POOL-", gsub(" ","",names), "FBP.png"), a1, base_asp = 1.6, ncol=1, nrow=1)
  return(a1)
  
}

## Create matrices of binned averages
files <- list.files(path="3.5-data/pool/", pattern="(CpG|CHG|CHH)", full.names=TRUE, recursive=FALSE)
names(files) <- sub('\\..*$', '', basename(list.files(path="3.5-data/pool/", pattern="(CpG|CHG|CHH)", full.names=TRUE, recursive=FALSE)))
regionObjs <- readRDS("3.5-data/regionObjs.rds")

genePoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$genes, outName="gene")
saveRDS(genePoolMat, file="3.5-data/pool/genePoolMat.rds")
allTEPoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$allTEs, outName="allTE")
saveRDS(allTEPoolMat, "3.5-data/pool/allTEPoolMat.rds")
rTEPoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$classI, outName="rTE")
saveRDS(rTEPoolMat, "3.5-data/pool/rTEPoolMat.rds")
dTEPoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$classII, outName="dTE")
saveRDS(dTEPoolMat, "3.5-data/pool/dTEPoolMat.rds")
ltrPoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$LTR, outName="ltr")
saveRDS(ltrPoolMat, "3.5-data/pool/ltrPoolMat.rds")
linePoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$LINE, outName="line")
saveRDS(linePoolMat, "3.5-data/pool/linePoolMat.rds")
sinePoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$SINE, outName="sine")
saveRDS(sinePoolMat, "3.5-data/pool/sinePoolMat.rds")
tirPoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$TIR, outName="tir")
saveRDS(tirPoolMat, "3.5-data/pool/tirPoolMat.rds")
rlcPoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$RLC, outName="rlc")
saveRDS(rlcPoolMat, "3.5-data/pool/rlcPoolMat.rds")
rlgPoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$RLG, outName="rlg")
saveRDS(rlgPoolMat, "3.5-data/pool/rlgPoolMat.rds")
rlxPoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$RLX, outName="rlx")
saveRDS(rlxPoolMat, "3.5-data/pool/rlxPoolMat.rds")
rilPoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$RIL, outName="ril")
saveRDS(rilPoolMat, "3.5-data/pool/rilPoolMat.rds")
ritPoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$RIT, outName="rit")
saveRDS(ritPoolMat, "3.5-data/pool/ritPoolMat.rds")
dtaPoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$DTA, outName="dta")
saveRDS(dtaPoolMat, "3.5-data/pool/dtaPoolMat.rds")
dtcPoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$DTC, outName="dtc")
saveRDS(dtcPoolMat, "3.5-data/pool/dtcPoolMat.rds")
dthPoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$DTH, outName="dth")
saveRDS(dthPoolMat, "3.5-data/pool/dthPoolMat.rds")
dtmPoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$DTM, outName="dtm")
saveRDS(dtmPoolMat, "3.5-data/pool/dtmPoolMat.rds")
dttPoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$DTT, outName="dtt")
saveRDS(dttPoolMat, "3.5-data/pool/dttPoolMat.rds")
dtxPoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$DTX, outName="dtx")
saveRDS(dtxPoolMat, "3.5-data/pool/dtxPoolMat.rds")
helPoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$helitron, outName="hel")
saveRDS(helPoolMat, "3.5-data/pool/helPoolMat.rds")
mitePoolMat <- lapply(files, pooledRegionMatrix, featRegion=regionObjs$MITE, outName="mite")
saveRDS(mitePoolMat, "3.5-data/pool/mitePoolMat.rds")

## Plot
matFiles <- list(genes="3.5-data/pool/genePoolMat.rds",
                 allTEs="3.5-data/pool/allTEPoolMat.rds",
                 classI="3.5-data/pool/rTEPoolMat.rds",
                 classII="3.5-data/pool/dTEPoolMat.rds",
                 LTR="3.5-data/pool/ltrPoolMat.rds",
                 LINE="3.5-data/pool/linePoolMat.rds",
                 SINE="3.5-data/pool/sinePoolMat.rds",
                 TIR="3.5-data/pool/tirPoolMat.rds",
                 helitron="3.5-data/pool/helPoolMat.rds",
                 RLC="3.5-data/pool/rlcPoolMat.rds",
                 RLG="3.5-data/pool/rlgPoolMat.rds",
                 RLX="3.5-data/pool/rlxPoolMat.rds",
                 RIL="3.5-data/pool/rilPoolMat.rds",
                 RIT="3.5-data/pool/ritPoolMat.rds",
                 DTA="3.5-data/pool/dtaPoolMat.rds",
                 DTC="3.5-data/pool/dtcPoolMat.rds",
                 DTH="3.5-data/pool/dthPoolMat.rds",
                 DTM="3.5-data/pool/dtmPoolMat.rds",
                 DTT="3.5-data/pool/dttPoolMat.rds",
                 DTX="3.5-data/pool/dtxPoolMat.rds",
                 MITE="3.5-data/pool/mitePoolMat.rds")
matNames <- list("Genes", "All transposable elements", "Retrotransposon", "DNA transposon", "LTR", "LINE", "SINE", "TIR",
                 "Helitron", "RLC", "RLG", "RLX", "RIL", "RIT",
                 "DTA", "DTC", "DTH", "DTM", "DTT", "DTX", "MITE")

outPlot <- Map(pooledPlotFeat, matFiles, matNames)

# Individual genotypes ----------------------------------------------------
indRegionMatrix <- function(percMeth, featRegion, outName) {
  name <- print(sub('\\..*$', '', basename(percMeth)))
  name
  
  methyl <- readRDS(percMeth) # read in percent methylation bedgraphs
  
  up <- flank(featRegion, width=2000, start=TRUE, both=FALSE, use.names=TRUE, ignore.strand=FALSE) # get upstream flank of feature region
  down <- flank(featRegion, width=2000, start=FALSE, both=FALSE, use.names=TRUE, ignore.strand=FALSE) # get downstream flank of feature region
  
  up=ScoreMatrixBin(target=methyl, windows=up, bin.op="mean", bin.num=20, is.noCovNA = TRUE, weight.col="perc.meth")
  bin=ScoreMatrixBin(target=methyl, windows=featRegion, bin.op="mean", bin.num=60, is.noCovNA = TRUE, weight.col="perc.meth")
  down=ScoreMatrixBin(target=methyl, windows=down, bin.op="mean", bin.num=20, is.noCovNA = TRUE, weight.col="perc.meth")
  
  upPlot <- plotMeta(up)
  binPlot <- plotMeta(bin)
  downPlot <- plotMeta(down)
  
  mat <- setNames(as.data.frame(t(cbind(upPlot, binPlot, downPlot))), name)
  mat$bin <- 1:100
  mat <- melt(mat, id.vars=c("bin"), variable.name="genotype", value.name="perc.meth")
  mat$background <- as.factor(toupper(substr(mat$genotype, 1, 1)))
  mat$status <- as.factor(substr(mat$genotype, 2, nchar(as.character(mat$genotype))-3))
  mat$context <- as.factor(substr(mat$genotype, nchar(as.character(mat$genotype))-2, nchar(as.character(mat$genotype))))
  
  return(mat)
  
}

indPlotFeat <- function(matrices, names) {
  matrix <- readRDS(matrices)
  
  df <- do.call(rbind, matrix)
  rownames(df) <- NULL
  
  df$status <- factor(df$status, levels = c("WT", "Mut"))
  df$context <- factor(df$context, levels = c("CpG", "CHG", "CHH"))
  
  plot <- ggplot(df) +
    geom_line(aes(bin, perc.meth, linetype=status, colour=context)) +
    labs(y="% mC") +
    ylim(0,100) +
    theme_classic() +
    geom_vline(xintercept=20, color="gray") +
    geom_vline(xintercept=80, color="gray") +
    scale_color_lancet(name="Context") +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank())  +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.box = "horizontal"
    ) +
    scale_x_continuous(breaks = c(0,20, 80, 100), labels=c("-2kb", "Start", "End", "+2kb")) +
    scale_linetype_discrete(name="Methylation Status")
  
  a <- plot %+% subset(df, background %in% c("A")) + ggtitle("A") + theme(legend.position="none")
  b <- plot %+% subset(df, background %in% c("B")) + ggtitle("B") + theme(legend.position="none")
  c <- plot %+% subset(df, background %in% c("C")) + ggtitle("C") + theme(legend.position="none")
  
  title <- ggdraw() +
    draw_label(paste0("Percent methylation across ", tolower(names), " bodies"), fontface = 'bold', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  l <- get_legend(plot)
  a1 <- plot_grid(e+annotate("text", x = 50, y = 99, label = paste0(names, " body")), g, h, ncol=1)
  a2 <- plot_grid(title, a1, l, labels=NULL, ncol=1, rel_heights = c(0.1, 1))
  save_plot(paste0("/home/nia/project/nia/results/3.5-featureBodyPlots/", gsub(" ","",names), "FBP.png"), a2, base_asp = 1.6, ncol=1, nrow=3)
  return(a2)
  
}

# Get matrices
files <- list.files(path="3.5-data/", pattern="(CpG|CHG|CHH)", full.names=TRUE, recursive=FALSE)
names(files) <- sub('\\..*$', '', basename(list.files(path="3.5-data/", pattern="(CpG|CHG|CHH)", full.names=TRUE, recursive=FALSE)))
regionObjs <- readRDS("3.5-data/regionObjs.rds")

geneIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$genes, outName="gene")
saveRDS(geneIndMat, file="3.5-data/ind/geneIndMat.rds")
allTEIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$allTEs, outName="allTE")
saveRDS(allTEIndMat, "3.5-data/ind/allTEIndMat.rds")
rTEIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$classI, outName="rTE")
saveRDS(rTEIndMat, "3.5-data/ind/rTEIndMat.rds")
dTEIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$classII, outName="dTE")
saveRDS(dTEIndMat, "3.5-data/ind/dTEIndMat.rds")
ltrIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$LTR, outName="ltr")
saveRDS(ltrIndMat, "3.5-data/ind/ltrIndMat.rds")
lineIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$LINE, outName="line")
saveRDS(lineIndMat, "3.5-data/ind/lineIndMat.rds")
sineIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$SINE, outName="sine")
saveRDS(sineIndMat, "3.5-data/ind/sineIndMat.rds")
tirIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$TIR, outName="tir")
saveRDS(tirIndMat, "3.5-data/ind/tirIndMat.rds")
rlcIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$RLC, outName="rlc")
saveRDS(rlcIndMat, "3.5-data/ind/rlcIndMat.rds")
rlgIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$RLG, outName="rlg")
saveRDS(rlgIndMat, "3.5-data/ind/rlgIndMat.rds")
rlxIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$RLX, outName="rlx")
saveRDS(rlxIndMat, "3.5-data/ind/rlxIndMat.rds")
rilIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$RIL, outName="ril")
saveRDS(rilIndMat, "3.5-data/ind/rilIndMat.rds")
ritIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$RIT, outName="rit")
saveRDS(ritIndMat, "3.5-data/ind/ritIndMat.rds")
dtaIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$DTA, outName="dta")
saveRDS(dtaIndMat, "3.5-data/ind/dtaIndMat.rds")
dtcIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$DTC, outName="dtc")
saveRDS(dtcIndMat, "3.5-data/ind/dtcIndMat.rds")
dthIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$DTH, outName="dth")
saveRDS(dthIndMat, "3.5-data/ind/dthIndMat.rds")
dtmIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$DTM, outName="dtm")
saveRDS(dtmIndMat, "3.5-data/ind/dtmIndMat.rds")
dttIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$DTT, outName="dtt")
saveRDS(dttIndMat, "3.5-data/ind/dttIndMat.rds")
dtxIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$DTX, outName="dtx")
saveRDS(dtxIndMat, "3.5-data/ind/dtxIndMat.rds")
helIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$helitron, outName="hel")
saveRDS(helIndMat, "3.5-data/ind/helIndMat.rds")
miteIndMat <- lapply(files, indRegionMatrix, featRegion=regionObjs$MITE, outName="mite")
saveRDS(miteIndMat, "3.5-data/ind/miteIndMat.rds")

#Plot
matFiles <- list(genes="3.5-data/ind/geneIndMat.rds",
                 allTEs="3.5-data/ind/allTEIndMat.rds",
                 classI="3.5-data/ind/rTEIndMat.rds",
                 classII="3.5-data/ind/dTEIndMat.rds",
                 LTR="3.5-data/ind/ltrIndMat.rds",
                 LINE="3.5-data/ind/lineIndMat.rds",
                 SINE="3.5-data/ind/sineIndMat.rds",
                 TIR="3.5-data/ind/tirIndMat.rds",
                 helitron="3.5-data/ind/helIndMat.rds",
                 RLC="3.5-data/ind/rlcIndMat.rds",
                 RLG="3.5-data/ind/rlgIndMat.rds",
                 RLX="3.5-data/ind/rlxIndMat.rds",
                 RIL="3.5-data/ind/rilIndMat.rds",
                 RIT="3.5-data/ind/ritIndMat.rds",
                 DTA="3.5-data/ind/dtaIndMat.rds",
                 DTC="3.5-data/ind/dtcIndMat.rds",
                 DTH="3.5-data/ind/dthIndMat.rds",
                 DTM="3.5-data/ind/dtmIndMat.rds",
                 DTT="3.5-data/ind/dttIndMat.rds",
                 DTX="3.5-data/ind/dtxIndMat.rds",
                 MITE="3.5-data/ind/miteIndMat.rds")
matNames <- list("Gene", "All transposable elements", "Retrotransposon", "DNA transposon", "LTR", "LINE", "SINE", "TIR",
                 "Helitron", "Copia", "Gypsy", "Unclassified LTR", "L1", "RTE",
                 "hAT", "CACTA", "PIF-Harbinger", "Mutator", "Tc1-Mariner", "Unclassified TIR", "MITE")
outPlot <- Map(indPlotFeat, matFiles, matNames)
outPlot <- indPlotFeat("3.5-data/ind/geneIndMat.rds", "Gene")
