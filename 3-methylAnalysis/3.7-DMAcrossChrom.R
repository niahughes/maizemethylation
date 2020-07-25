###################################################################
### RELATIVE LINE PLOTS OF ABSOLUTE AND DIFFERENTIAL METHYLATION 
### ACROSS GENOME AND INDIVIDUAL CHROMOSOMES
### R script 
###################################################################

# set up environment
setwd("/scratch/nia/manuFinal")
load(file=".RData")
library(GenomicRanges)
library(cowplot)
library(ggplot2)
library(ggsci)


# ***Differential methylation ---------------------------------------------
# Functions ---------------------------------------------------------------

methDiffMatrix <- function (obj, name, binsize){
  print(name)
  
  #preprocess data
  methyl <- as(obj, "GRanges")
  chrLengths <- c(307041717,150982314,244442276,235667834,246994605,223902240,174033170,182381542,181122637,159769782)
  seqlengths(methyl) <- chrLengths # chromosome lengths
  seqlevels(methyl) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
  
  #bin genome into chunks
  bins <- tileGenome(seqinfo(methyl), tilewidth=binsize,cut.last.tile.in.chrom=TRUE)
  GRbins <- as(bins, "GRanges")
  
  #add NA to regions without methylation counts
  allGaps <- gaps(methyl)
  impGaps <- allGaps[strand(allGaps) == "*"] # remove start-to-finish ranges (appear as + and - strands)
  impGaps$perc.meth <- NA
  
  #ensure that there is no overlap between NAs and percent methylation
  if (length(findOverlaps(impGaps, methyl)) != 0) {
    print("OVERLAPS APPEAR, FIX ME")
  }
  
  #merge and sort ranges
  fullRange <- c(impGaps, methyl)
  fullRangeSort <- sort(fullRange)
  
  #calculate Rle and bin the averages for each window
  score <- GenomicRanges::coverage(fullRangeSort, weight="meth.diff")
  chromMat <- GenomicRanges::binnedAverage(GRbins, score, "meth.diff", na.rm=TRUE)
  chromDF <- as.data.frame(chromMat)
  chromDF$bin <- 1:nrow(chromDF)
  chromDF$genotype <- name
  chromDF$background <- as.factor(toupper(substr(chromDF$genotype, 1, 1)))
  chromDF$context <- as.factor(substr(chromDF$genotype, 2, nchar(as.character(chromDF$genotype))))
  outDF = subset(chromDF, select = -c(start, end, width, strand))
  
  return(outDF)
}

plotDMAcrossChrom <- function(chrom, percValues, plotOutput=FALSE, binsize, rect=FALSE) {
  newMat2 <- percValues[percValues$seqnames==chrom,]
  newMat2$mb<-newMat2$bin-min(newMat2$bin)+1
  
  mop1pos <- 42155840/binsize
  
  plotChrom2<-ggplot(newMat2) +
    geom_line(aes(mb, meth.diff, colour=context)) +
    labs(y="% DM", x=paste("Position on chromosome", chrom, "in Mb")) +
    ylim(-10,10) +
    theme_classic() +
    scale_color_lancet(name="context") +
    theme(legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal") + 
    scale_x_continuous(expand = c(0, 0)) +
    {if(chrom==2) geom_vline(xintercept = mop1pos, colour="yellow")} +
    {if (chrom==2 && rect==TRUE) annotate("rect", xmin=330-308, xmax=490-308, ymin=-Inf, ymax=Inf, fill="blue", alpha=0.2)} 
  
  aChrom2 <- plotChrom2 %+% subset(newMat2, background %in% c("A")) + ggtitle("A") + theme(legend.position="none") + theme(axis.title.x=element_blank())
  bChrom2 <- plotChrom2 %+% subset(newMat2, background %in% c("B")) + ggtitle("B") + theme(legend.position="none") + theme(axis.title.x=element_blank())
  cChrom2 <- plotChrom2 %+% subset(newMat2, background %in% c("C")) + ggtitle("C") + theme(legend.position="none")
  
  title.Chrom2 <- ggdraw() + 
    draw_label(paste0("Percent DM across chromosome ", chrom), fontface = 'bold', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  l.Chrom2 <- get_legend(plotChrom2)
  a1.Chrom2 <- plot_grid(aChrom2, bChrom2, cChrom2, labels="AUTO", ncol=1)
  a2.Chrom2 <- plot_grid(title.Chrom2, a1.Chrom2, l.Chrom2, labels=NULL, ncol=1, rel_heights = c(0.1, 1))
  if (plotOutput==T && chrom==2 && rect==T) {
    save_plot(paste0("/home/nia/project/nia/results/3.7-DMAcrossChrom/DM-rect-chrom", chrom, ".png"), a2.Chrom2, base_asp = 2, ncol=1, nrow=3)
  } else if (plotOutput==T) {
    save_plot(paste0("/home/nia/project/nia/results/3.7-DMAcrossChrom/DM-chrom", chrom, ".png"), a2.Chrom2, base_asp = 2, ncol=1, nrow=3)
  }
  return(a2.Chrom2)
}

# Get matrices ------------------------------------------------------------
# Read in list
percList <- list(diff.aCpG.3x.100, diff.aCHG.3x.100, diff.aCHH.3x.100,
                 diff.bCpG.3x.100, diff.bCHG.3x.100, diff.bCHH.3x.100,
                 diff.cCpG.3x.100, diff.cCHG.3x.100, diff.cCHH.3x.100)
nameList <- list("aCpG", "aCHG", "aCHH", "bCpG", "bCHG", "bCHH","cCpG", "cCHG", "cCHH")
names(percList) <- nameList

wholeGenome1Mb <- Map(methDiffMatrix, obj=percList, name=nameList, binsize=1000000)

# Plot --------------------------------------------------------------------
percMatrix <- do.call(rbind, wholeGenome1Mb)
rownames(percMatrix) <- NULL
percMatrix$context <- factor(percMatrix$context, levels = c("CpG", "CHG", "CHH"))

# By chromosome
out <- lapply(1:10, plotDMAcrossChrom, percValues=percMatrix, plotOutput=T, binsize=1000000)
out2 <- plotDMAcrossChrom(2, percMatrix, plotOutput=F, rect=T, binsize=1000000)

# All genome
mop1pos <- 42155840/1000000+(percMatrix[percMatrix$seqnames == 2, "bin"][1]) #if you change binsize, you'll have to change me!
plotChrom<-ggplot(percMatrix) +
  geom_line(aes(bin, meth.diff, colour=context)) +
  labs(y="% DM") +
  ylim(-10,5) +
  scale_color_lancet(name="context") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank())  +
  theme(legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "horizontal") +
  geom_vline(xintercept = percMatrix[percMatrix$seqnames == 2, "bin"][1], colour="gray")+ # lines indicating chromosome boundaries
  geom_vline(xintercept = percMatrix[percMatrix$seqnames == 3, "bin"][1], colour="gray")+
  geom_vline(xintercept = percMatrix[percMatrix$seqnames == 4, "bin"][1], colour="gray")+
  geom_vline(xintercept = percMatrix[percMatrix$seqnames == 5, "bin"][1], colour="gray")+
  geom_vline(xintercept = percMatrix[percMatrix$seqnames == 6, "bin"][1], colour="gray")+
  geom_vline(xintercept = percMatrix[percMatrix$seqnames == 7, "bin"][1], colour="gray")+
  geom_vline(xintercept = percMatrix[percMatrix$seqnames == 8, "bin"][1], colour="gray")+
  geom_vline(xintercept = percMatrix[percMatrix$seqnames == 9, "bin"][1], colour="gray")+
  geom_vline(xintercept = percMatrix[percMatrix$seqnames == 10, "bin"][1], colour="gray")+
  geom_vline(xintercept = mop1pos, colour="yellow") +
  scale_x_continuous(breaks=c(percMatrix[percMatrix$seqnames == 2, "bin"][1]/2,
                              (percMatrix[percMatrix$seqnames == 2, "bin"][1]+percMatrix[percMatrix$seqnames == 3, "bin"][1])/2, # chromosome labels
                              (percMatrix[percMatrix$seqnames == 3, "bin"][1]+percMatrix[percMatrix$seqnames == 4, "bin"][1])/2,
                              (percMatrix[percMatrix$seqnames == 4, "bin"][1]+percMatrix[percMatrix$seqnames == 5, "bin"][1])/2,
                              (percMatrix[percMatrix$seqnames == 5, "bin"][1]+percMatrix[percMatrix$seqnames == 6, "bin"][1])/2,
                              (percMatrix[percMatrix$seqnames == 6, "bin"][1]+percMatrix[percMatrix$seqnames == 7, "bin"][1])/2,
                              (percMatrix[percMatrix$seqnames == 7, "bin"][1]+percMatrix[percMatrix$seqnames == 8, "bin"][1])/2,
                              (percMatrix[percMatrix$seqnames == 8, "bin"][1]+percMatrix[percMatrix$seqnames == 9, "bin"][1])/2,
                              (percMatrix[percMatrix$seqnames == 9, "bin"][1]+percMatrix[percMatrix$seqnames == 10, "bin"][1])/2,
                              (rev(percMatrix[percMatrix$seqnames == 10, "bin"])[1])-(rev(percMatrix[percMatrix$seqnames == 10, "bin"])[1]-rev(percMatrix[percMatrix$seqnames == 9, "bin"])[1])/2),
                     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),expand = c(0, 0))

aChrom <- plotChrom %+% subset(percMatrix, background %in% c("A")) + ggtitle("A") + theme(legend.position="none")
bChrom <- plotChrom %+% subset(percMatrix, background %in% c("B")) + ggtitle("B") + theme(legend.position="none")
cChrom <- plotChrom %+% subset(percMatrix, background %in% c("C")) + ggtitle("C") + theme(legend.position="none")

title.Chrom <- ggdraw() + 
  draw_label(paste0("Percent DM across chromosome "), fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
l.Chrom <- get_legend(plotChrom)
a1.Chrom <- plot_grid(aChrom, bChrom, cChrom, labels="AUTO", ncol=1)
a2.Chrom <- plot_grid(title.Chrom, a1.Chrom, l.Chrom, labels=NULL, ncol=1, rel_heights = c(0.1, 1))
save_plot(paste0("/home/nia/project/nia/results/3.7-DMAcrossChrom/DM-fullGenome.png"), a2.Chrom, base_asp = 2, ncol=1, nrow=3)

# ***Absolute methylation -------------------------------------------------
# Read in list
percList <- list("3.5-data/aMutCpG.rds","3.5-data/aWTCpG.rds","3.5-data/bMutCpG.rds","3.5-data/bWTCpG.rds","3.5-data/cMutCpG.rds","3.5-data/cWTCpG.rds",
                 "3.5-data/aMutCHG.rds","3.5-data/aWTCHG.rds","3.5-data/bMutCHG.rds","3.5-data/bWTCHG.rds","3.5-data/cMutCHG.rds","3.5-data/cWTCHG.rds",
                 "3.5-data/aMutCHH.rds","3.5-data/aWTCHH.rds","3.5-data/bMutCHH.rds","3.5-data/bWTCHH.rds","3.5-data/cMutCHH.rds","3.5-data/cWTCHH.rds")
nameList <- list("aMutCpG", "aWTCpG","bMutCpG", "bWTCpG","cMutCpG", "cWTCpG",
                 "aMutCHG", "aWTCHG","bMutCHG", "bWTCHG","cMutCHG", "cWTCHG",
                 "aMutCHH", "aWTCHH","bMutCHH", "bWTCHH","cMutCHH", "cWTCHH")
names(percList) <- nameList

#Get size of bins with which to divide chromosomes
binsize <- 1000

#Run function
percMatrix1k <- Map(function (x, y, z){
  print(y)
  
  #preprocess data
  raw <- readRDS(x)
  methyl = GRanges(seqnames(raw),
                   IRanges(start(raw), width=1),
                   perc.meth=mcols(raw)$perc.meth) #set width to 1 to avoid double counting
  chrLengths <- c(307041717,244442276,235667834,246994605,223902240,174033170,182381542,181122637,159769782,150982314)
  seqlengths(methyl) <- chrLengths # chromosome lengths
  
  #bin into even chunks
  bins <- IRangesList(lapply(seqlengths(methyl),
                             function(seqlen)
                               IRanges(breakInChunks(seqlen, z))))
  GRbins <- as(bins, "GRanges")
  
  #add NA to regions without methylation counts
  allGaps <- gaps(methyl)
  impGaps <- allGaps[strand(allGaps) == "*"] # remove start-to-finish ranges (appear as + and - strands)
  impGaps$perc.meth <- NA
  
  #ensure that there is no overlap between NAs and percent methylation
  if (length(findOverlaps(impGaps, methyl)) != 0) {
    print("OVERLAPS APPEAR, FIX ME")
  }
  
  #merge and sort ranges
  fullRange <- c(impGaps, methyl)
  fullRangeSort <- sort(fullRange)
  
  #calculate Rle and bin the averages for each window
  score <- GenomicRanges::coverage(fullRangeSort, weight="perc.meth")
  chromMat <- GenomicRanges::binnedAverage(GRbins, score, "perc.meth", na.rm=TRUE)
  chromDF <- as.data.frame(chromMat)
  chromDF$bin <- rep(1:1000, 10)
  chromDF$genotype <- y
  chromDF$background <- as.factor(toupper(substr(chromDF$genotype, 1, 1)))
  chromDF$status <- as.factor(substr(chromDF$genotype, 2, nchar(as.character(chromDF$genotype))-3))
  chromDF$context <- as.factor(substr(chromDF$genotype, nchar(as.character(chromDF$genotype))-2, nchar(as.character(chromDF$genotype))))
  outDF = subset(chromDF, select = -c(start, end, width, strand))
  
  return(outDF)
}, percList, nameList, binsize)

# Save matrix for later
saveRDS(percMatrix1k, "3.7-data/percMatrix.rds")


# Plot 
percMatrix1k <- readRDS("3.7-data/percMatrix.rds")

percMatrix <- do.call(rbind, percMatrix1k) # ADJUST ME AS NECESSARY - SUFFIX
rownames(percMatrix) <- NULL

percMatrix$status <- factor(percMatrix$status, levels = c("WT", "Mut"))
percMatrix$context <- factor(percMatrix$context, levels = c("CpG", "CHG", "CHH"))

centObj <- GRanges(seqnames=c("1","2", "3", "4", "5", "6", "7", "8", "9", "10"), 
                   ranges=IRanges(start=c(136770000,95510000,85780000,109070000,104540000,52300000,56380000,50530000,53750000,51390000),
                                  end=c(137120000, 97490000,86930000,110500000,106820000,53110000,56680000,52070000,57760000,52780000)))
chrLengths <- c(307041717,244442276,235667834,246994605,223902240,174033170,182381542,181122637,159769782,150982314)

for (i in 1:10) {
  newMat <- percMatrix[percMatrix$seqnames==i,]
  
  plotChrom<-ggplot(newMat) +
    geom_line(aes(bin, perc.meth, linetype=status, colour=context)) +
    labs(y="% mC") +
    ylim(0,100) +
    theme_classic() +
    scale_color_lancet(name="Context") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())  +
    theme(legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal"
    ) +
    scale_linetype_discrete(name="Methylation Status")  + 
    annotate("rect", xmin=start(centObj[i])/chrLengths[i]*1000, xmax=end(centObj[i])/chrLengths[i]*1000, ymin=0, ymax=Inf, alpha=0.8, fill="gray") 
  
  aChrom <- plotChrom %+% subset(newMat, background %in% c("A")) + ggtitle("A") + theme(legend.position="none")
  bChrom <- plotChrom %+% subset(newMat, background %in% c("B")) + ggtitle("B") + theme(legend.position="none")
  cChrom <- plotChrom %+% subset(newMat, background %in% c("C")) + ggtitle("C") + theme(legend.position="none")
  
  assign(paste("chr", i, ".a", sep=""),aChrom)
  assign(paste("chr", i, ".b", sep=""),bChrom)
  assign(paste("chr", i, ".c", sep=""),cChrom)
  
  title.Chrom <- ggdraw() + 
    draw_label(paste0("Percent methylation across chromosome ", i), fontface = 'bold', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  l.Chrom <- get_legend(plotChrom)
  a1.Chrom <- plot_grid(aChrom, bChrom, cChrom, labels="AUTO", ncol=1)
  a2.Chrom <- plot_grid(title.Chrom, a1.Chrom, l.Chrom, labels=NULL, ncol=1, rel_heights = c(0.1, 1))
  save_plot(paste0("/home/nia/project/nia/results/3.7-DMAcrossChrom/chr", i, ".png"), a2.Chrom, base_asp = 2, ncol=1, nrow=3)
}
