###################################################################
### HILBERT CURVES OF DIFFERENTIAL METHYLATION ACROSS GENOME
### R script 
###################################################################

# set up environment
setwd("/scratch/nia/manuFinal")
load(file=".RData")
library(methylKit)
library(GenomicRanges)
library(genomation)
library(reshape2)
library(ggplot2)
library(dplyr)
library(cowplot)
library(HilbertCurve)
library(circlize)
library(ComplexHeatmap)

# Function -----------------------------------------------------------------
getHil <- function (x, y, sig){
  print(y)
  x <- as(x, "GRanges")
  x <- renameSeqlevels(x, c("1"="chr1", "2"="chr2", "3"="chr3", "4"="chr4", "5"="chr5","6"="chr6", "7"="chr7", "8"="chr8", "9"="chr9", "10"="chr10"))
  
  if (sig == FALSE) {
    png(file=paste0("/home/nia/project/nia/results/3.6-hilbertCurve/new9-", gsub("\\s", "", y),".png"), width=1000, height=1000)
  } else {
    png(file=paste0("/home/nia/project/nia/results/3.6-hilbertCurve/new9-sig-", gsub("\\s", "", y),".png"), width=1000, height=1000)
  }
  hc = GenomicHilbertCurve(chr = paste0("chr",1:10), level = 9, mode = "pixel", title=paste0(y, " differential methylation"), legend = list(lgd))
  hc_layer(hc, x, col = col_fun(x$meth.diff), mean_mode = "absolute")
  hc_map(hc,add = TRUE, fill = NA, border = "grey", labels=1:10)
  dev.off()
  
  #Centromere overlay
  if (sig == FALSE) {
    png(file=paste0("/home/nia/project/nia/results/3.6-hilbertCurve/new9-cent-", gsub("\\s", "", y),".png"), width=1000, height=1000)
  } else {
    png(file=paste0("/home/nia/project/nia/results/3.6-hilbertCurve/new9-sigcent-", gsub("\\s", "", y),".png"), width=1000, height=1000)
  }
  hc = GenomicHilbertCurve(chr = paste0("chr",1:10), level = 9, mode = "pixel", title=paste0(y, " differential methylation"), legend = list(lgd, lgdCent))
  hc_layer(hc, x, col = col_fun(x$meth.diff), mean_mode = "absolute")
  hc_layer(hc, reduce(centObj), col = "#00000020")
  hc_map(hc,add = TRUE, fill = NA, border = "grey", labels=1:10)
  dev.off()
  
  return(NULL)
}

# Object prep -------------------------------------------------------------
print("Getting necessary objects")
# centromere positions
centObj <- GRanges(seqnames=c("chr1","chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr9", "chr10"), #
                   ranges=IRanges(start=c(136770000,95510000,85780000,109070000,104540000,52300000,56380000,50530000,53750000,57360000,51390000),
                                  end=c(137120000, 97490000,86930000,110500000,106820000,53110000,56680000,52070000,55390000,57760000,52780000)))
# legend setup
col_fun = colorRamp2(c(-100, 0, 100), c("red", "white", "blue"))
lgd = Legend(col_fun = col_fun, title = "Differential\nMethylation")
# centromere overlay
lgdCent = Legend(labels = c("Centromere"), legend_gp = gpar(fill = c("#00000020")), title = NULL)
# all DM levels as list
diffList <- list(diff.CpG.3x.100.oChi,diff.aCpG.3x.100,diff.bCpG.3x.100,diff.cCpG.3x.100,
                 diff.CHG.3x.100.oChi,diff.aCHG.3x.100,diff.bCHG.3x.100,diff.cCHG.3x.100,
                 diff.CHH.3x.100.oChi,diff.aCHH.3x.100,diff.bCHH.3x.100,diff.cCHH.3x.100)
diffNames <- list("Full model CpG", "A CpG", "B CpG", "C CpG",
                  "Full model CHG", "A CHG", "B CHG", "C CHG",
                  "Full model CHH", "A CHH", "B CHH", "C CHH")
names(diffList) <- diffNames
# only significantly DM levels as list
sigdiffList <- list(CpG.3x.100.oChi.50p.all,aCpG.3x.100.50p.all,bCpG.3x.100.50p.all,cCpG.3x.100.50p.all,
                    CHG.3x.100.oChi.50p.all,aCHG.3x.100.50p.all,bCHG.3x.100.50p.all,cCHG.3x.100.50p.all,
                    CHH.3x.100.oChi.50p.all,aCHH.3x.100.50p.all,bCHH.3x.100.50p.all,cCHH.3x.100.50p.all)
sigdiffNames <- list("Full model CpG", "A CpG", "B CpG", "C CpG",
                     "Full model CHG", "A CHG", "B CHG", "C CHG",
                     "Full model CHH", "A CHH", "B CHH", "C CHH")
names(sigdiffList) <- sigdiffNames

# Run ---------------------------------------------------------------------
Map(getHil, x=diffList, y=diffNames, sig=FALSE) # all DM levels
Map(getHil, x=sigdiffList, y=sigdiffNames, sig=TRUE) # only significantly DM tiles
