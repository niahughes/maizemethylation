###################################################################
### METHYLKIT BEDGRAPH GENERATION
### R script for bedgraph creation by methylKit: coverage, percent
### methylation per cytosine, differential methylation per 
### comparison
###################################################################

# set up environment
setwd("/scratch/nia/manuFinal")
load(file=".RData")
library(methylKit)
library(GenomicRanges)
library(genomation)


# Coverage/perc.meth ------------------------------------------------------
##NOTE this was not run again after subsetting; results were copied over from old results
##Subsetting was performed on data downstream from the input here and does not affect their values

print("Generating coverage and percent methylation bedgraphs")
fileNames<-c("aMut", "bMut", "cMut", "aWT", "bWT", "cWT")

inRaw <- list(CpG.raw, CpG.3x, CHG.raw, CHG.3x, CHH.raw, CHH.3x)
inNames <- list("raw.CpG", "3x.CpG", "raw.CHG", "3x.CHG", "raw.CHH", "3x.CHH")
names(inRaw) <- inNames

covOut <- Map(function(x,y) {
  print(y)
  for (i in 1:6){
    mypath1<-file.path("/home/nia/project/nia/results/bedgraphs/mKcoverage/", paste0(y, fileNames[i],".bedgraph"))
    bedgraph(x, file.name=mypath1, col.name="coverage") # output coverage bedgraph for each sample
    
    mypath2<-file.path("/home/nia/project/nia/results/bedgraphs/mKpercMeth/", paste(y, fileNames[i],".bedgraph", sep=""))
    bedgraph(x[[i]], file.name=mypath2, col.name="perc.meth") # output percent methylation bedgraph for each sample
  }
}, inRaw, inNames)


# All methylation differences ---------------------------------------------
print("Generating DM bedgraphs")
methCalls <- list(diff.CpG.3x.100.oChi, diff.aCpG.3x.100, diff.bCpG.3x.100, diff.cCpG.3x.100,
                  diff.CHG.3x.100.oChi, diff.aCHG.3x.100, diff.bCHG.3x.100, diff.cCHG.3x.100,
                  diff.CHH.3x.100.oChi, diff.aCHH.3x.100, diff.bCHH.3x.100, diff.cCHH.3x.100,
                  CpG.3x.100.oChi.50p.all, aCpG.3x.100.50p.all, bCpG.3x.100.50p.all, cCpG.3x.100.50p.all,
                  CHG.3x.100.oChi.50p.all, aCHG.3x.100.50p.all, bCHG.3x.100.50p.all, cCHG.3x.100.50p.all,
                  CHH.3x.100.oChi.50p.all, aCHH.3x.100.50p.all, bCHH.3x.100.50p.all, cCHH.3x.100.50p.all)
callNames <- list("CpG", "aCpG", "bCpG", "cCpG", "CHG", "aCHG", "bCHG", "cCHG", "CHH", "aCHH", "bCHH", "cCHH",
                  "sig-CpG", "sig-aCpG", "sig-bCpG", "sig-cCpG", "sig-CHG", "sig-aCHG", "sig-bCHG", "sig-cCHG", "sig-CHH", "sig-aCHH", "sig-bCHH", "sig-cCHH")
names(methCalls) <- callNames

out <- Map(function(x, y) {
  print(y)
  bedgraph(x, file.name=paste0("/home/nia/project/nia/results/bedgraphs/mKmethDiff/", y, ".bedgraph"), col.name="meth.diff") # output differential methylation bedgraph for each comparison
  return(NULL)
}, methCalls, callNames)

q(save="no")