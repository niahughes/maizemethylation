###################################################################
### METHYLKIT ANALYSES AND DIAGRAMS
### R script for methylKit comparative analyses and diagrams 
### Differential methylation by chromosome script adapted from 
### methylKit source code with help from Delvin So
###################################################################

# set up environment
setwd("/scratch/nia/manuFinal")
load(file=".RData")
library(methylKit)
library(GenomicRanges)
library(genomation)
library(reshape2)
library(ggplot2)


# Histograms --------------------------------------------------------------
##NOTE THAT THIS WAS NOT RE-RUN - histograms in target folder were copied from old results 
##They were generated before the subset, so were therefore unaffected

names<-c("A mutant", "B mutant", "C mutant", "A wild-type", "B wild-type", "C wild-type")
fileNames<-c("aMut", "bMut", "cMut", "aWT", "bWT", "cWT")

inRaw <- list(CpG.raw, CpG.3x, CHG.raw, CHG.3x, CHH.raw, CHH.3x)
inNames <- list("raw.CpG", "3x.CpG", "raw.CHG", "3x.CHG", "raw.CHH", "3x.CHH")
names(inRaw) <- inNames

covHist <- Map(function(x,y) { # coverage histograms
  print(y)
  for (i in 1:6){
    print(fileNames[i])
    mypath1<-file.path("/home/nia/project/nia/results/3.1-methylKit/histograms/", 
                      paste(y,"-Coverage-",fileNames[i],".pdf", sep=""))
    
    pdf(file=mypath1)
    getCoverageStats(x[[i]], plot=TRUE, both.strands=FALSE)
    dev.off()
    
    mypath2<-file.path("/home/nia/project/nia/results/3.1-methylKit/histograms/",
                       paste(y,"-Coverage-", fileNames[i], ".txt", sep=""))
    
    sink(mypath2,append=TRUE)
    cat(names[i], " ", y, " coverage stats\n")
    getCoverageStats(x[[i]], plot=FALSE, both.strands=FALSE)
    sink()
    
  }
}, inRaw, inNames)

methHist <- Map(function(x,y) { # methylation percentage histograms
  for (i in 1:6) {
    mypath1<-file.path("/home/nia/project/nia/results/3.1-methylKit/histograms/", 
                       paste(y,"-Methylation-",fileNames[i],".pdf", sep=""))
    
    pdf(file=mypath1)
    getMethylationStats(x[[i]], plot=TRUE, both.strands=FALSE)
    dev.off()
  
    mypath2<-file.path("/home/nia/project/nia/results/3.1-methylKit/histograms/",
                       paste(y,"-Methylation-", fileNames[i], ".txt", sep=""))
    
    sink(mypath2,append=TRUE)
    cat(names[i], " ", y, " raw CpG methylation stats\n")
    getMethylationStats(x[[i]], plot=FALSE, both.strands=FALSE)
  }
}, inRaw, inNames)


# Comparisons -------------------------------------------------------------
# comparative analyses on sample methylation profiles
compIn <- list(meth.CpG.3x.100, meth.CHG.3x.100, meth.CHH.3x.100)
compNames <- list("CpG", "CHG", "CHH")

compareSamples <- Map(function(x,y) {
  pdf(file=paste0("/home/nia/project/nia/results/3.1-methylKit/correlation/",y,"Corr.pdf")) # correlation plots
  getCorrelation(x, plot=TRUE)
  dev.off()
  
  pdf(file=paste0("/home/nia/project/nia/results/3.1-methylKit/clustering/",y,"Cluster.pdf")) # clustering
  clusterSamples(x, dist="correlation", method="ward", plot=TRUE)
  dev.off()
  
  pdf(file=paste0("/home/nia/project/nia/results/3.1-methylKit/PCA/",y,"PCA.pdf")) # PCA
  PCASamples(x)
  dev.off()
  
  pdf(file=paste0("/home/nia/project/nia/results/3.1-methylKit/PCA/",y,"Scree.pdf")) # screeplot
  PCASamples(x, screeplot=TRUE)
  dev.off()
}, compIn, compNames)

# DM by chromosome --------------------------------------------------------

diffList <- list(diff.CpG.3x.100.oChi,diff.aCpG.3x.100,diff.bCpG.3x.100,diff.cCpG.3x.100,
                 diff.CHG.3x.100.oChi,diff.aCHG.3x.100,diff.bCHG.3x.100,diff.cCHG.3x.100,
                 diff.CHH.3x.100.oChi,diff.aCHH.3x.100,diff.bCHH.3x.100,diff.cCHH.3x.100)
diffNames <- list("full CpG", "A CpG", "B CpG", "C CpG",
                  "full CHG", "A CHG", "B CHG", "C CHG",
                  "full CHH", "A CHH", "B CHH", "C CHH")
names(diffList) <- diffNames

diffList <- list(diff.aCpG.3x.100,diff.bCpG.3x.100,diff.cCpG.3x.100,
                 diff.aCHG.3x.100,diff.bCHG.3x.100,diff.cCHG.3x.100,
                 diff.aCHH.3x.100,diff.bCHH.3x.100,diff.cCHH.3x.100)
diffNames <- list("A CpG", "B CpG", "C CpG",
                  "A CHG", "B CHG", "C CHG",
                  "A CHH", "B CHH", "C CHH")
names(diffList) <- diffNames

qvalue.cutoff <- 0.01
meth.cutoff <- 50

# plot proportion of tiles on each chromosome considered significantly DM
dMBC <- Map(function(x, y, createPlot, returnPlot){
  print(y)
  x=getData(x)
  temp.hyper=x[x$qvalue < qvalue.cutoff & x$meth.diff >= meth.cutoff,]
  temp.hypo =x[x$qvalue < qvalue.cutoff & x$meth.diff <= -meth.cutoff,]
  
  dmc.hyper=100*nrow(temp.hyper)/nrow(x) # get percentages of hypo/ hyper
  dmc.hypo =100*nrow(temp.hypo )/nrow(x)
  
  all.hyper.hypo=data.frame(number.of.hypermethylated=nrow(temp.hyper),
                            percentage.of.hypermethylated=dmc.hyper,
                            number.of.hypomethylated=nrow(temp.hypo),
                            percentage.of.hypomethylated=dmc.hypo)
  
  # plot barplot for percentage of DMCs per chr
  dmc.hyper.chr=merge(as.data.frame(table(temp.hyper$chr)),
                      as.data.frame(table(x$chr)),by="Var1")
  dmc.hyper.chr=cbind(dmc.hyper.chr,
                      perc=100*dmc.hyper.chr[,2]/dmc.hyper.chr[,3])
  
  dmc.hypo.chr=merge(as.data.frame(table(temp.hypo$chr)),
                     as.data.frame(table(x$chr)),by="Var1")
  dmc.hypo.chr=cbind(dmc.hypo.chr,
                     perc=100*dmc.hypo.chr[,2]/dmc.hypo.chr[,3])
  
  # merge hyper hypo per chromosome
  dmc.hyper.hypo=merge(dmc.hyper.chr[,c(1,2,4)],
                       dmc.hypo.chr[,c(1,2,4)],by="Var1", all=TRUE)
  dmc.hyper.hypo[is.na(dmc.hyper.hypo)] <- 0
  dmc.hyper.hypo=dmc.hyper.hypo[order(as.numeric(sub("chr","",dmc.hyper.hypo$Var1))),] # order the chromosomes
  
  names(dmc.hyper.hypo)=c("chr","number.of.hypermethylated",
                          "percentage.of.hypermethylated",
                          "number.of.hypomethylated",
                          "percentage.of.hypomethylated")
  te <- melt(dmc.hyper.hypo, id.vars = c("chr", "number.of.hypermethylated", "number.of.hypomethylated"))
  te$chr_f <- factor(te$chr, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
  
  plot <- ggplot(data=te, aes(x=chr_f, y=value, fill=variable)) + 
    geom_bar(stat="identity") + 
    labs(title=paste0("Percent DM tiles in ", sub(".*? ", "", y), " context, ", sub(" .*", "", y), " comparison"), 
         x="Chromosome", y="% DM tiles") +
    scale_fill_manual(name="Methylation\nlevel",labels = c("Hypermethylated", "Hypomethylated"), values=c("#2b83ba", "#d7191c")) + 
    ylim(0, 0.6) +
    coord_flip() +
    theme(legend.direction = "horizontal", legend.position = "bottom")
  
  if (createPlot==TRUE) {
    png(paste0("/home/nia/project/nia/results/3.1-methylKit/DMByChrom/",gsub(" ", "", y, fixed = TRUE),"-dMBC.png"))
    print(plot)
    dev.off()
  }
  
  if (returnPlot==TRUE) {
    return(plot)
  } else {
    return(dmc.hyper.hypo)
  }
}, diffList, diffNames, createPlot=FALSE, returnPlot=TRUE)