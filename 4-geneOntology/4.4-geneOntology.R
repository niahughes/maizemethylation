###################################################################
### ONTOLOGY ANALYSIS OF DIFFERENTIALLY METHYLATED GENES
### Local R script 
###################################################################

# set up environment
setwd("D:/Nia/OneDrive - University of Guelph/Working/DMRanalysis")
source("D:/Nia/OneDrive - University of Guelph/Working/thesisScripts/4-geneOntology/4.1-dataPrep.R")
library(topGO)
library(Vennerable)
library(dplyr)
library(rtracklayer)

# Functions ---------------------------------------------------------------
# Convert mappings
getMappings <- function(inMap, outName) { # convert agriGO-provided mappings to topGO-compatible format
  
  GOraw <- read.table(inMap, sep="\t", header=FALSE, col.names=c("discard", "discard", "geneID", "GOterms"), row.names = NULL)[ ,3:4]
  GOmap <- GOraw %>% group_by(geneID) %>% summarise(GOterms=paste(GOterms, collapse=","))
  write.table(GOmap, file=paste0(outName, "-full.txt"), sep="\t",quote=FALSE, row.names=FALSE,col.names = FALSE) # output full mapping
  
  out <- anti_join(GOmap, namesToRemove, by="geneID")
  write.table(out, file=paste0(outName, "-noIntrogression.txt"), sep="\t",quote=FALSE, row.names=FALSE,col.names = FALSE) # output mapping without chr2 region
  
  return(NULL)
}

# gene ontology
geneOnt <- function(DMNames, mapping) {
  
  DM.list <- factor(as.integer(gamerNames %in% unique(DMNames))) #generate named list
  names(DM.list) <- gamerNames
  
  out <- lapply(ont, function(x) {
    data <- new("topGOdata", ontology=x, allGenes=DM.list, annot = annFUN.gene2GO, gene2GO = mapping, nodeSize=5) # construct topGO object
    result <- runTest(data, algorithm="weight01", statistic="fisher") #run GO analysis
    table <- GenTable(data, weight=result, topNodes=250) #generate table of top 250 terms
    
    coll <- list(data,result,table)
    names(coll) <- c("data", "result", "table")
    return(coll)
    
  })
  
  return(out)
}

geneOntElim <- function(DMNames, mapping) {
  
  DM.list <- factor(as.integer(gamerNames %in% unique(DMNames))) #generate named list
  names(DM.list) <- gamerNames
  
  out <- lapply(ont, function(x) {
    data <- new("topGOdata", ontology=x, allGenes=DM.list, annot = annFUN.gene2GO, gene2GO = mapping, nodeSize=5) # construct topGO object
    result <- runTest(data, algorithm="elim", statistic="fisher") #run GO analysis
    table <- GenTable(data, weight=result, topNodes=250) #generate table of top 250 terms
    
    coll <- list(data,result,table)
    names(coll) <- c("data", "result", "table")
    return(coll)
    
  })
  
  return(out)
}

# output significant tables
tableOut <- function(results, name, context, cutoffVal) {
  write.table(subset(results$BP$table[, c(1,6,3,4,5,2)], as.numeric(results$BP$table$weight)<cutoffVal), 
              file = paste0("4-Results/4.4-GO/", context, "/BP.", context, ".",name,"-", cutoffVal, ".txt"), sep="\t", row.names = F, quote = F)
  write.table(subset(results$CC$table[, c(1,6,3,4,5,2)], as.numeric(results$CC$table$weight)<cutoffVal), 
              file = paste0("4-Results/4.4-GO/", context, "/CC.", context, ".",name,"-", cutoffVal, ".txt"), sep="\t", row.names = F, quote = F)
  write.table(subset(results$MF$table[, c(1,6,3,4,5,2)], as.numeric(results$MF$table$weight)<cutoffVal), 
              file = paste0("4-Results/4.4-GO/", context, "/MF.", context, ".",name,"-", cutoffVal, ".txt"), sep="\t", row.names = F, quote = F)
}

# output diagram
diagOut <- function(results, name, context, numNodes) {
  pdf(file=paste0("4-Results/4.4-GO/",context,"diag/BP.", context, ".", name,".pdf"))
    showSigOfNodes(results$BP$data, score(results$BP$result), firstSigNodes=numNodes, useInfo='all')
  dev.off()
  
  pdf(file=paste0("4-Results/4.4-GO/",context,"diag/MF.", context, ".", name,".pdf"))
    showSigOfNodes(results$MF$data, score(results$MF$result), firstSigNodes=numNodes, useInfo='all')
  dev.off()
  
  pdf(file=paste0("4-Results/4.4-GO/",context,"diag/CC.", context, ".", name,".pdf"))
    showSigOfNodes(results$CC$data, score(results$CC$result), firstSigNodes=numNodes, useInfo='all')
  dev.off()
}

# Preparing GO mappings ---------------------------------------------------
## run once

#Download maize-GAMER and Ensembl mappings from agriGO v2 - included for posterity
download.file(url='http://systemsbiology.cau.edu.cn/agriGOv2/download/940_slimGO', destfile = "940_slimGO.txt") # GAMER
download.file(url="http://systemsbiology.cau.edu.cn/agriGOv2/download/958_slimGO", destfile = "958_slimGO.txt") # Ensembl

geneObj=import("D:/Nia/IGV/geneAnnotationb73.bed", format="BED")
removalRegion <- GRanges(seqnames=c("2"), ranges=IRanges(start=c(22000000), end=c(182000000)), strand="*")
genesToRemove <- subsetByOverlaps(geneObj, removalRegion)
namesToRemove <- data.frame(geneID=unique(gsub("transcript:|_.*$", "", genesToRemove$name)))

getMappings(inMap="940_slimGO.txt", outName="gamerMap")
getMappings(inMap="958_slimGO.txt", outName="ensemblMap")

# Run gene ontologies -----------------------------------------------------
# Get appropriate objects
gamerMap<-readMappings(file="gamerMap-noIntrogression.txt", sep = "\t", IDsep = ",") 
gamerNames<-names(gamerMap) 

ensemblMap<-readMappings(file="ensemblMap-noIntrogression.txt", sep = "\t", IDsep = ",")
ensemblNames<-names(ensemblMap)

names <- c("f.all", "a.all", "b.all", "c.all", "f.hypo", "a.hypo", "b.hypo", "c.hypo", "f.hyper", "a.hyper", "b.hyper", "c.hyper")
ont <- c(BP="BP", MF="MF", CC="CC")

# CpG ---------------------------------------------------------------------
## Each comparison individually
CpG.list <- list(complete.CpG$name, complete.aCpG$name, complete.bCpG$name, complete.cCpG$name,
                hyp.CpG$hypo, hyp.aCpG$hypo, hyp.bCpG$hypo, hyp.cCpG$hypo,
                hyp.CpG$hyper, hyp.aCpG$hyper, hyp.bCpG$hyper, hyp.cCpG$hyper)
names(CpG.list) <- names

CpG.gamer <- lapply(CpG.list, geneOnt, gamerMap)
Map(tableOut, CpG.gamer, names, "CpG", 0.005)
# Map(diagOut, CpG.gamer, names, "CpG", 5)
saveRDS(CpG.gamer, file = "CpGgamer.rds")

## Create Venn for various subsets/combinations
allVennCpG <- Venn(list(fullAll=complete.CpG$name, AAll=complete.aCpG$name, BAll=complete.bCpG$name, CAll=complete.cCpG$name))
hypoVennCpG <- Venn(list(fullHypo=hyp.CpG$hypo,AHypo=hyp.aCpG$hypo,BHypo=hyp.bCpG$hypo,CHypo=hyp.cCpG$hypo))
hyperVennCpG <- Venn(list(fullhyper=hyp.CpG$hyper,Ahyper=hyp.aCpG$hyper,Bhyper=hyp.bCpG$hyper,Chyper=hyp.cCpG$hyper))

## Solo
CpG.solo <- list(solo.a.all=allVennCpG@IntersectionSets$`0100`, solo.b.all=allVennCpG@IntersectionSets$`0010`, solo.c.all=allVennCpG@IntersectionSets$`0001`,
                 solo.a.hypo=hypoVennCpG@IntersectionSets$`0100`, solo.b.hypo=hypoVennCpG@IntersectionSets$`0010`, solo.c.hypo=hypoVennCpG@IntersectionSets$`0001`,
                 solo.a.hyper=hyperVennCpG@IntersectionSets$`0100`, solo.b.hyper=hyperVennCpG@IntersectionSets$`0010`, solo.c.hyper=hyperVennCpG@IntersectionSets$`0001`)
solo.CpG.gamer <- lapply(CpG.solo, geneOnt, gamerMap)
Map(tableOut, solo.CpG.gamer, names(CpG.solo), "CpG", 0.005)
saveRDS(solo.CpG.gamer, file = "soloCpGgamer.rds")

# CHG ---------------------------------------------------------------------
## Each comparison individually
CHG.list <- list(complete.CHG$name, complete.aCHG$name, complete.bCHG$name, complete.cCHG$name,
                 hyp.CHG$hypo, hyp.aCHG$hypo, hyp.bCHG$hypo, hyp.cCHG$hypo,
                 hyp.CHG$hyper, hyp.aCHG$hyper, hyp.bCHG$hyper, hyp.cCHG$hyper)
names(CHG.list) <- names

CHG.gamer <- lapply(CHG.list, geneOnt, gamerMap)
Map(tableOut, CHG.gamer, names, "CHG", 0.005)
# Map(diagOut, CHG.gamer, names, "CHG", 5)
saveRDS(CHG.gamer, file = "CHGgamer.rds")

## Create Venn for various subsets/combinations
allVennCHG <- Venn(list(fullAll=complete.CHG$name, AAll=complete.aCHG$name, BAll=complete.bCHG$name, CAll=complete.cCHG$name))
hypoVennCHG <- Venn(list(fullHypo=hyp.CHG$hypo,AHypo=hyp.aCHG$hypo,BHypo=hyp.bCHG$hypo,CHypo=hyp.cCHG$hypo))
hyperVennCHG <- Venn(list(fullhyper=hyp.CHG$hyper,Ahyper=hyp.aCHG$hyper,Bhyper=hyp.bCHG$hyper,Chyper=hyp.cCHG$hyper))

## Solo
CHG.solo <- list(solo.a.all=allVennCHG@IntersectionSets$`0100`, solo.b.all=allVennCHG@IntersectionSets$`0010`, solo.c.all=allVennCHG@IntersectionSets$`0001`,
                 solo.a.hypo=hypoVennCHG@IntersectionSets$`0100`, solo.b.hypo=hypoVennCHG@IntersectionSets$`0010`, solo.c.hypo=hypoVennCHG@IntersectionSets$`0001`,
                 solo.a.hyper=hyperVennCHG@IntersectionSets$`0100`, solo.b.hyper=hyperVennCHG@IntersectionSets$`0010`, solo.c.hyper=hyperVennCHG@IntersectionSets$`0001`)
solo.CHG.gamer <- lapply(CHG.solo, geneOnt, gamerMap)
Map(tableOut, solo.CHG.gamer, names(CHG.solo), "CHG", 0.005)
saveRDS(solo.CHG.gamer, file = "soloCHGgamer.rds")

# CHH ---------------------------------------------------------------------
## Each comparison individually
CHH.list <- list(complete.CHH$name, complete.aCHH$name, complete.bCHH$name, complete.cCHH$name,
                 hyp.CHH$hypo, hyp.aCHH$hypo, hyp.bCHH$hypo, hyp.cCHH$hypo,
                 hyp.CHH$hyper, hyp.aCHH$hyper, hyp.bCHH$hyper, hyp.cCHH$hyper)
names(CHH.list) <- names

CHH.gamer <- lapply(CHH.list, geneOnt, gamerMap)
Map(tableOut, CHH.gamer, names, "CHH", 0.005)
Map(diagOut, CHH.gamer, names, "CHH", 5)
saveRDS(CHH.gamer, file = "CHHgamer.rds")

## Create Venn for various subsets/combinations
allVennCHH <- Venn(list(fullAll=complete.CHH$name, AAll=complete.aCHH$name, BAll=complete.bCHH$name, CAll=complete.cCHH$name))
hypoVennCHH <- Venn(list(fullHypo=hyp.CHH$hypo,AHypo=hyp.aCHH$hypo,BHypo=hyp.bCHH$hypo,CHypo=hyp.cCHH$hypo))
hyperVennCHH <- Venn(list(fullhyper=hyp.CHH$hyper,Ahyper=hyp.aCHH$hyper,Bhyper=hyp.bCHH$hyper,Chyper=hyp.cCHH$hyper))

## Solo
CHH.solo <- list(solo.a.all=allVennCHH@IntersectionSets$`0100`, solo.b.all=allVennCHH@IntersectionSets$`0010`, solo.c.all=allVennCHH@IntersectionSets$`0001`,
                 solo.a.hypo=hypoVennCHH@IntersectionSets$`0100`, solo.b.hypo=hypoVennCHH@IntersectionSets$`0010`, solo.c.hypo=hypoVennCHH@IntersectionSets$`0001`,
                 solo.a.hyper=hyperVennCHH@IntersectionSets$`0100`, solo.b.hyper=hyperVennCHH@IntersectionSets$`0010`, solo.c.hyper=hyperVennCHH@IntersectionSets$`0001`)
solo.CHH.gamer <- lapply(CHH.solo, geneOnt, gamerMap)
Map(tableOut, solo.CHH.gamer, names(CHH.solo), "CHH", 0.005)
saveRDS(solo.CHH.gamer, file = "soloCHHgamer.rds")