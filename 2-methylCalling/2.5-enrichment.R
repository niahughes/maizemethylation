###################################################################
### ENRICHMENT TESTS OF DMR ANNOTATION
### R script to test if DMRs preferentially map to feature types
### using the chi-squared test
###################################################################

# set up environment
setwd("/scratch/nia/manuFinal")
load(file=".RData")
library(methylKit)
library(GenomicRanges)
library(genomation)
library(rtracklayer)
library(reshape2)
library(ggplot2)
library(ggsci)
library(cowplot)

# Enrichment tests --------------------------------------------------------
# list of all tiles with coverage in each comparison
bg.tiles <- list(diff.CpG.3x.100.oChi, diff.aCpG.3x.100, diff.bCpG.3x.100, diff.cCpG.3x.100,
                 diff.CHG.3x.100.oChi, diff.aCHG.3x.100, diff.bCHG.3x.100, diff.hCHG.3x.100,
                 diff.CHH.3x.100.oChi, diff.aCHH.3x.100, diff.bCHH.3x.100, diff.cCHH.3x.100)
bg.names <- list("fullCpG", "aCpG", "bCpG", "cCpG",
                 "fullCHG", "aCHG", "bCHG", "hCHG",
                 "fullCHH", "aCHH", "bCHH", "cCHH") 
names(bg.tiles) <- bg.names
# list of significantly DM tiles in each comparison
dm.tiles <- list(CpG.3x.100.oChi.50p.all, aCpG.3x.100.50p.all, bCpG.3x.100.50p.all, cCpG.3x.100.50p.all,
                 CHG.3x.100.oChi.50p.all, aCHG.3x.100.50p.all, bCHG.3x.100.50p.all, hCHG.3x.100.50p.all,
                 CHH.3x.100.oChi.50p.all, aCHH.3x.100.50p.all, bCHH.3x.100.50p.all, cCHH.3x.100.50p.all)
# annotation object
feature <- readRDS("3.8-data/regionObjs.rds")
feature$genesAndFlanks <- resize(feature$genes, width=width(feature$genes)+2000, fix="center")
featureType <- names(feature)

enrichment <- Map(function(n1, n2, n3){
  print(paste(n3, "-------------------------"))
  
  bg=as(n1, "GRanges")  ## All tiles with coverage in this model
  dm=as(n2, "GRanges") ## Tiles meeting DM cutoffs in this model
  hypo=dm[dm$meth.diff<0] ## Sig. hypomethylated tiles
  hyper=dm[dm$meth.diff>0] ## Sig. hypermethylated tiles
  
  A <- length(bg) ## Total number of tiles to be considered
  B <- length(hypo) ## Number of hypo tiles
  C <- length(hyper) ## Number of hyper tiles
  
  #Run tests for each feature type
  outList <- Map(function(feature, name) {
    print(name)
    
    bg.feat=mergeByOverlaps(bg, feature)
    hypo.feat=mergeByOverlaps(hypo, feature) ## Hypo tiles mapped to feature
    hyper.feat=mergeByOverlaps(hyper, feature) ## Hyper tiles mapped to feature
    
    D <- length(bg.feat$bg[!duplicated(bg.feat$bg)]) # Total number of unique tiles mapping to feature
    E <- length(hypo.feat$hypo[!duplicated(hypo.feat$hypo)]) # Number of hypo tiles mapped to feature
    F <- length(hyper.feat$hyper[!duplicated(hyper.feat$hyper)]) # Number of hyper tiles mapped to feature
    
    table.feat <- as.table(rbind(c(E, F, (D-(E+F))), c((B-E), (C-F), ((A-D)-((B-E)+(C-F))))))
    dimnames(table.feat) <- list(mapping=c("feature", "non-feature"),state=c("hypo", "hyper", "non-DM"))
    
    out.feat <- chisq.test(table.feat)
    
    return(out.feat)
    
  }, feature, featureType)
  
  return(outList)
  
}, bg.tiles, dm.tiles, bg.names)

saveRDS(enrichment, "/home/nia/project/nia/results/2.5-enrichment/enrichment.rds")


# Hypomethylation ---------------------------------------------------------
hypoEnr <- Map(function(n1, n2, n3){
  print(paste(n3, " hypo --------------------"))
  
  bg=as(n1, "GRanges")  ## All tiles with coverage in this model
  dm=as(n2, "GRanges") ## Tiles meeting DM cutoffs in this model
  hypo=dm[dm$meth.diff<0] ## Sig. hypomethylated tiles
  
  A <- length(bg) ## Total number of tiles to be considered
  B <- length(hypo) ## Number of hypo tiles
  
  #Run tests for each feature type
  outList <- Map(function(feature, name) {
    print(name)
    
    bg.feat=mergeByOverlaps(bg, feature)
    hypo.feat=mergeByOverlaps(hypo, feature) ## Hypo tiles mapped to feature
    
    D <- length(bg.feat$bg[!duplicated(bg.feat$bg)]) # Total number of unique tiles mapping to feature
    E <- length(hypo.feat$hypo[!duplicated(hypo.feat$hypo)]) # Number of hypo tiles mapped to feature
    
    table.feat <- as.table(rbind(c(E, (D-E)), c((B-E), ((A-B)-(D-E)))))
    dimnames(table.feat) <- list(mapping=c("feature", "non-feature"),state=c("hypo",  "non-hypo"))
    
    out.feat <- chisq.test(table.feat)
    
    return(out.feat)
    
  }, feature, featureType)
  
  return(outList)
  
}, bg.tiles, dm.tiles, bg.names)

saveRDS(hypoEnr, "/home/nia/project/nia/results/2.5-enrichment/hypoEnr.rds")


# Hypermethylation --------------------------------------------------------
hyperEnr <- Map(function(n1, n2, n3){
  print(paste(n3, " hyper --------------------"))
  
  bg=as(n1, "GRanges")  ## All tiles with coverage in this model
  dm=as(n2, "GRanges") ## Tiles meeting DM cutoffs in this model
  hyper=dm[dm$meth.diff<0] ## Sig. hypermethylated tiles
  
  A <- length(bg) ## Total number of tiles to be considered
  B <- length(hyper) ## Number of hyper tiles
  
  #Run tests for each feature type
  outList <- Map(function(feature, name) {
    print(name)
    
    bg.feat=mergeByOverlaps(bg, feature)
    hyper.feat=mergeByOverlaps(hyper, feature) ## hyper tiles mapped to feature
    
    D <- length(bg.feat$bg[!duplicated(bg.feat$bg)]) # Total number of unique tiles mapping to feature
    E <- length(hyper.feat$hyper[!duplicated(hyper.feat$hyper)]) # Number of hyper tiles mapped to feature
    
    table.feat <- as.table(rbind(c(E, (D-E)), c((B-E), ((A-B)-(D-E)))))
    dimnames(table.feat) <- list(mapping=c("feature", "non-feature"),state=c("hyper",  "non-hyper"))
    
    out.feat <- chisq.test(table.feat)
    
    return(out.feat)
    
  }, feature, featureType)
  
  return(outList)
  
}, bg.tiles, dm.tiles, bg.names)

saveRDS(hyperEnr, "/home/nia/project/nia/results/2.5-enrichment/hyperEnr.rds")

hypoEnr <- readRDS("/home/nia/project/nia/results/2.5-enrichment/hypoEnr.rds")

# Plot --------------------------------------------------------------------
methyl = GRanges(seqnames= c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
                 IRanges(start(rep(1,10)), 
                         end=c(307041717,244442276,235667834,246994605,223902240,174033170,182381542,181122637,159769782,150982314)))
methyl2 = tile(methyl, width=100)
methyl3 <- unlist(methyl2)
features=c(feature$genes, feature$promoters, feature$allTEs)
features2 <- reduce(features)
features2 <- features2[seqnames(features2) != "Pt"]
features2 <- features2[seqnames(features2) != "Mt"]
(sum(width(methyl))-(sum(width(features2)))/sum(width(methyl)))*100
sum(width(subsetByOverlaps( allObj,methyl)))

justGenes <- feature$genes
geneObj <- resize(feature$genes, width=width(feature$genes)+2000, fix="center")
promObj <- feature$promoters
retObj <- feature$classI
dnaObj <- feature$classII
allObj <- c(feature$genes, feature$allTEs, feature$promoters)
bg.names <- list("full model CpG", "A CpG", "B CpG", "C CpG",
                 "full model CHG", "A CHG", "B CHG", "C CHG",
                 "full model CHH", "A CHH", "B CHH", "C CHH") 
names(dm.tiles) <- bg.names

out <- lapply(dm.tiles, function(x) {
  y <- as(x, "GRanges")
  hypo <- y[y$meth.diff < 0]
  hyper <- y[y$meth.diff > 0]
  tiles <- list(hypo=hypo, hyper=hyper)
  out2 <- lapply(tiles, function(x) {
    
    
    gene <- length(subsetByOverlaps(x, justGenes))/length(x)
    prom <- length(subsetByOverlaps(x, promObj))/length(x)
    dna <- length(subsetByOverlaps(x, dnaObj))/length(x)
    ret <- length(subsetByOverlaps(x, retObj))/length(x)
    none <- length(subsetByOverlaps(x, allObj, invert=TRUE))/length(x)
    
    outList <- list(gene=gene, prom=prom, ret=ret, dna=dna, none=none)
    
    return(outList)
  })
  return(out2)
})


outMelt <- melt(out)
outMelt$context <- factor(substr(outMelt$L1, nchar(outMelt$L1)-2, nchar(outMelt$L1)), levels=c("CpG", "CHG", "CHH"))
outMelt$model <- factor(substr(outMelt$L1, 1, nchar(outMelt$L1)-4), levels=c("full model", "A", "B", "C"))
outMelt$L3 <- factor(outMelt$L3, levels=c("gene", "prom", "ret", "dna", "none"))
outMelt <- subset(outMelt, model !="full model")

ggplot(outMelt, aes(x = context, y = value, fill = L3)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_grid(~model) +
  scale_fill_lancet(name="Feature", labels=c("Genes", "Promoters", "Retrotransposons", "DNA transposons", "No feature")) +
  labs(title="Percentage of DM tiles mapping to features", 
       x ="Cytosine sequence context", y = "% of tiles")

hypo <- ggplot(outMelt[outMelt$L2 == "hypo",], aes(x = context, y = value, fill = L3)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_grid(~model) +
  scale_fill_lancet() +
  labs(title="Percentage of hypomethylated tiles mapping to features", 
       x =NULL, y = "% of tiles") +
  theme(legend.position = "none")

hyper <- ggplot(outMelt[outMelt$L2 == "hyper",], aes(x = context, y = value, fill = L3)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_grid(~model)+
  scale_fill_lancet(name="Feature", labels=c("Genes", "Promoters", "Retrotransposons", "DNA transposons", "No feature")) +
  labs(title="Percentage of hypermethylated tiles mapping to features", 
       x ="Cytosine sequence context", y = "% of tiles") +
  theme(legend.direction = "horizontal", legend.position = "bottom")


all <- ggplot(outMelt2, aes(x=all, y=value, fill=L3)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_lancet(name="Feature", labels=c("Genes", "Promoters", "Retrotransposons", "DNA transposons", "No feature"))+ 
  theme(legend.position="none") +
  labs(title="Genome", x ="All tiles", y = "% of tiles") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.text.x=element_blank())

l <- get_legend(hyper)

a <- plot_grid(
  hypo, hyper + theme(legend.position="none"),
  labels = "AUTO", ncol = 1
)
c <- plot_grid(a, all, ncol=2, rel_widths = c(4, 1), labels = c('', 'C'))
d <- plot_grid(c, l, ncol=1, rel_heights = c(1, 0.05))
d
