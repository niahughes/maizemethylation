###################################################################
### QUANTIFYING ANNOTATED TILES
### Local R script 
###################################################################

# set up environment
setwd("D:/Nia/OneDrive - University of Guelph/Working/DMRanalysis")
source("D:/Nia/OneDrive - University of Guelph/Working/thesisScripts/4-geneOntology/4.1-dataPrep.R")
library(ggplot2)

# CpG ---------------------------------------------------------------------

# Quantifying raw numbers -------------------------------------------------
# Feature counts by genotype, all hits
CpG.tile <- data.frame (
  comp=c(rep("full",6), rep("A",6), rep("B",6), rep("C",6)),
  feature=c(rep("promoter",2), rep("exon",2), rep("intron",2)),
  level=c("hyper","hypo"),
  amount=c(length(which(tile.CpG$promoters$meth.diff>0)), length(which(tile.CpG$promoters$meth.diff<0)), 
           length(which(tile.CpG$exons$meth.diff>0)), length(which(tile.CpG$exons$meth.diff<0)),
           length(which(tile.CpG$introns$meth.diff>0)), length(which(tile.CpG$introns$meth.diff<0)),
           length(which(tile.aCpG$promoters$meth.diff>0)), length(which(tile.aCpG$promoters$meth.diff<0)), 
           length(which(tile.aCpG$exons$meth.diff>0)), length(which(tile.aCpG$exons$meth.diff<0)),
           length(which(tile.aCpG$introns$meth.diff>0)), length(which(tile.aCpG$introns$meth.diff<0)),
           length(which(tile.bCpG$promoters$meth.diff>0)), length(which(tile.bCpG$promoters$meth.diff<0)), 
           length(which(tile.bCpG$exons$meth.diff>0)), length(which(tile.bCpG$exons$meth.diff<0)),
           length(which(tile.bCpG$introns$meth.diff>0)), length(which(tile.bCpG$introns$meth.diff<0)),
           length(which(tile.cCpG$promoters$meth.diff>0)), length(which(tile.cCpG$promoters$meth.diff<0)), 
           length(which(tile.cCpG$exons$meth.diff>0)), length(which(tile.cCpG$exons$meth.diff<0)),
           length(which(tile.cCpG$introns$meth.diff>0)), length(which(tile.cCpG$introns$meth.diff<0)))
)
CpG.tile$comp_f=factor(CpG.tile$comp, levels=c("full", "A", "B", "C"))
CpG.tile$feature_f=factor(CpG.tile$feature, levels=c("promoter", "exon", "intron"))

pdf(file="4-Results/4.2-quant/CpG-1-AllTilesByFeat.pdf")
ggplot(data=CpG.tile, aes(x=feature_f, y=amount, fill=level)) + 
  geom_bar(stat="identity") + 
  facet_grid(~comp_f) +
  labs(title="Number of DM tiles in CpG context by comparison and feature", x="Feature type", y="Number of DM tiles") +
  scale_fill_manual(values=c("#2b83ba", "#d7191c"))
dev.off()

# Promoter vs. genic, all hits
CpG.allhits <- data.frame (
  comp=c(rep("full",4), rep("A",4), rep("B",4), rep("C",4)),
  feature=c(rep("promoter",2), rep("genic",2)),
  level=c("hyper","hypo"),
  amount=c(nrow(split3.CpG$hyperPromoters), nrow(split3.CpG$hypoPromoters), 
           nrow(split3.CpG$hyperGenic), nrow(split3.CpG$hypoGenic),
           nrow(split3.aCpG$hyperPromoters), nrow(split3.aCpG$hypoPromoters), 
           nrow(split3.aCpG$hyperGenic), nrow(split3.aCpG$hypoGenic),
           nrow(split3.bCpG$hyperPromoters), nrow(split3.bCpG$hypoPromoters), 
           nrow(split3.bCpG$hyperGenic), nrow(split3.bCpG$hypoGenic),
           nrow(split3.cCpG$hyperPromoters), nrow(split3.cCpG$hypoPromoters), 
           nrow(split3.cCpG$hyperGenic), nrow(split3.cCpG$hypoGenic))
)
CpG.allhits$comp_f=factor(CpG.allhits$comp, levels=c("full", "A", "B", "C"))
CpG.allhits$feature_f=factor(CpG.allhits$feature, levels=c("promoter", "genic"))

pdf(file="4-Results/4.2-quant/CpG-2-AllTilesPromvGenic.pdf")
ggplot(data=CpG.allhits, aes(x=feature_f, y=amount, fill=level)) + 
  geom_bar(stat="identity") + 
  facet_grid(~comp_f) +
  labs(title="Number of tiles in CpG context by comparison and feature", subtitle = "All hits",
       x="Feature type", y="Number of DM tiles") +
  scale_fill_manual(values=c("#2b83ba", "#d7191c"))
dev.off()

# Promoter vs. genic, genic tiles unique
CpG.feat <- data.frame (
  comp=c(rep("full",4), rep("A",4), rep("B",4), rep("C",4)),
  feature=c(rep("promoter",2), rep("genic",2)),
  level=c("hyper","hypo"),
  amount=c(length(which(feat.CpG$promoters$meth.diff>0)), length(which(feat.CpG$promoters$meth.diff<0)), 
           length(which(feat.CpG$genic$meth.diff>0)), length(which(feat.CpG$genic$meth.diff<0)),
           length(which(feat.aCpG$promoters$meth.diff>0)), length(which(feat.aCpG$promoters$meth.diff<0)), 
           length(which(feat.aCpG$genic$meth.diff>0)), length(which(feat.aCpG$genic$meth.diff<0)),
           length(which(feat.bCpG$promoters$meth.diff>0)), length(which(feat.bCpG$promoters$meth.diff<0)), 
           length(which(feat.bCpG$genic$meth.diff>0)), length(which(feat.bCpG$genic$meth.diff<0)),
           length(which(feat.cCpG$promoters$meth.diff>0)), length(which(feat.cCpG$promoters$meth.diff<0)), 
           length(which(feat.cCpG$genic$meth.diff>0)), length(which(feat.cCpG$genic$meth.diff<0)))
)
CpG.feat$comp_f=factor(CpG.feat$comp, levels=c("full", "A", "B", "C"))
CpG.feat$feature_f=factor(CpG.feat$feature, levels=c("promoter", "genic"))

pdf(file="4-Results/4.2-quant/CpG-3-PromvUniqueGenic.pdf")
ggplot(data=CpG.feat, aes(x=feature_f, y=amount, fill=level)) + 
  geom_bar(stat="identity") + 
  facet_grid(~comp_f) +
  labs(title="Number of DM tiles in CpG context by comparison and feature", subtitle = "Genic tiles unique",
       x="Feature type", y="Number of DM tiles") +
  scale_fill_manual(values=c("#2b83ba", "#d7191c"))
dev.off()

# Promoter vs. genic, non-overlapping tiles only
CpG.unique <- data.frame (
  comp=c(rep("full",4), rep("A",4), rep("B",4), rep("C",4)),
  feature=c(rep("promoter",2), rep("genic",2)),
  level=c("hyper","hypo"),
  amount=c(nrow(unique.CpG$hyperPromoters), nrow(unique.CpG$hypoPromoters), 
           nrow(unique.CpG$hyperGenic), nrow(unique.CpG$hypoGenic),
           nrow(unique.aCpG$hyperPromoters), nrow(unique.aCpG$hypoPromoters), 
           nrow(unique.aCpG$hyperGenic), nrow(unique.aCpG$hypoGenic),
           nrow(unique.bCpG$hyperPromoters), nrow(unique.bCpG$hypoPromoters), 
           nrow(unique.bCpG$hyperGenic), nrow(unique.bCpG$hypoGenic),
           nrow(unique.cCpG$hyperPromoters), nrow(unique.cCpG$hypoPromoters), 
           nrow(unique.cCpG$hyperGenic), nrow(unique.cCpG$hypoGenic))
)
CpG.unique$comp_f=factor(CpG.unique$comp, levels=c("full", "A", "B", "C"))
CpG.unique$feature_f=factor(CpG.unique$feature, levels=c("promoter", "genic"))

pdf(file="4-Results/4.2-quant/CpG-4-NonOverlappingTiles.pdf")
ggplot(data=CpG.unique, aes(x=feature_f, y=amount, fill=level)) + 
  geom_bar(stat="identity") + 
  facet_grid(~comp_f) +
  labs(title="Number of unique DM tiles in CpG context by comparison and feature", subtitle = "Only tiles mapping to a single feature",
       x="Feature type", y="Number of DM tiles") +
  scale_fill_manual(values=c("#2b83ba", "#d7191c"))
dev.off()

# Tile overlap between gene and promoter regions --------------------------
# CpG
# Between exons and introns
CpG.eiOverlap.table <- data.frame(
  comparison = c("full", "A","B","C"),
  exons=c(nrow(tile.CpG[[1]]), nrow(tile.aCpG[[1]]), nrow(tile.bCpG[[1]]), nrow(tile.cCpG[[1]])),
  introns=c(nrow(tile.CpG[[2]]), nrow(tile.aCpG[[2]]), nrow(tile.bCpG[[2]]), nrow(tile.cCpG[[2]])),
  overlap=c(nrow(intersect(tile.CpG[[1]][,1:2], tile.CpG[[2]][,1:2])),
            nrow(intersect(tile.aCpG[[1]][,1:2], tile.aCpG[[2]][,1:2])),
            nrow(intersect(tile.bCpG[[1]][,1:2], tile.bCpG[[2]][,1:2])),
            nrow(intersect(tile.cCpG[[1]][,1:2], tile.cCpG[[2]][,1:2])))
)

CpG.eiOverlap <- data.frame(
  comp = c(rep("full",3), rep("A",3),rep("B",3),rep("C",3)),
  feature = c("exons", "introns", "overlap"),
  amount = c(nrow(tile.CpG$exons), nrow(tile.CpG$introns), nrow(intersect(tile.CpG$exons[,1:2], tile.CpG$introns[,1:2])),
             nrow(tile.aCpG$exons), nrow(tile.aCpG$introns), nrow(intersect(tile.aCpG$exons[,1:2], tile.aCpG$introns[,1:2])),
             nrow(tile.bCpG$exons), nrow(tile.bCpG$introns), nrow(intersect(tile.bCpG$exons[,1:2], tile.bCpG$introns[,1:2])),
             nrow(tile.cCpG$exons), nrow(tile.cCpG$introns), nrow(intersect(tile.cCpG$exons[,1:2], tile.cCpG$introns[,1:2])))
)
CpG.eiOverlap$comp_f=factor(CpG.eiOverlap$comp, levels=c("full", "A", "B", "C"))
CpG.eiOverlap$feature_f=factor(CpG.eiOverlap$feature, levels=c("exons", "introns","overlap"))
pdf(file="4-Results/4.2-quant/CpG-5-ExonIntronOverlap.pdf")
ggplot(data=CpG.eiOverlap, aes(x=feature_f, y=amount, fill=feature_f)) + 
  geom_bar(stat="identity",position="dodge") + 
  facet_grid(~comp_f) +
  labs(title="DM CpG tiles overlapping both introns and exons",
       x="Feature type", y="Number of DM tiles") 
dev.off()

#Promoters and genic
CpG.pgOverlap.table <- data.frame(
  comparison = c("full", "A","B","C"),
  promoters=c(nrow(feat.CpG[[1]]), nrow(feat.aCpG[[1]]), nrow(feat.bCpG[[1]]), nrow(feat.cCpG[[1]])),
  genic=c(nrow(feat.CpG[[2]]), nrow(feat.aCpG[[2]]), nrow(feat.bCpG[[2]]), nrow(feat.cCpG[[2]])),
  overlap=c(nrow(intersect(feat.CpG[[1]][,1:2], feat.CpG[[2]][,1:2])),
            nrow(intersect(feat.aCpG[[1]][,1:2], feat.aCpG[[2]][,1:2])),
            nrow(intersect(feat.bCpG[[1]][,1:2], feat.bCpG[[2]][,1:2])),
            nrow(intersect(feat.cCpG[[1]][,1:2], feat.cCpG[[2]][,1:2])))
)
CpG.pgOverlap <- data.frame(
  comp = c(rep("full",3), rep("A",3),rep("B",3),rep("C",3)),
  feature = c("promoters", "genic", "overlap"),
  amount = c(nrow(feat.CpG$promoters), nrow(feat.CpG$genic), nrow(intersect(feat.CpG$promoters[,1:2], feat.CpG$genic[,1:2])),
            nrow(feat.aCpG$promoters), nrow(feat.aCpG$genic), nrow(intersect(feat.aCpG$promoters[,1:2], feat.aCpG$genic[,1:2])),
            nrow(feat.bCpG$promoters), nrow(feat.bCpG$genic), nrow(intersect(feat.bCpG$promoters[,1:2], feat.bCpG$genic[,1:2])),
            nrow(feat.cCpG$promoters), nrow(feat.cCpG$genic), nrow(intersect(feat.cCpG$promoters[,1:2], feat.cCpG$genic[,1:2])))
  )
CpG.pgOverlap$comp_f=factor(CpG.pgOverlap$comp, levels=c("full", "A", "B", "C"))
CpG.pgOverlap$feature_f=factor(CpG.pgOverlap$feature, levels=c("promoters", "genic","overlap"))
pdf(file="4-Results/4.2-quant/CpG-6-PromoterGenicOverlap.pdf")
ggplot(data=CpG.pgOverlap, aes(x=feature_f, y=amount, fill=feature_f)) + 
  geom_bar(stat="identity",position="dodge") + 
  facet_grid(~comp_f) +
  labs(title="DM CpG tiles overlapping both genic and promoter regions",
       x="Feature type", y="Number of DM tiles") 
dev.off()

# CHG ---------------------------------------------------------------------

# Quantifying raw numbers -------------------------------------------------
# Feature counts by genotype, all hits
CHG.tile <- data.frame (
  comp=c(rep("full",6), rep("A",6), rep("B",6), rep("C",6)),
  feature=c(rep("promoter",2), rep("exon",2), rep("intron",2)),
  level=c("hyper","hypo"),
  amount=c(length(which(tile.CHG$promoters$meth.diff>0)), length(which(tile.CHG$promoters$meth.diff<0)), 
           length(which(tile.CHG$exons$meth.diff>0)), length(which(tile.CHG$exons$meth.diff<0)),
           length(which(tile.CHG$introns$meth.diff>0)), length(which(tile.CHG$introns$meth.diff<0)),
           length(which(tile.aCHG$promoters$meth.diff>0)), length(which(tile.aCHG$promoters$meth.diff<0)), 
           length(which(tile.aCHG$exons$meth.diff>0)), length(which(tile.aCHG$exons$meth.diff<0)),
           length(which(tile.aCHG$introns$meth.diff>0)), length(which(tile.aCHG$introns$meth.diff<0)),
           length(which(tile.bCHG$promoters$meth.diff>0)), length(which(tile.bCHG$promoters$meth.diff<0)), 
           length(which(tile.bCHG$exons$meth.diff>0)), length(which(tile.bCHG$exons$meth.diff<0)),
           length(which(tile.bCHG$introns$meth.diff>0)), length(which(tile.bCHG$introns$meth.diff<0)),
           length(which(tile.cCHG$promoters$meth.diff>0)), length(which(tile.cCHG$promoters$meth.diff<0)), 
           length(which(tile.cCHG$exons$meth.diff>0)), length(which(tile.cCHG$exons$meth.diff<0)),
           length(which(tile.cCHG$introns$meth.diff>0)), length(which(tile.cCHG$introns$meth.diff<0)))
)
CHG.tile$comp_f=factor(CHG.tile$comp, levels=c("full", "A", "B", "C"))
CHG.tile$feature_f=factor(CHG.tile$feature, levels=c("promoter", "exon", "intron"))

pdf(file="4-Results/4.2-quant/CHG-1-AllTilesByFeat.pdf")
ggplot(data=CHG.tile, aes(x=feature_f, y=amount, fill=level)) + 
  geom_bar(stat="identity") + 
  facet_grid(~comp_f) +
  labs(title="Number of DM tiles in CHG context by comparison and feature", x="Feature type", y="Number of DM tiles") +
  scale_fill_manual(values=c("#2b83ba", "#d7191c"))
dev.off()

# Promoter vs. genic, all hits
CHG.allhits <- data.frame (
  comp=c(rep("full",4), rep("A",4), rep("B",4), rep("C",4)),
  feature=c(rep("promoter",2), rep("genic",2)),
  level=c("hyper","hypo"),
  amount=c(nrow(split3.CHG$hyperPromoters), nrow(split3.CHG$hypoPromoters), 
           nrow(split3.CHG$hyperGenic), nrow(split3.CHG$hypoGenic),
           nrow(split3.aCHG$hyperPromoters), nrow(split3.aCHG$hypoPromoters), 
           nrow(split3.aCHG$hyperGenic), nrow(split3.aCHG$hypoGenic),
           nrow(split3.bCHG$hyperPromoters), nrow(split3.bCHG$hypoPromoters), 
           nrow(split3.bCHG$hyperGenic), nrow(split3.bCHG$hypoGenic),
           nrow(split3.cCHG$hyperPromoters), nrow(split3.cCHG$hypoPromoters), 
           nrow(split3.cCHG$hyperGenic), nrow(split3.cCHG$hypoGenic))
)
CHG.allhits$comp_f=factor(CHG.allhits$comp, levels=c("full", "A", "B", "C"))
CHG.allhits$feature_f=factor(CHG.allhits$feature, levels=c("promoter", "genic"))

pdf(file="4-Results/4.2-quant/CHG-2-AllTilesPromvGenic.pdf")
ggplot(data=CHG.allhits, aes(x=feature_f, y=amount, fill=level)) + 
  geom_bar(stat="identity") + 
  facet_grid(~comp_f) +
  labs(title="Number of tiles in CHG context by comparison and feature", subtitle = "All hits",
       x="Feature type", y="Number of DM tiles") +
  scale_fill_manual(values=c("#2b83ba", "#d7191c"))
dev.off()

# Promoter vs. genic, genic tiles unique
CHG.feat <- data.frame (
  comp=c(rep("full",4), rep("A",4), rep("B",4), rep("C",4)),
  feature=c(rep("promoter",2), rep("genic",2)),
  level=c("hyper","hypo"),
  amount=c(length(which(feat.CHG$promoters$meth.diff>0)), length(which(feat.CHG$promoters$meth.diff<0)), 
           length(which(feat.CHG$genic$meth.diff>0)), length(which(feat.CHG$genic$meth.diff<0)),
           length(which(feat.aCHG$promoters$meth.diff>0)), length(which(feat.aCHG$promoters$meth.diff<0)), 
           length(which(feat.aCHG$genic$meth.diff>0)), length(which(feat.aCHG$genic$meth.diff<0)),
           length(which(feat.bCHG$promoters$meth.diff>0)), length(which(feat.bCHG$promoters$meth.diff<0)), 
           length(which(feat.bCHG$genic$meth.diff>0)), length(which(feat.bCHG$genic$meth.diff<0)),
           length(which(feat.cCHG$promoters$meth.diff>0)), length(which(feat.cCHG$promoters$meth.diff<0)), 
           length(which(feat.cCHG$genic$meth.diff>0)), length(which(feat.cCHG$genic$meth.diff<0)))
)
CHG.feat$comp_f=factor(CHG.feat$comp, levels=c("full", "A", "B", "C"))
CHG.feat$feature_f=factor(CHG.feat$feature, levels=c("promoter", "genic"))

pdf(file="4-Results/4.2-quant/CHG-3-PromvUniqueGenic.pdf")
ggplot(data=CHG.feat, aes(x=feature_f, y=amount, fill=level)) + 
  geom_bar(stat="identity") + 
  facet_grid(~comp_f) +
  labs(title="Number of DM tiles in CHG context by comparison and feature", subtitle = "Genic tiles unique",
       x="Feature type", y="Number of DM tiles") +
  scale_fill_manual(values=c("#2b83ba", "#d7191c"))
dev.off()

# Promoter vs. genic, non-overlapping tiles only
CHG.unique <- data.frame (
  comp=c(rep("full",4), rep("A",4), rep("B",4), rep("C",4)),
  feature=c(rep("promoter",2), rep("genic",2)),
  level=c("hyper","hypo"),
  amount=c(nrow(unique.CHG$hyperPromoters), nrow(unique.CHG$hypoPromoters), 
           nrow(unique.CHG$hyperGenic), nrow(unique.CHG$hypoGenic),
           nrow(unique.aCHG$hyperPromoters), nrow(unique.aCHG$hypoPromoters), 
           nrow(unique.aCHG$hyperGenic), nrow(unique.aCHG$hypoGenic),
           nrow(unique.bCHG$hyperPromoters), nrow(unique.bCHG$hypoPromoters), 
           nrow(unique.bCHG$hyperGenic), nrow(unique.bCHG$hypoGenic),
           nrow(unique.cCHG$hyperPromoters), nrow(unique.cCHG$hypoPromoters), 
           nrow(unique.cCHG$hyperGenic), nrow(unique.cCHG$hypoGenic))
)
CHG.unique$comp_f=factor(CHG.unique$comp, levels=c("full", "A", "B", "C"))
CHG.unique$feature_f=factor(CHG.unique$feature, levels=c("promoter", "genic"))

pdf(file="4-Results/4.2-quant/CHG-4-NonOverlappingTiles.pdf")
ggplot(data=CHG.unique, aes(x=feature_f, y=amount, fill=level)) + 
  geom_bar(stat="identity") + 
  facet_grid(~comp_f) +
  labs(title="Number of unique DM tiles in CHG context by comparison and feature", subtitle = "Only tiles mapping to a single feature",
       x="Feature type", y="Number of DM tiles") +
  scale_fill_manual(values=c("#2b83ba", "#d7191c"))
dev.off()

# Tile overlap between gene and promoter regions --------------------------
# CHG
# Between exons and introns
CHG.eiOverlap.table <- data.frame(
  comparison = c("full", "A","B","C"),
  exons=c(nrow(tile.CHG[[1]]), nrow(tile.aCHG[[1]]), nrow(tile.bCHG[[1]]), nrow(tile.cCHG[[1]])),
  introns=c(nrow(tile.CHG[[2]]), nrow(tile.aCHG[[2]]), nrow(tile.bCHG[[2]]), nrow(tile.cCHG[[2]])),
  overlap=c(nrow(intersect(tile.CHG[[1]][,1:2], tile.CHG[[2]][,1:2])),
            nrow(intersect(tile.aCHG[[1]][,1:2], tile.aCHG[[2]][,1:2])),
            nrow(intersect(tile.bCHG[[1]][,1:2], tile.bCHG[[2]][,1:2])),
            nrow(intersect(tile.cCHG[[1]][,1:2], tile.cCHG[[2]][,1:2])))
)

CHG.eiOverlap <- data.frame(
  comp = c(rep("full",3), rep("A",3),rep("B",3),rep("C",3)),
  feature = c("exons", "introns", "overlap"),
  amount = c(nrow(tile.CHG$exons), nrow(tile.CHG$introns), nrow(intersect(tile.CHG$exons[,1:2], tile.CHG$introns[,1:2])),
             nrow(tile.aCHG$exons), nrow(tile.aCHG$introns), nrow(intersect(tile.aCHG$exons[,1:2], tile.aCHG$introns[,1:2])),
             nrow(tile.bCHG$exons), nrow(tile.bCHG$introns), nrow(intersect(tile.bCHG$exons[,1:2], tile.bCHG$introns[,1:2])),
             nrow(tile.cCHG$exons), nrow(tile.cCHG$introns), nrow(intersect(tile.cCHG$exons[,1:2], tile.cCHG$introns[,1:2])))
)
CHG.eiOverlap$comp_f=factor(CHG.eiOverlap$comp, levels=c("full", "A", "B", "C"))
CHG.eiOverlap$feature_f=factor(CHG.eiOverlap$feature, levels=c("exons", "introns","overlap"))
pdf(file="4-Results/4.2-quant/CHG-5-ExonIntronOverlap.pdf")
ggplot(data=CHG.eiOverlap, aes(x=feature_f, y=amount, fill=feature_f)) + 
  geom_bar(stat="identity",position="dodge") + 
  facet_grid(~comp_f) +
  labs(title="DM CHG tiles overlapping both introns and exons",
       x="Feature type", y="Number of DM tiles") 
dev.off()

#Promoters and genic
CHG.pgOverlap.table <- data.frame(
  comparison = c("full", "A","B","C"),
  promoters=c(nrow(feat.CHG[[1]]), nrow(feat.aCHG[[1]]), nrow(feat.bCHG[[1]]), nrow(feat.cCHG[[1]])),
  genic=c(nrow(feat.CHG[[2]]), nrow(feat.aCHG[[2]]), nrow(feat.bCHG[[2]]), nrow(feat.cCHG[[2]])),
  overlap=c(nrow(intersect(feat.CHG[[1]][,1:2], feat.CHG[[2]][,1:2])),
            nrow(intersect(feat.aCHG[[1]][,1:2], feat.aCHG[[2]][,1:2])),
            nrow(intersect(feat.bCHG[[1]][,1:2], feat.bCHG[[2]][,1:2])),
            nrow(intersect(feat.cCHG[[1]][,1:2], feat.cCHG[[2]][,1:2])))
)
CHG.pgOverlap <- data.frame(
  comp = c(rep("full",3), rep("A",3),rep("B",3),rep("C",3)),
  feature = c("promoters", "genic", "overlap"),
  amount = c(nrow(feat.CHG$promoters), nrow(feat.CHG$genic), nrow(intersect(feat.CHG$promoters[,1:2], feat.CHG$genic[,1:2])),
             nrow(feat.aCHG$promoters), nrow(feat.aCHG$genic), nrow(intersect(feat.aCHG$promoters[,1:2], feat.aCHG$genic[,1:2])),
             nrow(feat.bCHG$promoters), nrow(feat.bCHG$genic), nrow(intersect(feat.bCHG$promoters[,1:2], feat.bCHG$genic[,1:2])),
             nrow(feat.cCHG$promoters), nrow(feat.cCHG$genic), nrow(intersect(feat.cCHG$promoters[,1:2], feat.cCHG$genic[,1:2])))
)
CHG.pgOverlap$comp_f=factor(CHG.pgOverlap$comp, levels=c("full", "A", "B", "C"))
CHG.pgOverlap$feature_f=factor(CHG.pgOverlap$feature, levels=c("promoters", "genic","overlap"))
pdf(file="4-Results/4.2-quant/CHG-6-PromoterGenicOverlap.pdf")
ggplot(data=CHG.pgOverlap, aes(x=feature_f, y=amount, fill=feature_f)) + 
  geom_bar(stat="identity",position="dodge") + 
  facet_grid(~comp_f) +
  labs(title="DM CHG tiles overlapping both genic and promoter regions",
       x="Feature type", y="Number of DM tiles") 
dev.off()
  
# CHH ---------------------------------------------------------------------

# Quantifying raw numbers -------------------------------------------------
# Feature counts by genotype, all hits
CHH.tile <- data.frame (
  comp=c(rep("full",6), rep("A",6), rep("B",6), rep("C",6)),
  feature=c(rep("promoter",2), rep("exon",2), rep("intron",2)),
  level=c("hyper","hypo"),
  amount=c(length(which(tile.CHH$promoters$meth.diff>0)), length(which(tile.CHH$promoters$meth.diff<0)), 
           length(which(tile.CHH$exons$meth.diff>0)), length(which(tile.CHH$exons$meth.diff<0)),
           length(which(tile.CHH$introns$meth.diff>0)), length(which(tile.CHH$introns$meth.diff<0)),
           length(which(tile.aCHH$promoters$meth.diff>0)), length(which(tile.aCHH$promoters$meth.diff<0)), 
           length(which(tile.aCHH$exons$meth.diff>0)), length(which(tile.aCHH$exons$meth.diff<0)),
           length(which(tile.aCHH$introns$meth.diff>0)), length(which(tile.aCHH$introns$meth.diff<0)),
           length(which(tile.bCHH$promoters$meth.diff>0)), length(which(tile.bCHH$promoters$meth.diff<0)), 
           length(which(tile.bCHH$exons$meth.diff>0)), length(which(tile.bCHH$exons$meth.diff<0)),
           length(which(tile.bCHH$introns$meth.diff>0)), length(which(tile.bCHH$introns$meth.diff<0)),
           length(which(tile.cCHH$promoters$meth.diff>0)), length(which(tile.cCHH$promoters$meth.diff<0)), 
           length(which(tile.cCHH$exons$meth.diff>0)), length(which(tile.cCHH$exons$meth.diff<0)),
           length(which(tile.cCHH$introns$meth.diff>0)), length(which(tile.cCHH$introns$meth.diff<0)))
)
CHH.tile$comp_f=factor(CHH.tile$comp, levels=c("full", "A", "B", "C"))
CHH.tile$feature_f=factor(CHH.tile$feature, levels=c("promoter", "exon", "intron"))

pdf(file="4-Results/4.2-quant/CHH-1-AllTilesByFeat.pdf")
ggplot(data=CHH.tile, aes(x=feature_f, y=amount, fill=level)) + 
  geom_bar(stat="identity") + 
  facet_grid(~comp_f) +
  labs(title="Number of DM tiles in CHH context by comparison and feature", x="Feature type", y="Number of DM tiles") +
  scale_fill_manual(values=c("#2b83ba", "#d7191c"))
dev.off()

# Promoter vs. genic, all hits
CHH.allhits <- data.frame (
  comp=c(rep("full",4), rep("A",4), rep("B",4), rep("C",4)),
  feature=c(rep("promoter",2), rep("genic",2)),
  level=c("hyper","hypo"),
  amount=c(nrow(split3.CHH$hyperPromoters), nrow(split3.CHH$hypoPromoters), 
           nrow(split3.CHH$hyperGenic), nrow(split3.CHH$hypoGenic),
           nrow(split3.aCHH$hyperPromoters), nrow(split3.aCHH$hypoPromoters), 
           nrow(split3.aCHH$hyperGenic), nrow(split3.aCHH$hypoGenic),
           nrow(split3.bCHH$hyperPromoters), nrow(split3.bCHH$hypoPromoters), 
           nrow(split3.bCHH$hyperGenic), nrow(split3.bCHH$hypoGenic),
           nrow(split3.cCHH$hyperPromoters), nrow(split3.cCHH$hypoPromoters), 
           nrow(split3.cCHH$hyperGenic), nrow(split3.cCHH$hypoGenic))
)
CHH.allhits$comp_f=factor(CHH.allhits$comp, levels=c("full", "A", "B", "C"))
CHH.allhits$feature_f=factor(CHH.allhits$feature, levels=c("promoter", "genic"))

pdf(file="4-Results/4.2-quant/CHH-2-AllTilesPromvGenic.pdf")
ggplot(data=CHH.allhits, aes(x=feature_f, y=amount, fill=level)) + 
  geom_bar(stat="identity") + 
  facet_grid(~comp_f) +
  labs(title="Number of tiles in CHH context by comparison and feature", subtitle = "All hits",
       x="Feature type", y="Number of DM tiles") +
  scale_fill_manual(values=c("#2b83ba", "#d7191c"))
dev.off()

# Promoter vs. genic, genic tiles unique
CHH.feat <- data.frame (
  comp=c(rep("full",4), rep("A",4), rep("B",4), rep("C",4)),
  feature=c(rep("promoter",2), rep("genic",2)),
  level=c("hyper","hypo"),
  amount=c(length(which(feat.CHH$promoters$meth.diff>0)), length(which(feat.CHH$promoters$meth.diff<0)), 
           length(which(feat.CHH$genic$meth.diff>0)), length(which(feat.CHH$genic$meth.diff<0)),
           length(which(feat.aCHH$promoters$meth.diff>0)), length(which(feat.aCHH$promoters$meth.diff<0)), 
           length(which(feat.aCHH$genic$meth.diff>0)), length(which(feat.aCHH$genic$meth.diff<0)),
           length(which(feat.bCHH$promoters$meth.diff>0)), length(which(feat.bCHH$promoters$meth.diff<0)), 
           length(which(feat.bCHH$genic$meth.diff>0)), length(which(feat.bCHH$genic$meth.diff<0)),
           length(which(feat.cCHH$promoters$meth.diff>0)), length(which(feat.cCHH$promoters$meth.diff<0)), 
           length(which(feat.cCHH$genic$meth.diff>0)), length(which(feat.cCHH$genic$meth.diff<0)))
)
CHH.feat$comp_f=factor(CHH.feat$comp, levels=c("full", "A", "B", "C"))
CHH.feat$feature_f=factor(CHH.feat$feature, levels=c("promoter", "genic"))

pdf(file="4-Results/4.2-quant/CHH-3-PromvUniqueGenic.pdf")
ggplot(data=CHH.feat, aes(x=feature_f, y=amount, fill=level)) + 
  geom_bar(stat="identity") + 
  facet_grid(~comp_f) +
  labs(title="Number of DM tiles in CHH context by comparison and feature", subtitle = "Genic tiles unique",
       x="Feature type", y="Number of DM tiles") +
  scale_fill_manual(values=c("#2b83ba", "#d7191c"))
dev.off()

# Promoter vs. genic, non-overlapping tiles only
CHH.unique <- data.frame (
  comp=c(rep("full",4), rep("A",4), rep("B",4), rep("C",4)),
  feature=c(rep("promoter",2), rep("genic",2)),
  level=c("hyper","hypo"),
  amount=c(nrow(unique.CHH$hyperPromoters), nrow(unique.CHH$hypoPromoters), 
           nrow(unique.CHH$hyperGenic), nrow(unique.CHH$hypoGenic),
           nrow(unique.aCHH$hyperPromoters), nrow(unique.aCHH$hypoPromoters), 
           nrow(unique.aCHH$hyperGenic), nrow(unique.aCHH$hypoGenic),
           nrow(unique.bCHH$hyperPromoters), nrow(unique.bCHH$hypoPromoters), 
           nrow(unique.bCHH$hyperGenic), nrow(unique.bCHH$hypoGenic),
           nrow(unique.cCHH$hyperPromoters), nrow(unique.cCHH$hypoPromoters), 
           nrow(unique.cCHH$hyperGenic), nrow(unique.cCHH$hypoGenic))
)
CHH.unique$comp_f=factor(CHH.unique$comp, levels=c("full", "A", "B", "C"))
CHH.unique$feature_f=factor(CHH.unique$feature, levels=c("promoter", "genic"))

pdf(file="4-Results/4.2-quant/CHH-4-NonOverlappingTiles.pdf")
ggplot(data=CHH.unique, aes(x=feature_f, y=amount, fill=level)) + 
  geom_bar(stat="identity") + 
  facet_grid(~comp_f) +
  labs(title="Number of unique DM tiles in CHH context by comparison and feature", subtitle = "Only tiles mapping to a single feature",
       x="Feature type", y="Number of DM tiles") +
  scale_fill_manual(values=c("#2b83ba", "#d7191c"))
dev.off()

# Tile overlap between gene and promoter regions --------------------------
# CHH
# Between exons and introns
CHH.eiOverlap.table <- data.frame(
  comparison = c("full", "A","B","C"),
  exons=c(nrow(tile.CHH[[1]]), nrow(tile.aCHH[[1]]), nrow(tile.bCHH[[1]]), nrow(tile.cCHH[[1]])),
  introns=c(nrow(tile.CHH[[2]]), nrow(tile.aCHH[[2]]), nrow(tile.bCHH[[2]]), nrow(tile.cCHH[[2]])),
  overlap=c(nrow(intersect(tile.CHH[[1]][,1:2], tile.CHH[[2]][,1:2])),
            nrow(intersect(tile.aCHH[[1]][,1:2], tile.aCHH[[2]][,1:2])),
            nrow(intersect(tile.bCHH[[1]][,1:2], tile.bCHH[[2]][,1:2])),
            nrow(intersect(tile.cCHH[[1]][,1:2], tile.cCHH[[2]][,1:2])))
)

CHH.eiOverlap <- data.frame(
  comp = c(rep("full",3), rep("A",3),rep("B",3),rep("C",3)),
  feature = c("exons", "introns", "overlap"),
  amount = c(nrow(tile.CHH$exons), nrow(tile.CHH$introns), nrow(intersect(tile.CHH$exons[,1:2], tile.CHH$introns[,1:2])),
             nrow(tile.aCHH$exons), nrow(tile.aCHH$introns), nrow(intersect(tile.aCHH$exons[,1:2], tile.aCHH$introns[,1:2])),
             nrow(tile.bCHH$exons), nrow(tile.bCHH$introns), nrow(intersect(tile.bCHH$exons[,1:2], tile.bCHH$introns[,1:2])),
             nrow(tile.cCHH$exons), nrow(tile.cCHH$introns), nrow(intersect(tile.cCHH$exons[,1:2], tile.cCHH$introns[,1:2])))
)
CHH.eiOverlap$comp_f=factor(CHH.eiOverlap$comp, levels=c("full", "A", "B", "C"))
CHH.eiOverlap$feature_f=factor(CHH.eiOverlap$feature, levels=c("exons", "introns","overlap"))
pdf(file="4-Results/4.2-quant/CHH-5-ExonIntronOverlap.pdf")
ggplot(data=CHH.eiOverlap, aes(x=feature_f, y=amount, fill=feature_f)) + 
  geom_bar(stat="identity",position="dodge") + 
  facet_grid(~comp_f) +
  labs(title="DM CHH tiles overlapping both introns and exons",
       x="Feature type", y="Number of DM tiles") 
dev.off()

#Promoters and genic
CHH.pgOverlap.table <- data.frame(
  comparison = c("full", "A","B","C"),
  promoters=c(nrow(feat.CHH[[1]]), nrow(feat.aCHH[[1]]), nrow(feat.bCHH[[1]]), nrow(feat.cCHH[[1]])),
  genic=c(nrow(feat.CHH[[2]]), nrow(feat.aCHH[[2]]), nrow(feat.bCHH[[2]]), nrow(feat.cCHH[[2]])),
  overlap=c(nrow(intersect(feat.CHH[[1]][,1:2], feat.CHH[[2]][,1:2])),
            nrow(intersect(feat.aCHH[[1]][,1:2], feat.aCHH[[2]][,1:2])),
            nrow(intersect(feat.bCHH[[1]][,1:2], feat.bCHH[[2]][,1:2])),
            nrow(intersect(feat.cCHH[[1]][,1:2], feat.cCHH[[2]][,1:2])))
)
CHH.pgOverlap <- data.frame(
  comp = c(rep("full",3), rep("A",3),rep("B",3),rep("C",3)),
  feature = c("promoters", "genic", "overlap"),
  amount = c(nrow(feat.CHH$promoters), nrow(feat.CHH$genic), nrow(intersect(feat.CHH$promoters[,1:2], feat.CHH$genic[,1:2])),
             nrow(feat.aCHH$promoters), nrow(feat.aCHH$genic), nrow(intersect(feat.aCHH$promoters[,1:2], feat.aCHH$genic[,1:2])),
             nrow(feat.bCHH$promoters), nrow(feat.bCHH$genic), nrow(intersect(feat.bCHH$promoters[,1:2], feat.bCHH$genic[,1:2])),
             nrow(feat.cCHH$promoters), nrow(feat.cCHH$genic), nrow(intersect(feat.cCHH$promoters[,1:2], feat.cCHH$genic[,1:2])))
)
CHH.pgOverlap$comp_f=factor(CHH.pgOverlap$comp, levels=c("full", "A", "B", "C"))
CHH.pgOverlap$feature_f=factor(CHH.pgOverlap$feature, levels=c("promoters", "genic","overlap"))
pdf(file="4-Results/4.2-quant/CHH-6-PromoterGenicOverlap.pdf")
ggplot(data=CHH.pgOverlap, aes(x=feature_f, y=amount, fill=feature_f)) + 
  geom_bar(stat="identity",position="dodge") + 
  facet_grid(~comp_f) +
  labs(title="DM CHH tiles overlapping both genic and promoter regions",
       x="Feature type", y="Number of DM tiles") 
dev.off()