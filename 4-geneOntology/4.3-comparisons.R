###################################################################
### VENN DIAGRAMS OF ANNOTATED GENES
### Local R script 
###################################################################

# set up environment
setwd("D:/Nia/OneDrive - University of Guelph/Working/DMRanalysis")
source("D:/Nia/OneDrive - University of Guelph/Working/thesisScripts/4-geneOntology/4.1-dataPrep.R")
library(UpSetR)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(cowplot)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# CpG ---------------------------------------------------------------------

# Hypo vs Hyper -----------------------------------------------------------
pdf(file="4-Results/4.3-comparisons/hyperHypo/CpG-hyperHypo-upset.pdf")
upset(fromList(list(
  fullHyper=hyp.CpG$hyper,
  AHyper=hyp.aCpG$hyper,
  BHyper=hyp.bCpG$hyper,
  CHyper=hyp.cCpG$hyper,
  fullHypo=hyp.CpG$hypo,
  AHypo=hyp.aCpG$hypo,
  BHypo=hyp.bCpG$hypo,
  CHypo=hyp.cCpG$hypo)), 
  mainbar.y.label = "Gene Intersections", sets.x.label = "DM genes per genotype", order.by="freq", number.angles=30,  text.scale = c(1.3, 1.3, 1, 1, 1, 0.75),
  sets=c("CHypo", "BHypo", "AHypo", "fullHypo","CHyper", "BHyper", "AHyper", "fullHyper"), keep.order=TRUE)
grid.text("Hyper- vs. hypomethylated genes across comparisons, CpG context",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

grid.newpage()
venn.diagram(x=list(unique(hyp.CpG$hyper),unique(hyp.cCpG$hyper),unique(hyp.aCpG$hyper),unique(hyp.bCpG$hyper)),
             category.names = c("full","C","A","B"), main="Hypermethylated genes across comparisons",
             filename="4-Results/4.3-comparisons/hyperHypo/CpG-hypermethylated.tiff",imagetype="tiff",
             fill=brewer.pal(4,"Set1"))
grid.newpage()
venn.diagram(x=list(unique(hyp.CpG$hypo),unique(hyp.cCpG$hypo),unique(hyp.aCpG$hypo),unique(hyp.bCpG$hypo)),
             category.names = c("full","C","A","B"), main="Hypomethylated genes across comparisons",
             filename="4-Results/4.3-comparisons/hyperHypo/CpG-hypomethylated.tiff",imagetype="tiff",          
             fill=brewer.pal(4,"Set1"))

grid.newpage()
venn.diagram(x=list(unique(hyp.CpG$hyper),unique(hyp.CpG$hypo)),
             category.names = c("hyper","hypo"), main="Hyper- vs hypomethylated genes in the full model",
             filename="4-Results/4.3-comparisons/hyperHypo/CpG-full-levelcomp.tiff",imagetype="tiff",
             fill=c("red", "blue"))
grid.newpage()
venn.diagram(x=list(unique(hyp.aCpG$hyper),unique(hyp.aCpG$hypo)),
             category.names = c("hyper","hypo"), main="Hyper- vs hypomethylated genes in A background",
             filename="4-Results/4.3-comparisons/hyperHypo/CpG-A-levelcomp.tiff",imagetype="tiff",
             fill=c("red", "blue"))
grid.newpage()
venn.diagram(x=list(unique(hyp.bCpG$hyper),unique(hyp.bCpG$hypo)),
             category.names = c("hyper","hypo"), main="Hyper- vs hypomethylated genes in B background",
             filename="4-Results/4.3-comparisons/hyperHypo/CpG-B-levelcomp.tiff",imagetype="tiff",
             fill=c("red", "blue"))
grid.newpage()
venn.diagram(x=list(unique(hyp.cCpG$hyper),unique(hyp.cCpG$hypo)),
             category.names = c("hyper","hypo"), main="Hyper- vs hypomethylated genes in C background",
             filename="4-Results/4.3-comparisons/hyperHypo/CpG-C-levelcomp.tiff",imagetype="tiff",
             fill=c("red", "blue"))


# Comparing genotype-specific DMRs ----------------------------------------
# Hypomethylated promoters
grid.newpage()
venn.diagram(x=list(gene.CpG$hypoPromoters, gene.cCpG$hypoPromoters,gene.aCpG$hypoPromoters,gene.bCpG$hypoPromoters),
             category.names = c("full","C","A","B"), main="Hypomethylated CpG promoters across comparisons",
             filename="4-Results/4.3-comparisons/genotypeSpecific/CpG-hypoProm-venn.tiff",imagetype="tiff",
             fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/genotypeSpecific/CpG-hypoProm-upset.pdf")
upset(fromList(list(full=gene.CpG$hypoPromoters, 
                    A=gene.aCpG$hypoPromoters, 
                    B=gene.bCpG$hypoPromoters, 
                    C=gene.cCpG$hypoPromoters)),
      order.by="freq", sets=c("C","B","A", "full"), keep.order=TRUE)
grid.text("Hypomethylated CpG promoters across comparisons",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# Hypermethylated promoters
grid.newpage()
venn.diagram(x=list(gene.CpG$hyperPromoters, gene.cCpG$hyperPromoters,gene.aCpG$hyperPromoters,gene.bCpG$hyperPromoters),
                       category.names = c("full","C","A","B"), main="Hypermethylated CpG promoters across comparisons",
                       filename="4-Results/4.3-comparisons/genotypeSpecific/CpG-hyperProm-venn.tiff",imagetype="tiff", fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/genotypeSpecific/CpG-hyperProm-upset.pdf")
upset(fromList(list(full=gene.CpG$hyperPromoters, 
                    A=gene.aCpG$hyperPromoters, 
                    B=gene.bCpG$hyperPromoters, 
                    C=gene.cCpG$hyperPromoters)),
      order.by="freq", sets=c("C","B","A", "full"), keep.order=TRUE)
grid.text("Hypermethylated CpG promoters across comparisons",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# Hypomethylated genic
grid.newpage()
venn.diagram(x=list(gene.CpG$hypoGenic, gene.cCpG$hypoGenic,gene.aCpG$hypoGenic,gene.bCpG$hypoGenic),
                       category.names = c("full","C","A","B"), main="Hypomethylated CpG genic regions across comparisons",
                       filename="4-Results/4.3-comparisons/genotypeSpecific/CpG-hypoGenic-venn.tiff",imagetype="tiff", fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/genotypeSpecific/CpG-hypoGenic-upset.pdf")
upset(fromList(list(full=gene.CpG$hypoGenic, 
                    A=gene.aCpG$hypoGenic, 
                    B=gene.bCpG$hypoGenic, 
                    C=gene.cCpG$hypoGenic)),
      order.by="freq", sets=c("C","B","A", "full"), keep.order=TRUE)
grid.text("Hypomethylated CpG genic regions across comparisons",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# Hypermethylated genic
grid.newpage()
venn.diagram(x=list(gene.CpG$hyperGenic, gene.cCpG$hyperGenic,gene.aCpG$hyperGenic,gene.bCpG$hyperGenic),
                       category.names = c("full","C","A","B"), main="Hypermethylated CpG genic regions across comparisons",
                       filename="4-Results/4.3-comparisons/genotypeSpecific/CpG-hyperGenic-venn.tiff",imagetype="tiff",fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/genotypeSpecific/CpG-hyperGenic-upset.pdf")
upset(fromList(list(full=gene.CpG$hyperGenic, 
                    A=gene.aCpG$hyperGenic, 
                    B=gene.bCpG$hyperGenic, 
                    C=gene.cCpG$hyperGenic)),
      order.by="freq", sets=c("C","B","A", "full"), keep.order=TRUE)
grid.text("Hypermethylated CpG genic regions across comparisons",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# Promoter vs. genic ------------------------------------------------------

# full
grid.newpage()
venn.diagram(x=list(hypoPromoters=gene.CpG$hypoPromoters, 
                              hyperPromoters=gene.CpG$hyperPromoters, 
                              hypoGenic=gene.CpG$hypoGenic, 
                              hyperGenic=gene.CpG$hyperGenic),
                       category.names = c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), 
                       main="Full model comparison of CpG hypo vs. hyper promoter vs. genic regions",
                       filename="4-Results/4.3-comparisons/promoterGenic/CpG-full-venn.tiff",imagetype="tiff",fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/promoterGenic/CpG-full-upset.pdf")
upset(fromList(list(hypoPromoters=gene.CpG$hypoPromoters, 
                    hyperPromoters=gene.CpG$hyperPromoters, 
                    hypoGenic=gene.CpG$hypoGenic, 
                    hyperGenic=gene.CpG$hyperGenic)),
      order.by="freq", sets=c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), keep.order=TRUE)
grid.text("Full model comparison of CpG hypo vs. hyper promoter vs. genic regions",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# aCpG
grid.newpage()
venn.diagram(x=list(hypoPromoters=gene.aCpG$hypoPromoters, 
                              hyperPromoters=gene.aCpG$hyperPromoters, 
                              hypoGenic=gene.aCpG$hypoGenic, 
                              hyperGenic=gene.aCpG$hyperGenic),
                       category.names = c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), 
                       main="Comparison of aCpG hypo vs. hyper promoter vs. genic regions",
             filename="4-Results/4.3-comparisons/promoterGenic/CpG-E-venn.tiff",imagetype="tiff",fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/promoterGenic/CpG-E-upset.pdf")
upset(fromList(list(hypoPromoters=gene.aCpG$hypoPromoters, 
                    hyperPromoters=gene.aCpG$hyperPromoters, 
                    hypoGenic=gene.aCpG$hypoGenic, 
                    hyperGenic=gene.aCpG$hyperGenic)),
      order.by="freq", sets=c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), keep.order=TRUE)
grid.text("Comparison of aCpG hypo vs. hyper promoter vs. genic regions",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# bCpG
grid.newpage()
venn.diagram(x=list(hypoPromoters=gene.bCpG$hypoPromoters, 
                              hyperPromoters=gene.bCpG$hyperPromoters, 
                              hypoGenic=gene.bCpG$hypoGenic, 
                              hyperGenic=gene.bCpG$hyperGenic),
                       category.names = c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), 
                       main="Comparison of bCpG hypo vs. hyper promoter vs. genic regions",
             filename="4-Results/4.3-comparisons/promoterGenic/CpG-G-venn.tiff",imagetype="tiff",fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/promoterGenic/CpG-G-upset.pdf")
upset(fromList(list(hypoPromoters=gene.bCpG$hypoPromoters, 
                    hyperPromoters=gene.bCpG$hyperPromoters, 
                    hypoGenic=gene.bCpG$hypoGenic, 
                    hyperGenic=gene.bCpG$hyperGenic)),
      order.by="freq", sets=c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), keep.order=TRUE)
grid.text("Comparison of bCpG hypo vs. hyper promoter vs. genic regions",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# cCpG
grid.newpage()
venn.diagram(x=list(hypoPromoters=gene.cCpG$hypoPromoters, 
                              hyperPromoters=gene.cCpG$hyperPromoters, 
                              hypoGenic=gene.cCpG$hypoGenic, 
                              hyperGenic=gene.cCpG$hyperGenic),
                       category.names = c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), 
                       main="Comparison of cCpG hypo vs. hyper promoter vs. genic regions",
             filename="4-Results/4.3-comparisons/promoterGenic/CpG-H-venn.tiff",imagetype="tiff",fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/promoterGenic/CpG-H-upset.pdf")
upset(fromList(list(hypoPromoters=gene.cCpG$hypoPromoters, 
                    hyperPromoters=gene.cCpG$hyperPromoters, 
                    hypoGenic=gene.cCpG$hypoGenic, 
                    hyperGenic=gene.cCpG$hyperGenic)),
      order.by="freq", sets=c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), keep.order=TRUE)
grid.text("Comparison of cCpG hypo vs. hyper promoter vs. genic regions",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()


# Histogram of gene:tile mappings -----------------------------------------
# Full model
ggplot(freq.CpG, aes(x=counts)) + 
  geom_histogram(binwidth=1, fill="white", color="black") +
  stat_bin(binwidth=1, geom='text', color='black', aes(label=..count..),position=position_dodge(0), vjust=-.6) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=max(freq.CpG$counts))) +
  labs(x="Number of tiles mapped to gene", y="Number of genes", title = "Frequency of tiles mapping to individual genes, full model, CpG context")
ggsave(filename="4-Results/4.3-comparisons/geneTileMap/fullCpGMap.pdf", plot=last_plot(),device="pdf")
ggplot(freq.aCpG, aes(x=counts)) + 
  geom_histogram(binwidth=1, fill="white", color="black") +
  stat_bin(binwidth=1, geom='text', color='black', aes(label=..count..),position=position_dodge(0), vjust=-.6) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=max(freq.aCpG$counts))) +
  labs(x="Number of tiles mapped to gene", y="Number of genes", title = "Frequency of tiles mapping to individual genes, A, CpG context")
ggsave(filename="4-Results/4.3-comparisons/geneTileMap/aCpGMap.pdf", plot=last_plot(),device="pdf")
ggplot(freq.bCpG, aes(x=counts)) + 
  geom_histogram(binwidth=1, fill="white", color="black") +
  stat_bin(binwidth=1, geom='text', color='black', aes(label=..count..),position=position_dodge(0), vjust=-.6) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=max(freq.bCpG$counts))) +
  labs(x="Number of tiles mapped to gene", y="Number of genes", title = "Frequency of tiles mapping to individual genes, B, CpG context")
ggsave(filename="4-Results/4.3-comparisons/geneTileMap/bCpGMap.pdf", plot=last_plot(),device="pdf")
ggplot(freq.cCpG, aes(x=counts)) + 
  geom_histogram(binwidth=1, fill="white", color="black") +
  stat_bin(binwidth=1, geom='text', color='black', aes(label=..count..),position=position_dodge(0), vjust=-.6) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=max(freq.cCpG$counts))) +
  labs(x="Number of tiles mapped to gene", y="Number of genes", title = "Frequency of tiles mapping to individual genes, C, CpG context")
ggsave(filename="4-Results/4.3-comparisons/geneTileMap/cCpGMap.pdf", plot=last_plot(),device="pdf")

# CHG ---------------------------------------------------------------------

# Hypo vs Hyper -----------------------------------------------------------
pdf(file="4-Results/4.3-comparisons/hyperHypo/CHG-hyperHypo-upset.pdf")
upset(fromList(list(
  fullHyper=hyp.CHG$hyper,
  AHyper=hyp.aCHG$hyper,
  BHyper=hyp.bCHG$hyper,
  CHyper=hyp.cCHG$hyper,
  fullHypo=hyp.CHG$hypo,
  AHypo=hyp.aCHG$hypo,
  BHypo=hyp.bCHG$hypo,
  CHypo=hyp.cCHG$hypo)), 
  mainbar.y.label = "Gene Intersections", sets.x.label = "DM genes per genotype", order.by="freq", number.angles=30,  text.scale = c(1.3, 1.3, 1, 1, 1, 0.75),
  sets=c("CHypo", "BHypo", "AHypo", "fullHypo","CHyper", "BHyper", "AHyper", "fullHyper"), keep.order=TRUE)
grid.text("Hyper- vs. hypomethylated genes across comparisons, CHG context",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

grid.newpage()
venn.diagram(x=list(unique(hyp.CHG$hyper),unique(hyp.cCHG$hyper),unique(hyp.aCHG$hyper),unique(hyp.bCHG$hyper)),
             category.names = c("full","C","A","B"), main="Hypermethylated genes across comparisons",
             filename="4-Results/4.3-comparisons/hyperHypo/CHG-hypermethylated.tiff",imagetype="tiff",
             fill=brewer.pal(4,"Set1"))
grid.newpage()
venn.diagram(x=list(unique(hyp.CHG$hypo),unique(hyp.cCHG$hypo),unique(hyp.aCHG$hypo),unique(hyp.bCHG$hypo)),
             category.names = c("full","C","A","B"), main="Hypomethylated genes across comparisons",
             filename="4-Results/4.3-comparisons/hyperHypo/CHG-hypomethylated.tiff",imagetype="tiff",          
             fill=brewer.pal(4,"Set1"))

grid.newpage()
venn.diagram(x=list(unique(hyp.CHG$hyper),unique(hyp.CHG$hypo)),
             category.names = c("hyper","hypo"), main="Hyper- vs hypomethylated genes in the full model",
             filename="4-Results/4.3-comparisons/hyperHypo/CHG-full-levelcomp.tiff",imagetype="tiff",
             fill=c("red", "blue"))
grid.newpage()
venn.diagram(x=list(unique(hyp.aCHG$hyper),unique(hyp.aCHG$hypo)),
             category.names = c("hyper","hypo"), main="Hyper- vs hypomethylated genes in A background",
             filename="4-Results/4.3-comparisons/hyperHypo/CHG-A-levelcomp.tiff",imagetype="tiff",
             fill=c("red", "blue"))
grid.newpage()
venn.diagram(x=list(unique(hyp.bCHG$hyper),unique(hyp.bCHG$hypo)),
             category.names = c("hyper","hypo"), main="Hyper- vs hypomethylated genes in B background",
             filename="4-Results/4.3-comparisons/hyperHypo/CHG-B-levelcomp.tiff",imagetype="tiff",
             fill=c("red", "blue"))
grid.newpage()
venn.diagram(x=list(unique(hyp.cCHG$hyper),unique(hyp.cCHG$hypo)),
             category.names = c("hyper","hypo"), main="Hyper- vs hypomethylated genes in C background",
             filename="4-Results/4.3-comparisons/hyperHypo/CHG-C-levelcomp.tiff",imagetype="tiff",
             fill=c("red", "blue"))

# Comparing genotype-specific DMRs ----------------------------------------
# Hypomethylated promoters
grid.newpage()
venn.diagram(x=list(gene.CHG$hypoPromoters, gene.cCHG$hypoPromoters,gene.aCHG$hypoPromoters,gene.bCHG$hypoPromoters),
             category.names = c("full","C","A","B"), main="Hypomethylated CHG promoters across comparisons",
             filename="4-Results/4.3-comparisons/genotypeSpecific/CHG-hypoProm-venn.tiff",imagetype="tiff",
             fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/genotypeSpecific/CHG-hypoProm-upset.pdf")
upset(fromList(list(full=gene.CHG$hypoPromoters, 
                    A=gene.aCHG$hypoPromoters, 
                    B=gene.bCHG$hypoPromoters, 
                    C=gene.cCHG$hypoPromoters)),
      order.by="freq", sets=c("C","B","A", "full"), keep.order=TRUE)
grid.text("Hypomethylated CHG promoters across comparisons",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# Hypermethylated promoters
grid.newpage()
venn.diagram(x=list(gene.CHG$hyperPromoters, gene.cCHG$hyperPromoters,gene.aCHG$hyperPromoters,gene.bCHG$hyperPromoters),
             category.names = c("full","C","A","B"), main="Hypermethylated CHG promoters across comparisons",
             filename="4-Results/4.3-comparisons/genotypeSpecific/CHG-hyperProm-venn.tiff",imagetype="tiff", fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/genotypeSpecific/CHG-hyperProm-upset.pdf")
upset(fromList(list(full=gene.CHG$hyperPromoters, 
                    A=gene.aCHG$hyperPromoters, 
                    B=gene.bCHG$hyperPromoters, 
                    C=gene.cCHG$hyperPromoters)),
      order.by="freq", sets=c("C","B","A", "full"), keep.order=TRUE)
grid.text("Hypermethylated CHG promoters across comparisons",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# Hypomethylated genic
grid.newpage()
venn.diagram(x=list(gene.CHG$hypoGenic, gene.cCHG$hypoGenic,gene.aCHG$hypoGenic,gene.bCHG$hypoGenic),
             category.names = c("full","C","A","B"), main="Hypomethylated CHG genic regions across comparisons",
             filename="4-Results/4.3-comparisons/genotypeSpecific/CHG-hypoGenic-venn.tiff",imagetype="tiff", fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/genotypeSpecific/CHG-hypoGenic-upset.pdf")
upset(fromList(list(full=gene.CHG$hypoGenic, 
                    A=gene.aCHG$hypoGenic, 
                    B=gene.bCHG$hypoGenic, 
                    C=gene.cCHG$hypoGenic)),
      order.by="freq", sets=c("C","B","A", "full"), keep.order=TRUE)
grid.text("Hypomethylated CHG genic regions across comparisons",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# Hypermethylated genic
grid.newpage()
venn.diagram(x=list(gene.CHG$hyperGenic, gene.cCHG$hyperGenic,gene.aCHG$hyperGenic,gene.bCHG$hyperGenic),
             category.names = c("full","C","A","B"), main="Hypermethylated CHG genic regions across comparisons",
             filename="4-Results/4.3-comparisons/genotypeSpecific/CHG-hyperGenic-venn.tiff",imagetype="tiff",fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/genotypeSpecific/CHG-hyperGenic-upset.pdf")
upset(fromList(list(full=gene.CHG$hyperGenic, 
                    A=gene.aCHG$hyperGenic, 
                    B=gene.bCHG$hyperGenic, 
                    C=gene.cCHG$hyperGenic)),
      order.by="freq", sets=c("C","B","A", "full"), keep.order=TRUE)
grid.text("Hypermethylated CHG genic regions across comparisons",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# Promoter vs. genic ------------------------------------------------------

# full
grid.newpage()
venn.diagram(x=list(hypoPromoters=gene.CHG$hypoPromoters, 
                    hyperPromoters=gene.CHG$hyperPromoters, 
                    hypoGenic=gene.CHG$hypoGenic, 
                    hyperGenic=gene.CHG$hyperGenic),
             category.names = c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), 
             main="Full model comparison of CHG hypo vs. hyper promoter vs. genic regions",
             filename="4-Results/4.3-comparisons/promoterGenic/CHG-full-venn.tiff",imagetype="tiff",fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/promoterGenic/CHG-full-upset.pdf")
upset(fromList(list(hypoPromoters=gene.CHG$hypoPromoters, 
                    hyperPromoters=gene.CHG$hyperPromoters, 
                    hypoGenic=gene.CHG$hypoGenic, 
                    hyperGenic=gene.CHG$hyperGenic)),
      order.by="freq", sets=c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), keep.order=TRUE)
grid.text("Full model comparison of CHG hypo vs. hyper promoter vs. genic regions",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# aCHG
grid.newpage()
venn.diagram(x=list(hypoPromoters=gene.aCHG$hypoPromoters, 
                    hyperPromoters=gene.aCHG$hyperPromoters, 
                    hypoGenic=gene.aCHG$hypoGenic, 
                    hyperGenic=gene.aCHG$hyperGenic),
             category.names = c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), 
             main="Comparison of aCHG hypo vs. hyper promoter vs. genic regions",
             filename="4-Results/4.3-comparisons/promoterGenic/CHG-E-venn.tiff",imagetype="tiff",fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/promoterGenic/CHG-E-upset.pdf")
upset(fromList(list(hypoPromoters=gene.aCHG$hypoPromoters, 
                    hyperPromoters=gene.aCHG$hyperPromoters, 
                    hypoGenic=gene.aCHG$hypoGenic, 
                    hyperGenic=gene.aCHG$hyperGenic)),
      order.by="freq", sets=c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), keep.order=TRUE)
grid.text("Comparison of aCHG hypo vs. hyper promoter vs. genic regions",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# bCHG
grid.newpage()
venn.diagram(x=list(hypoPromoters=gene.bCHG$hypoPromoters, 
                    hyperPromoters=gene.bCHG$hyperPromoters, 
                    hypoGenic=gene.bCHG$hypoGenic, 
                    hyperGenic=gene.bCHG$hyperGenic),
             category.names = c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), 
             main="Comparison of bCHG hypo vs. hyper promoter vs. genic regions",
             filename="4-Results/4.3-comparisons/promoterGenic/CHG-G-venn.tiff",imagetype="tiff",fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/promoterGenic/CHG-G-upset.pdf")
upset(fromList(list(hypoPromoters=gene.bCHG$hypoPromoters, 
                    hyperPromoters=gene.bCHG$hyperPromoters, 
                    hypoGenic=gene.bCHG$hypoGenic, 
                    hyperGenic=gene.bCHG$hyperGenic)),
      order.by="freq", sets=c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), keep.order=TRUE)
grid.text("Comparison of bCHG hypo vs. hyper promoter vs. genic regions",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# cCHG
grid.newpage()
venn.diagram(x=list(hypoPromoters=gene.cCHG$hypoPromoters, 
                    hyperPromoters=gene.cCHG$hyperPromoters, 
                    hypoGenic=gene.cCHG$hypoGenic, 
                    hyperGenic=gene.cCHG$hyperGenic),
             category.names = c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), 
             main="Comparison of cCHG hypo vs. hyper promoter vs. genic regions",
             filename="4-Results/4.3-comparisons/promoterGenic/CHG-H-venn.tiff",imagetype="tiff",fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/promoterGenic/CHG-H-upset.pdf")
upset(fromList(list(hypoPromoters=gene.cCHG$hypoPromoters, 
                    hyperPromoters=gene.cCHG$hyperPromoters, 
                    hypoGenic=gene.cCHG$hypoGenic, 
                    hyperGenic=gene.cCHG$hyperGenic)),
      order.by="freq", sets=c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), keep.order=TRUE)
grid.text("Comparison of cCHG hypo vs. hyper promoter vs. genic regions",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()


# Histogram of gene:tile mappings -----------------------------------------
# Full model
ggplot(freq.CHG, aes(x=counts)) + 
  geom_histogram(binwidth=1, fill="white", color="black") +
  stat_bin(binwidth=1, geom='text', color='black', aes(label=..count..),position=position_dodge(0), vjust=-.6) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=max(freq.CHG$counts))) +
  labs(x="Number of tiles mapped to gene", y="Number of genes", title = "Frequency of tiles mapping to individual genes, full model, CHG context")
ggsave(filename="4-Results/4.3-comparisons/geneTileMap/fullCHGMap.pdf", plot=last_plot(),device="pdf")
ggplot(freq.aCHG, aes(x=counts)) + 
  geom_histogram(binwidth=1, fill="white", color="black") +
  stat_bin(binwidth=1, geom='text', color='black', aes(label=..count..),position=position_dodge(0), vjust=-.6) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=max(freq.aCHG$counts))) +
  labs(x="Number of tiles mapped to gene", y="Number of genes", title = "Frequency of tiles mapping to individual genes, A, CHG context")
ggsave(filename="4-Results/4.3-comparisons/geneTileMap/aCHGMap.pdf", plot=last_plot(),device="pdf")
ggplot(freq.bCHG, aes(x=counts)) + 
  geom_histogram(binwidth=1, fill="white", color="black") +
  stat_bin(binwidth=1, geom='text', color='black', aes(label=..count..),position=position_dodge(0), vjust=-.6) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=max(freq.bCHG$counts))) +
  labs(x="Number of tiles mapped to gene", y="Number of genes", title = "Frequency of tiles mapping to individual genes, B, CHG context")
ggsave(filename="4-Results/4.3-comparisons/geneTileMap/bCHGMap.pdf", plot=last_plot(),device="pdf")
ggplot(freq.cCHG, aes(x=counts)) + 
  geom_histogram(binwidth=1, fill="white", color="black") +
  stat_bin(binwidth=1, geom='text', color='black', aes(label=..count..),position=position_dodge(0), vjust=-.6) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=max(freq.cCHG$counts))) +
  labs(x="Number of tiles mapped to gene", y="Number of genes", title = "Frequency of tiles mapping to individual genes, C, CHG context")
ggsave(filename="4-Results/4.3-comparisons/geneTileMap/cCHGMap.pdf", plot=last_plot(),device="pdf")

# CHH ---------------------------------------------------------------------

# Hypo vs Hyper -----------------------------------------------------------
pdf(file="4-Results/4.3-comparisons/hyperHypo/CHH-hyperHypo-upset.pdf")
upset(fromList(list(
  AHyper=hyp.aCHH$hyper,
  BHyper=hyp.bCHH$hyper,
  CHyper=hyp.cCHH$hyper,
  fullHypo=hyp.CHH$hypo,
  AHypo=hyp.aCHH$hypo,
  BHypo=hyp.bCHH$hypo,
  CHypo=hyp.cCHH$hypo)), 
  mainbar.y.label = "Gene Intersections", sets.x.label = "DM genes per genotype", order.by="freq", number.angles=30,  text.scale = c(1.3, 1.3, 1, 1, 1, 0.75),
  sets=c("CHypo", "BHypo", "AHypo", "fullHypo","CHyper", "BHyper", "AHyper"), keep.order=TRUE)
grid.text("Hyper- vs. hypomethylated genes across comparisons, CHH context",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

grid.newpage()
venn.diagram(x=list(unique(hyp.CHH$hyper),unique(hyp.cCHH$hyper),unique(hyp.aCHH$hyper),unique(hyp.bCHH$hyper)),
             category.names = c("full","C","A","B"), main="Hypermethylated genes across comparisons",
             filename="4-Results/4.3-comparisons/hyperHypo/CHH-hypermethylated.tiff",imagetype="tiff",
             fill=brewer.pal(4,"Set1"))
grid.newpage()
venn.diagram(x=list(unique(hyp.CHH$hypo),unique(hyp.cCHH$hypo),unique(hyp.aCHH$hypo),unique(hyp.bCHH$hypo)),
             category.names = c("full","C","A","B"), main="Hypomethylated genes across comparisons",
             filename="4-Results/4.3-comparisons/hyperHypo/CHH-hypomethylated.tiff",imagetype="tiff",          
             fill=brewer.pal(4,"Set1"))

grid.newpage()
venn.diagram(x=list(unique(hyp.CHH$hyper),unique(hyp.CHH$hypo)),
             category.names = c("hyper","hypo"), main="Hyper- vs hypomethylated genes in the full model",
             filename="4-Results/4.3-comparisons/hyperHypo/CHH-full-levelcomp.tiff",imagetype="tiff",
             fill=c("red", "blue"))
grid.newpage()
venn.diagram(x=list(unique(hyp.aCHH$hyper),unique(hyp.aCHH$hypo)),
             category.names = c("hyper","hypo"), main="Hyper- vs hypomethylated genes in A background",
             filename="4-Results/4.3-comparisons/hyperHypo/CHH-A-levelcomp.tiff",imagetype="tiff",
             fill=c("red", "blue"))
grid.newpage()
venn.diagram(x=list(unique(hyp.bCHH$hyper),unique(hyp.bCHH$hypo)),
             category.names = c("hyper","hypo"), main="Hyper- vs hypomethylated genes in B background",
             filename="4-Results/4.3-comparisons/hyperHypo/CHH-B-levelcomp.tiff",imagetype="tiff",
             fill=c("red", "blue"))
grid.newpage()
venn.diagram(x=list(unique(hyp.cCHH$hyper),unique(hyp.cCHH$hypo)),
             category.names = c("hyper","hypo"), main="Hyper- vs hypomethylated genes in C background",
             filename="4-Results/4.3-comparisons/hyperHypo/CHH-C-levelcomp.tiff",imagetype="tiff",
             fill=c("red", "blue"))

# Comparing genotype-specific DMRs ----------------------------------------
# Hypomethylated promoters
grid.newpage()
venn.diagram(x=list(gene.CHH$hypoPromoters, gene.cCHH$hypoPromoters,gene.aCHH$hypoPromoters,gene.bCHH$hypoPromoters),
             category.names = c("full","C","A","B"), main="Hypomethylated CHH promoters across comparisons",
             filename="4-Results/4.3-comparisons/genotypeSpecific/CHH-hypoProm-venn.tiff",imagetype="tiff",
             fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/genotypeSpecific/CHH-hypoProm-upset.pdf")
upset(fromList(list(full=gene.CHH$hypoPromoters, 
                    A=gene.aCHH$hypoPromoters, 
                    B=gene.bCHH$hypoPromoters, 
                    C=gene.cCHH$hypoPromoters)),
      order.by="freq", sets=c("C","B","A", "full"), keep.order=TRUE)
grid.text("Hypomethylated CHH promoters across comparisons",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# Hypermethylated promoters
grid.newpage()
venn.diagram(x=list(gene.CHH$hyperPromoters, gene.cCHH$hyperPromoters,gene.aCHH$hyperPromoters,gene.bCHH$hyperPromoters),
             category.names = c("full","C","A","B"), main="Hypermethylated CHH promoters across comparisons",
             filename="4-Results/4.3-comparisons/genotypeSpecific/CHH-hyperProm-venn.tiff",imagetype="tiff", fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/genotypeSpecific/CHH-hyperProm-upset.pdf")
upset(fromList(list(A=gene.aCHH$hyperPromoters, 
                    B=gene.bCHH$hyperPromoters, 
                    C=gene.cCHH$hyperPromoters)),
      order.by="freq", sets=c("C","B","A"), keep.order=TRUE)
grid.text("Hypermethylated CHH promoters across comparisons",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# Hypomethylated genic
grid.newpage()
venn.diagram(x=list(gene.CHH$hypoGenic, gene.cCHH$hypoGenic,gene.aCHH$hypoGenic,gene.bCHH$hypoGenic),
             category.names = c("full","C","A","B"), main="Hypomethylated CHH genic regions across comparisons",
             filename="4-Results/4.3-comparisons/genotypeSpecific/CHH-hypoGenic-venn.tiff",imagetype="tiff", fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/genotypeSpecific/CHH-hypoGenic-upset.pdf")
upset(fromList(list(full=gene.CHH$hypoGenic, 
                    A=gene.aCHH$hypoGenic, 
                    B=gene.bCHH$hypoGenic, 
                    C=gene.cCHH$hypoGenic)),
      order.by="freq", sets=c("C","B","A", "full"), keep.order=TRUE)
grid.text("Hypomethylated CHH genic regions across comparisons",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# Hypermethylated genic
grid.newpage()
venn.diagram(x=list(gene.CHH$hyperGenic, gene.cCHH$hyperGenic,gene.aCHH$hyperGenic,gene.bCHH$hyperGenic),
             category.names = c("full","C","A","B"), main="Hypermethylated CHH genic regions across comparisons",
             filename="4-Results/4.3-comparisons/genotypeSpecific/CHH-hyperGenic-venn.tiff",imagetype="tiff",fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/genotypeSpecific/CHH-hyperGenic-upset.pdf")
upset(fromList(list(A=gene.aCHH$hyperGenic, 
                    B=gene.bCHH$hyperGenic, 
                    C=gene.cCHH$hyperGenic)),
      order.by="freq", sets=c("C","B","A"), keep.order=TRUE)
grid.text("Hypermethylated CHH genic regions across comparisons",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# Promoter vs. genic ------------------------------------------------------

# full
grid.newpage()
venn.diagram(x=list(hypoPromoters=gene.CHH$hypoPromoters, 
                    hyperPromoters=gene.CHH$hyperPromoters, 
                    hypoGenic=gene.CHH$hypoGenic, 
                    hyperGenic=gene.CHH$hyperGenic),
             category.names = c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), 
             main="Full model comparison of CHH hypo vs. hyper promoter vs. genic regions",
             filename="4-Results/4.3-comparisons/promoterGenic/CHH-full-venn.tiff",imagetype="tiff",fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/promoterGenic/CHH-full-upset.pdf")
upset(fromList(list(hypoPromoters=gene.CHH$hypoPromoters, 
                    hyperPromoters=gene.CHH$hyperPromoters, 
                    hypoGenic=gene.CHH$hypoGenic, 
                    hyperGenic=gene.CHH$hyperGenic)),
      order.by="freq", sets=c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), keep.order=TRUE)
grid.text("Full model comparison of CHH hypo vs. hyper promoter vs. genic regions",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# aCHH
grid.newpage()
venn.diagram(x=list(hypoPromoters=gene.aCHH$hypoPromoters, 
                    hyperPromoters=gene.aCHH$hyperPromoters, 
                    hypoGenic=gene.aCHH$hypoGenic, 
                    hyperGenic=gene.aCHH$hyperGenic),
             category.names = c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), 
             main="Comparison of aCHH hypo vs. hyper promoter vs. genic regions",
             filename="4-Results/4.3-comparisons/promoterGenic/CHH-E-venn.tiff",imagetype="tiff",fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/promoterGenic/CHH-E-upset.pdf")
upset(fromList(list(hypoPromoters=gene.aCHH$hypoPromoters, 
                    hyperPromoters=gene.aCHH$hyperPromoters, 
                    hypoGenic=gene.aCHH$hypoGenic, 
                    hyperGenic=gene.aCHH$hyperGenic)),
      order.by="freq", sets=c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), keep.order=TRUE)
grid.text("Comparison of aCHH hypo vs. hyper promoter vs. genic regions",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# bCHH
grid.newpage()
venn.diagram(x=list(hypoPromoters=gene.bCHH$hypoPromoters, 
                    hyperPromoters=gene.bCHH$hyperPromoters, 
                    hypoGenic=gene.bCHH$hypoGenic, 
                    hyperGenic=gene.bCHH$hyperGenic),
             category.names = c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), 
             main="Comparison of bCHH hypo vs. hyper promoter vs. genic regions",
             filename="4-Results/4.3-comparisons/promoterGenic/CHH-G-venn.tiff",imagetype="tiff",fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/promoterGenic/CHH-G-upset.pdf")
upset(fromList(list(hypoPromoters=gene.bCHH$hypoPromoters, 
                    hyperPromoters=gene.bCHH$hyperPromoters, 
                    hypoGenic=gene.bCHH$hypoGenic, 
                    hyperGenic=gene.bCHH$hyperGenic)),
      order.by="freq", sets=c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), keep.order=TRUE)
grid.text("Comparison of bCHH hypo vs. hyper promoter vs. genic regions",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

# cCHH
grid.newpage()
venn.diagram(x=list(hypoPromoters=gene.cCHH$hypoPromoters, 
                    hyperPromoters=gene.cCHH$hyperPromoters, 
                    hypoGenic=gene.cCHH$hypoGenic, 
                    hyperGenic=gene.cCHH$hyperGenic),
             category.names = c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), 
             main="Comparison of cCHH hypo vs. hyper promoter vs. genic regions",
             filename="4-Results/4.3-comparisons/promoterGenic/CHH-H-venn.tiff",imagetype="tiff",fill=brewer.pal(4,"Set1"))

pdf(file="4-Results/4.3-comparisons/promoterGenic/CHH-H-upset.pdf")
upset(fromList(list(hypoPromoters=gene.cCHH$hypoPromoters, 
                    hyperPromoters=gene.cCHH$hyperPromoters, 
                    hypoGenic=gene.cCHH$hypoGenic, 
                    hyperGenic=gene.cCHH$hyperGenic)),
      order.by="freq", sets=c("hypoPromoters", "hyperPromoters","hypoGenic","hyperGenic"), keep.order=TRUE)
grid.text("Comparison of cCHH hypo vs. hyper promoter vs. genic regions",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()


# Histogram of gene:tile mappings -----------------------------------------
# Full model
ggplot(freq.CHH, aes(x=counts)) + 
  geom_histogram(binwidth=1, fill="white", color="black") +
  stat_bin(binwidth=1, geom='text', color='black', aes(label=..count..),position=position_dodge(0), vjust=-.6) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=max(freq.CHH$counts))) +
  labs(x="Number of tiles mapped to gene", y="Number of genes", title = "Frequency of tiles mapping to individual genes, full model, CHH context")
ggsave(filename="4-Results/4.3-comparisons/geneTileMap/fullCHHMap.pdf", plot=last_plot(),device="pdf")
ggplot(freq.aCHH, aes(x=counts)) + 
  geom_histogram(binwidth=1, fill="white", color="black") +
  stat_bin(binwidth=1, geom='text', color='black', aes(label=..count..),position=position_dodge(0), vjust=-.6) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=max(freq.aCHH$counts))) +
  labs(x="Number of tiles mapped to gene", y="Number of genes", title = "Frequency of tiles mapping to individual genes, A, CHH context")
ggsave(filename="4-Results/4.3-comparisons/geneTileMap/aCHHMap.pdf", plot=last_plot(),device="pdf")
ggplot(freq.bCHH, aes(x=counts)) + 
  geom_histogram(binwidth=1, fill="white", color="black") +
  stat_bin(binwidth=1, geom='text', color='black', aes(label=..count..),position=position_dodge(0), vjust=-.6) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=max(freq.bCHH$counts))) +
  labs(x="Number of tiles mapped to gene", y="Number of genes", title = "Frequency of tiles mapping to individual genes, B, CHH context")
ggsave(filename="4-Results/4.3-comparisons/geneTileMap/bCHHMap.pdf", plot=last_plot(),device="pdf")
ggplot(freq.cCHH, aes(x=counts)) + 
  geom_histogram(binwidth=1, fill="white", color="black") +
  stat_bin(binwidth=1, geom='text', color='black', aes(label=..count..),position=position_dodge(0), vjust=-.6) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=max(freq.cCHH$counts))) +
  labs(x="Number of tiles mapped to gene", y="Number of genes", title = "Frequency of tiles mapping to individual genes, C, CHH context")
ggsave(filename="4-Results/4.3-comparisons/geneTileMap/cCHHMap.pdf", plot=last_plot(),device="pdf")

# Overlap between genotypes ----------------------------------------------
grid.newpage()
venn.diagram(x=list(hyp.aCpG$hyper,hyp.bCpG$hyper,hyp.cCpG$hyper), force.unique=TRUE,
             category.names = c("A","B","C"), 
             filename="4-Results/4.3-comparisons/hyperHypo/3hyperCpG.tiff",imagetype="tiff",
             fill=brewer.pal(3,"Set1"),print.mode=c("raw","percent"))
grid.newpage()
venn.diagram(x=list(hyp.aCpG$hypo,hyp.bCpG$hypo,hyp.cCpG$hypo), force.unique=TRUE,
             category.names =  c("A","B","C"), 
             filename="4-Results/4.3-comparisons/hyperHypo/3hypoCpG.tiff",imagetype="tiff",          
             fill=brewer.pal(3,"Set1"),print.mode=c("raw","percent"))
grid.newpage()

venn.diagram(x=list(hyp.aCHG$hyper,hyp.bCHG$hyper,hyp.cCHG$hyper), force.unique=TRUE,
             category.names = c("A","B","C"), 
             filename="4-Results/4.3-comparisons/hyperHypo/3hyperCHG.tiff",imagetype="tiff",
             fill=brewer.pal(3,"Set1"),print.mode=c("raw","percent"))
grid.newpage()
venn.diagram(x=list(hyp.aCHG$hypo,hyp.bCHG$hypo,hyp.cCHG$hypo), force.unique=TRUE,
             category.names =  c("A","B","C"), 
             filename="4-Results/4.3-comparisons/hyperHypo/3hypoCHG.tiff",imagetype="tiff",          
             fill=brewer.pal(3,"Set1"),print.mode=c("raw","percent"))
grid.newpage()

venn.diagram(x=list(hyp.aCHH$hyper,hyp.bCHH$hyper,hyp.cCHH$hyper), force.unique=TRUE,
             category.names = c("A","B","C"), 
             filename="4-Results/4.3-comparisons/hyperHypo/3hyperCHH.tiff",imagetype="tiff",
             fill=brewer.pal(3,"Set1"),print.mode=c("raw","percent"))
grid.newpage()
venn.diagram(x=list(hyp.aCHH$hypo,hyp.bCHH$hypo,hyp.cCHH$hypo), force.unique=TRUE,
             category.names =  c("A","B","C"), 
             filename="4-Results/4.3-comparisons/hyperHypo/3hypoCHH.tiff",imagetype="tiff",          
             fill=brewer.pal(3,"Set1"),print.mode=c("raw","percent"))
grid.newpage()

# Gene:tile maps ----------------------------------------------------------
freqList <- list(freq.aCpG, freq.bCpG, freq.cCpG, freq.aCHG, freq.bCHG, freq.cCHG, freq.aCHH, freq.bCHH, freq.cCHH)
freqNames <- list("A CpG", "B CpG", "C CpG", "A CHG", "B CHG", "C CHG", "A CHH", "B CHH", "C CHH")
freqPlot <- Map(function(x,y) {
  
  ggplot(x, aes(x=counts)) + 
    geom_histogram(binwidth=1, fill="white", color="black") +
    stat_bin(binwidth=1, geom='text', color='black', aes(label=..count..), angle=90, hjust=-0.3) +
    ylim(0,8000) +
    labs(x="# tiles", y="# genes", title = y)+
    coord_cartesian(xlim = c(1, 14))
  
}, freqList, freqNames)

plot_grid(freqPlot[[1]], freqPlot[[2]], freqPlot[[3]], freqPlot[[4]], freqPlot[[5]], freqPlot[[6]], freqPlot[[7]], freqPlot[[8]], freqPlot[[9]],
          ncol=3, labels="AUTO")
