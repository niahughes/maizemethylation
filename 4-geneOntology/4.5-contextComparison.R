###################################################################
### COMPARISONS ACROSS CONTEXTS
### Local R script 
###################################################################

# set up environment
setwd("D:/Nia/OneDrive - University of Guelph/Working/DMRanalysis")
source("D:/Nia/OneDrive - University of Guelph/Working/thesisScripts/4-geneOntology/4.1-dataPrep.R")
library(UpSetR)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")


# Comparison-specific -----------------------------------------------------
# Full model
pdf(file="4-Results/4.5-contextComparison/comparisonSpecific/fullContextUpset.pdf")
upset(fromList(list(
  fullHyperCpG=hyp.CpG$hyper,
  fullHyperCHG=hyp.CHG$hyper,
  fullHyperCHH=hyp.CHH$hyper,
  fullHypoCpG=hyp.CpG$hypo,
  fullHypoCHG=hyp.CHG$hypo,
  fullHypoCHH=hyp.CHH$hypo
)), 
  mainbar.y.label = "Gene Intersections", sets.x.label = "DM genes per context", order.by="freq", number.angles=30,  text.scale = c(1.3, 1.3, 1, 1, 1, 0.75),
  sets=c("fullHypoCHH","fullHypoCHG","fullHypoCpG", "fullHyperCHH","fullHyperCHG","fullHyperCpG"), keep.order=TRUE)
grid.text("Hyper- vs. hypomethylated genes across sequence contexts, full model",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

#hypomethylation
grid.newpage()
venn.diagram(x=list(hyp.CpG$hypo, hyp.CHG$hypo, hyp.CHH$hypo),
             category.names = c("CpG","CHG","CHH"), main="Hypomethylated genes across contexts, full model",
             filename="4-Results/4.5-contextComparison/comparisonSpecific/fullContextHypoVenn.tiff",imagetype="tiff",
             fill=brewer.pal(3,"Set1"))
#hypermethylation
grid.newpage()
venn.diagram(x=list(hyp.CpG$hyper, hyp.CHG$hyper, hyp.CHH$hyper),
                       category.names = c("CpG","CHG","CHH"), main="Hypermethylated genes across contexts, full model",
                       filename="4-Results/4.5-contextComparison/comparisonSpecific/fullContextHyperVenn.tiff",imagetype="tiff",
                       fill=brewer.pal(3,"Set1"))

# A
pdf(file="4-Results/4.5-contextComparison/comparisonSpecific/eContextUpset.pdf")
upset(fromList(list(
  eHyperCpG=hyp.aCpG$hyper,
  eHyperCHG=hyp.aCHG$hyper,
  eHyperCHH=hyp.aCHH$hyper,
  eHypoCpG=hyp.aCpG$hypo,
  eHypoCHG=hyp.aCHG$hypo,
  eHypoCHH=hyp.aCHH$hypo
)), 
mainbar.y.label = "Gene Intersections", sets.x.label = "DM genes per context", order.by="freq", number.angles=30,  text.scale = c(1.3, 1.3, 1, 1, 1, 0.75),
sets=c("eHypoCHH","eHypoCHG","eHypoCpG", "eHyperCHH","eHyperCHG","eHyperCpG"), keep.order=TRUE)
grid.text("Hyper- vs. hypomethylated genes across sequence contexts, A",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

#hypomethylation
grid.newpage()
venn.diagram(x=list(hyp.aCpG$hypo, hyp.aCHG$hypo, hyp.aCHH$hypo),
                       category.names = c("CpG","CHG","CHH"), main="Hypomethylated genes across contexts, A",
                       filename="4-Results/4.5-contextComparison/comparisonSpecific/eContextHypoVenn.tiff",imagetype="tiff",
                       fill=brewer.pal(3,"Set1"))
#hypermethylation
grid.newpage()
venn.diagram(x=list(hyp.aCpG$hyper, hyp.aCHG$hyper, hyp.aCHH$hyper),
                       category.names = c("CpG","CHG","CHH"), main="Hypermethylated genes across contexts, A",
                       filename="4-Results/4.5-contextComparison/comparisonSpecific/eContextHyperVenn.tiff",imagetype="tiff",
                       fill=brewer.pal(3,"Set1"))

# B
pdf(file="4-Results/4.5-contextComparison/comparisonSpecific/gContextUpset.pdf")
upset(fromList(list(
  gHyperCpG=hyp.bCpG$hyper,
  gHyperCHG=hyp.bCHG$hyper,
  gHyperCHH=hyp.bCHH$hyper,
  gHypoCpG=hyp.bCpG$hypo,
  gHypoCHG=hyp.bCHG$hypo,
  gHypoCHH=hyp.bCHH$hypo
)), 
mainbar.y.label = "Gene Intersections", sets.x.label = "DM genes per context", order.by="freq", number.angles=30,  text.scale = c(1.3, 1.3, 1, 1, 1, 0.75),
sets=c("gHypoCHH","gHypoCHG","gHypoCpG", "gHyperCHH","gHyperCHG","gHyperCpG"), keep.order=TRUE)
grid.text("Hyper- vs. hypomethylated genes across sequence contexts, B",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

#hypomethylation
grid.newpage()
venn.diagram(x=list(hyp.bCpG$hypo, hyp.bCHG$hypo, hyp.bCHH$hypo),
                       category.names = c("CpG","CHG","CHH"), main="Hypomethylated genes across contexts, B",
                       filename="4-Results/4.5-contextComparison/comparisonSpecific/gContextHypoVenn.tiff",imagetype="tiff",
                       fill=brewer.pal(3,"Set1"))
#hypermethylation
grid.newpage()
venn.diagram(x=list(hyp.bCpG$hyper, hyp.bCHG$hyper, hyp.bCHH$hyper),
                       category.names = c("CpG","CHG","CHH"), main="Hypermethylated genes across contexts, B",
                       filename="4-Results/4.5-contextComparison/comparisonSpecific/gContextHyperVenn.tiff",imagetype="tiff",
                       fill=brewer.pal(3,"Set1"))

# C
pdf(file="4-Results/4.5-contextComparison/comparisonSpecific/hContextUpset.pdf")
upset(fromList(list(
  hHyperCpG=hyp.cCpG$hyper,
  hHyperCHG=hyp.cCHG$hyper,
  hHyperCHH=hyp.cCHH$hyper,
  hHypoCpG=hyp.cCpG$hypo,
  hHypoCHG=hyp.cCHG$hypo,
  hHypoCHH=hyp.cCHH$hypo
)), 
mainbar.y.label = "Gene Intersections", sets.x.label = "DM genes per context", order.by="freq", number.angles=30,  text.scale = c(1.3, 1.3, 1, 1, 1, 0.75),
sets=c("hHypoCHH","hHypoCHG","hHypoCpG", "hHyperCHH","hHyperCHG","hHyperCpG"), keep.order=TRUE)
grid.text("Hyper- vs. hypomethylated genes across sequence contexts, C",x = 0.65, y=0.95, gp=gpar(fontsize=8))
dev.off()

#hypomethylation
grid.newpage()
venn.diagram(x=list(hyp.cCpG$hypo, hyp.cCHG$hypo, hyp.cCHH$hypo),
                       category.names = c("CpG","CHG","CHH"), main="Hypomethylated genes across contexts, C",
                       filename="4-Results/4.5-contextComparison/comparisonSpecific/hContextHypoVenn.tiff",imagetype="tiff",
                       fill=brewer.pal(3,"Set1"))
#hypermethylation
grid.newpage()
venn.diagram(x=list(hyp.cCpG$hyper, hyp.cCHG$hyper, hyp.cCHH$hyper),
                       category.names = c("CpG","CHG","CHH"), main="Hypermethylated genes across contexts, C",
                       filename="4-Results/4.5-contextComparison/comparisonSpecific/hContextHyperVenn.tiff",imagetype="tiff",
                       fill=brewer.pal(3,"Set1"))



# Pooled ------------------------------------------------------------------
# hypomethylation
allHypoCpG <- unique(do.call(c, list(as.character(hyp.CpG$hypo), as.character(hyp.aCpG$hypo),
                                  as.character(hyp.bCpG$hypo), as.character(hyp.cCpG$hypo))))
allHypoCHG <- unique(do.call(c, list(as.character(hyp.CHG$hypo), as.character(hyp.aCHG$hypo),
                                  as.character(hyp.bCHG$hypo), as.character(hyp.cCHG$hypo))))
allHypoCHH <- unique(do.call(c, list(as.character(hyp.CHH$hypo), as.character(hyp.aCHH$hypo),
                                  as.character(hyp.bCHH$hypo), as.character(hyp.cCHH$hypo))))
grid.newpage()
venn.diagram(x=list(allHypoCpG, allHypoCHG, allHypoCHH),
                       category.names = c("CpG","CHG","CHH"), main="Hypomethylated genes across sequence contexts",
                       filename="4-Results/4.5-contextComparison/pooled/hypoVenn.tiff",imagetype="tiff",
                       fill=brewer.pal(3,"Set1"))

# hypermethylation
allHyperCpG <- unique(do.call(c, list(as.character(hyp.CpG$hyper), as.character(hyp.aCpG$hyper),
                                     as.character(hyp.bCpG$hyper), as.character(hyp.cCpG$hyper))))
allHyperCHG <- unique(do.call(c, list(as.character(hyp.CHG$hyper), as.character(hyp.aCHG$hyper),
                                     as.character(hyp.bCHG$hyper), as.character(hyp.cCHG$hyper))))
allHyperCHH <- unique(do.call(c, list(as.character(hyp.CHH$hyper), as.character(hyp.aCHH$hyper),
                                     as.character(hyp.bCHH$hyper), as.character(hyp.cCHH$hyper))))
grid.newpage()
venn.diagram(x=list(allHyperCpG, allHyperCHG, allHyperCHH),
                       category.names = c("CpG","CHG","CHH"), main="Hypermethylated genes across sequence contexts",
                       filename="4-Results/4.5-contextComparison/pooled/hyperVenn.tiff",imagetype="tiff",
                       fill=brewer.pal(3,"Set1"))