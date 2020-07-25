###################################################################
### QUANTIFYING TILED DATA
### R script for miscellaneous quantifications of tiled data
###################################################################

# set up environment
setwd("/scratch/nia/manuFinal")
load(file=".RData")
library(methylKit)
library(GenomicRanges)
library(genomation)
library(reshape2)

# Overdispersion correction: how does this affect the number of DM tiles?
# CpG
tileCountCpG <- c(CpG.3x.100.nocorr.50p.all@num.records, CpG.3x.100.oChi.50p.all@num.records, CpG.3x.100.oF.50p.all@num.records)
barplot(tileCountCpG, main = "Number of DM tiles by overdispersion correction, CpG context", sub="CpG context, >|50|% methylation difference, q < 0.01",
xlab = "Overdispersion Correction Type",
ylab = "Number of DM tiles",
names.arg = c("uncorrected", "Chi-squared", "F-test"),
col = "darkred")
# CHG
tileCountCHG <- c(CHG.3x.100.nocorr.50p.all@num.records, CHG.3x.100.oChi.50p.all@num.records, CHG.3x.100.oF.50p.all@num.records)
barplot(tileCountCHG, main = "Number of DM tiles by overdispersion correction, CHG context", sub="CHG context, >|50|% methylation difference, q < 0.01",
xlab = "Overdispersion Correction Type",
ylab = "Number of DM tiles",
names.arg = c("uncorrected", "Chi-squared", "F-test"),
col = "darkred")
# CHH
tileCountCHH <- c(CHH.3x.100.nocorr.50p.all@num.records, CHH.3x.100.oChi.50p.all@num.records, CHH.3x.100.oF.50p.all@num.records)
barplot(tileCountCHH, main = "Number of DM tiles by overdispersion correction, CHH context", sub="CHH context, >|50|% methylation difference, q < 0.01",
xlab = "Overdispersion Correction Type",
ylab = "Number of DM tiles",
names.arg = c("uncorrected", "Chi-squared", "F-test"),
col = "darkred")

#What percent of tiles are differentially methylated?
#CpG
CpG <- getData(CpG.3x.100.oChi.50p.all)
aCpG <- getData(aCpG.3x.100.50p.all)
bCpG <- getData(bCpG.3x.100.50p.all)
cCpG <- getData(cCpG.3x.100.50p.all)

WTvMutCpG <- data.frame("model"=c("full model", "A", "B", "C"),
                        "totaltiles"=c(diff.CpG.3x.100.oChi@num.records, diff.aCpG.3x.100@num.records, diff.bCpG.3x.100@num.records, diff.cCpG.3x.100@num.records),
                        "DMtiles"=c(CpG.3x.100.oChi.50p.all@num.records, aCpG.3x.100.50p.all@num.records, bCpG.3x.100.50p.all@num.records, cCpG.3x.100.50p.all@num.records),
                        "hypo"=c(nrow(CpG[CpG$meth.diff<0,]), nrow(aCpG[aCpG$meth.diff<0,]),nrow(bCpG[bCpG$meth.diff<0,]), nrow(cCpG[cCpG$meth.diff<0,])),
                        "hyper"=c(nrow(CpG[CpG$meth.diff>0,]), nrow(aCpG[aCpG$meth.diff>0,]),nrow(bCpG[bCpG$meth.diff>0,]), nrow(cCpG[cCpG$meth.diff>0,])))

WTvMutCpG$percent <- round(WTvMutCpG$DMtiles/WTvMutCpG$totaltiles*100,3)
WTvMutCpG$hypoPerc <- round(WTvMutCpG$hypo/WTvMutCpG$totaltiles*100,3)
WTvMutCpG$hyperPerc <- round(WTvMutCpG$hyper/WTvMutCpG$totaltiles*100,3)
print(WTvMutCpG)
#CHG
CHG <- getData(CHG.3x.100.oChi.50p.all)
aCHG <- getData(aCHG.3x.100.50p.all)
bCHG <- getData(bCHG.3x.100.50p.all)
cCHG <- getData(cCHG.3x.100.50p.all)

WTvMutCHG <- data.frame("model"=c("full model","A", "B", "C"),
                        "totaltiles"=c(diff.CHG.3x.100.oChi@num.records, diff.aCHG.3x.100@num.records, diff.bCHG.3x.100@num.records, diff.cCHG.3x.100@num.records),
                        "DMtiles"=c(CHG.3x.100.oChi.50p.all@num.records, aCHG.3x.100.50p.all@num.records, bCHG.3x.100.50p.all@num.records, cCHG.3x.100.50p.all@num.records),
                        "hypo"=c(nrow(CHG[CHG$meth.diff<0,]), nrow(aCHG[aCHG$meth.diff<0,]),nrow(bCHG[bCHG$meth.diff<0,]), nrow(cCHG[cCHG$meth.diff<0,])),
                        "hyper"=c(nrow(CHG[CHG$meth.diff>0,]), nrow(aCHG[aCHG$meth.diff>0,]),nrow(bCHG[bCHG$meth.diff>0,]), nrow(cCHG[cCHG$meth.diff>0,])))

WTvMutCHG$percent <- round(WTvMutCHG$DMtiles/WTvMutCHG$totaltiles*100,3)
WTvMutCHG$hypoPerc <- round(WTvMutCHG$hypo/WTvMutCHG$totaltiles*100,3)
WTvMutCHG$hyperPerc <- round(WTvMutCHG$hyper/WTvMutCHG$totaltiles*100,3)
print(WTvMutCHG)
#CHH
CHH <- getData(CHH.3x.100.oChi.50p.all)
aCHH <- getData(aCHH.3x.100.50p.all)
bCHH <- getData(bCHH.3x.100.50p.all)
cCHH <- getData(cCHH.3x.100.50p.all)

WTvMutCHH <- data.frame("model"=c("full model", "A", "B", "C"),
                        "totaltiles"=c(diff.CHH.3x.100.oChi@num.records, diff.aCHH.3x.100@num.records, diff.bCHH.3x.100@num.records, diff.cCHH.3x.100@num.records),
                        "DMtiles"=c(CHH.3x.100.oChi.50p.all@num.records, aCHH.3x.100.50p.all@num.records, bCHH.3x.100.50p.all@num.records, cCHH.3x.100.50p.all@num.records),
                        "hypo"=c(nrow(CHH[CHH$meth.diff<0,]), nrow(aCHH[aCHH$meth.diff<0,]),nrow(bCHH[bCHH$meth.diff<0,]), nrow(cCHH[cCHH$meth.diff<0,])),
                        "hyper"=c(nrow(CHH[CHH$meth.diff>0,]), nrow(aCHH[aCHH$meth.diff>0,]),nrow(bCHH[bCHH$meth.diff>0,]), nrow(cCHH[cCHH$meth.diff>0,])))

WTvMutCHH$percent <- round(WTvMutCHH$DMtiles/WTvMutCHH$totaltiles*100,3)
WTvMutCHH$hypoPerc <- round(WTvMutCHH$hypo/WTvMutCHH$totaltiles*100,3)
WTvMutCHH$hyperPerc <- round(WTvMutCHH$hyper/WTvMutCHH$totaltiles*100,3)
print(WTvMutCHH)

# plot percent of all tiles considered DM
WTvMutCHH$context="CHH"
WTvMutCHG$context="CHG"
WTvMutCpG$context="CpG"
new <- rbind(WTvMutCHG, WTvMutCHH, WTvMutCpG)
new[2:6] <- list(NULL)
new2 <- melt(new)
new2$model_f <- factor(new2$model,levels = c("full model", "A", "B", "C"))
new2$variable_f <- factor(new2$variable, levels=c("hyperPerc", "hypoPerc"))
new2$context_f <- factor(new2$context, levels=c("CpG", "CHG", "CHH"))
new2 <- subset(new2, model!= "full model")

perc <- ggplot(data=new2, aes(x=context_f, y=value, fill=variable_f)) + 
  geom_bar(stat="identity") + 
  facet_grid(~model_f) +
  scale_fill_lancet(name="Type of DM", labels=c("Hypermethylated","Hypomethylated")) +
  labs(title="Percentage of tiles with methylation difference +/- 50% and q-value < 0.01", 
       x ="Cytosine sequence context", y = "% of tiles considered significantly DM") 

ggsave(filename = "nofull-percentDM.tiff", device="tiff", path="/home/nia/project/nia/results/3.3-quantification/", plot=perc)

# plot size of raw files, average coverage per sample
df <- data.frame("geno"=c("A", "A", "B", "B", "C", "C"),
                 "status"=c("Wild-type", "Mutant","Wild-type", "Mutant","Wild-type", "Mutant"),
                 "pairs"=c(195494645, 173618601, 197454925, 50218874, 187434448,152910311),
                 "name"=c("A wild-type", "A Mutant", "B wild-type", "B Mutant", "C wild-type", "C Mutant"),
                 "coverage"=c(7.1916, 6.6481, 6.9162, 2.2664, 6.3997, 5.909))

df$`mop1 status` <- factor(df$status, levels=c("Wild-type", "Mutant"))
# size of files in # of paired-end reads
ggplot(df, aes(fill=`mop1 status`, y=pairs, x=geno)) + 
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_lancet() + 
  scale_y_continuous(labels = scales::comma) + 
  labs(title="Number of paired-end reads per sample",
       x="Genotype",
       y="Number of read pairs")
# average coverage
ggplot(df, aes(fill=`mop1 status`, y=coverage, x=geno)) + 
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_lancet() + 
  scale_y_continuous(labels = scales::comma) + 
  labs(title="Mean genomic coverage",
       x="Genotype",
       y="Mean coverage")
