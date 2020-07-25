###################################################################
### VOLCANO PLOTS OF TILES BY SIGNIFICANCE AND METHYLATION DIFF
### R script
###################################################################

# set up environment
setwd("/scratch/nia/manuFinal")
load(file=".RData")
library(methylKit)
library(GenomicRanges)
library(genomation)
library(ggplot2)
library(dplyr)
library(cowplot)
library(wesanderson)
library(RColorBrewer)
library(ggsci)

# Overdispersion correction -----------------------------------------------
# CpG
CpG.nocorr.10 <- getData(diff.CpG.3x.100.nocorr)[getData(diff.CpG.3x.100.nocorr)$chr=="10",]
CpG.nocorr.10 <- CpG.nocorr.10 %>% mutate(threshold=ifelse(meth.diff <= -50 & qvalue <= 0.01 , "A",
                                                           ifelse(meth.diff >= 50 & qvalue <= 0.01 , "B",
                                                                  ifelse(abs(meth.diff) >= 50, "C",
                                                                         ifelse(qvalue <= 0.01,"D","E")))))
CpG.oChi.10 <- getData(diff.CpG.3x.100.oChi)[getData(diff.CpG.3x.100.oChi)$chr=="10",]
CpG.oChi.10 <- CpG.oChi.10 %>% mutate(threshold=ifelse(meth.diff <= -50 & qvalue <= 0.01 , "A",
                                                       ifelse(meth.diff >= 50 & qvalue <= 0.01 , "B",
                                                              ifelse(abs(meth.diff) >= 50, "C",
                                                                     ifelse(qvalue <= 0.01,"D","E")))))
CpG.oF.10 <- getData(diff.CpG.3x.100.oF)[getData(diff.CpG.3x.100.oF)$chr=="10",]
CpG.oF.10 <- CpG.oF.10 %>% mutate(threshold=ifelse(meth.diff <= -50 & qvalue <= 0.01 , "A",
                                                   ifelse(meth.diff >= 50 & qvalue <= 0.01 , "B",
                                                          ifelse(abs(meth.diff) >= 50, "C",
                                                                 ifelse(qvalue <= 0.01,"D","E")))))
CpG.numrec <- data.frame(Correction=c("Uncorrected", "Chi-squared", "F-test"),
                         Tiles=c(CpG.3x.100.nocorr.50p.all@num.records, CpG.3x.100.oChi.50p.all@num.records, CpG.3x.100.oF.50p.all@num.records))

pl.CpG.nocorr = ggplot(data=CpG.nocorr.10, aes(x=meth.diff, y=-log10(qvalue))) +
  geom_point(alpha=0.5, size=1,aes(colour=threshold)) + 
  theme(legend.position="none") +
  ylim(0,-log10(min(CpG.nocorr.10$qvalue))) +
  ggtitle("Uncorrected") +
  xlab("Percent methylation difference") + ylab("-log10(qvalue)") +
  scale_colour_manual(values=c(A="red", B="green",C="gray59",D="gray59",E="gray59")) +
  geom_vline(xintercept = 50, linetype="dashed", color = "blue", size=.25) +
  geom_vline(xintercept = -50, linetype="dashed", color = "blue", size=.25) +
  geom_hline(yintercept = -log10(0.01), linetype="dashed", color = "blue", size=.25)


pl.CpG.oChi = ggplot(data=CpG.oChi.10, aes(x=meth.diff, y=-log10(qvalue))) +
  geom_point(alpha=0.5, size=1,aes(colour=threshold)) + 
  theme(legend.position="none") +
  ylim(0,-log10(min(CpG.nocorr.10$qvalue))) +
  ggtitle("Chi-squared test") +
  xlab("Percent methylation difference") + ylab("-log10(qvalue)") +
  scale_colour_manual(values=c(A="red", B="green",C="gray59",D="gray59",E="gray59")) +
  geom_vline(xintercept = 50, linetype="dashed", color = "blue", size=.25) +
  geom_vline(xintercept = -50, linetype="dashed", color = "blue", size=.25) +
  geom_hline(yintercept = -log10(0.01), linetype="dashed", color = "blue", size=.25)


pl.CpG.oF = ggplot(data=CpG.oF.10, aes(x=meth.diff, y=-log10(qvalue))) +
  geom_point(alpha=0.5, size=1,aes(colour=threshold)) + 
  theme(legend.position="none") +
  ylim(0,-log10(min(CpG.nocorr.10$qvalue))) +
  ggtitle("F-test") +
  xlab("Percent methylation difference") + ylab("-log10(qvalue)") +
  scale_colour_manual(values=c(A="red", B="green",C="gray59",D="gray59",E="gray59")) +
  geom_vline(xintercept = 50, linetype="dashed", color = "blue", size=.25) +
  geom_vline(xintercept = -50, linetype="dashed", color = "blue", size=.25) +
  geom_hline(yintercept = -log10(0.01), linetype="dashed", color = "blue", size=.25)

pl.CpG.numrec <- ggplot(data=CpG.numrec, aes(x=Correction,y=Tiles,fill=Correction)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=Tiles), vjust=2, size=3.5) +
  ggtitle("Number of DM tiles") +
  theme(legend.position="none") +
  scale_fill_npg()
  

pl.CpG.all <- plot_grid(pl.CpG.nocorr, pl.CpG.oChi, pl.CpG.oF, pl.CpG.numrec, nrow=2, labels="AUTO")
save_plot("/home/nia/project/nia/results/3.4-volcanoPlots/volcanoCpG.png", pl.CpG.all, ncol=2, nrow=2)

# CHG
CHG.nocorr.10 <- getData(diff.CHG.3x.100.nocorr)[getData(diff.CHG.3x.100.nocorr)$chr=="10",]
CHG.nocorr.10 <- CHG.nocorr.10 %>% mutate(threshold=ifelse(meth.diff <= -50 & qvalue <= 0.01 , "A",
                                                           ifelse(meth.diff >= 50 & qvalue <= 0.01 , "B",
                                                                  ifelse(abs(meth.diff) >= 50, "C",
                                                                         ifelse(qvalue <= 0.01,"D","E")))))
CHG.oChi.10 <- getData(diff.CHG.3x.100.oChi)[getData(diff.CHG.3x.100.oChi)$chr=="10",]
CHG.oChi.10 <- CHG.oChi.10 %>% mutate(threshold=ifelse(meth.diff <= -50 & qvalue <= 0.01 , "A",
                                                       ifelse(meth.diff >= 50 & qvalue <= 0.01 , "B",
                                                              ifelse(abs(meth.diff) >= 50, "C",
                                                                     ifelse(qvalue <= 0.01,"D","E")))))
CHG.oF.10 <- getData(diff.CHG.3x.100.oF)[getData(diff.CHG.3x.100.oF)$chr=="10",]
CHG.oF.10 <- CHG.oF.10 %>% mutate(threshold=ifelse(meth.diff <= -50 & qvalue <= 0.01 , "A",
                                                   ifelse(meth.diff >= 50 & qvalue <= 0.01 , "B",
                                                          ifelse(abs(meth.diff) >= 50, "C",
                                                                 ifelse(qvalue <= 0.01,"D","E")))))
CHG.numrec <- data.frame(Correction=c("Uncorrected", "Chi-squared", "F-test"),
                         Tiles=c(CHG.3x.100.nocorr.50p.all@num.records, CHG.3x.100.oChi.50p.all@num.records, CHG.3x.100.oF.50p.all@num.records))

pl.CHG.nocorr = ggplot(data=CHG.nocorr.10, aes(x=meth.diff, y=-log10(qvalue))) +
  geom_point(alpha=0.5, size=1,aes(colour=threshold)) + 
  theme(legend.position="none") +
  ylim(0,-log10(min(CHG.nocorr.10$qvalue))) +
  ggtitle("Uncorrected") +
  xlab("Percent methylation difference") + ylab("-log10(qvalue)") +
  scale_colour_manual(values=c(A="red", B="green",C="gray59",D="gray59",E="gray59")) +
  geom_vline(xintercept = 50, linetype="dashed", color = "blue", size=.25) +
  geom_vline(xintercept = -50, linetype="dashed", color = "blue", size=.25) +
  geom_hline(yintercept = -log10(0.01), linetype="dashed", color = "blue", size=.25)


pl.CHG.oChi = ggplot(data=CHG.oChi.10, aes(x=meth.diff, y=-log10(qvalue))) +
  geom_point(alpha=0.5, size=1,aes(colour=threshold)) + 
  theme(legend.position="none") +
  ylim(0,-log10(min(CHG.nocorr.10$qvalue))) +
  ggtitle("Chi-squared test") +
  xlab("Percent methylation difference") + ylab("-log10(qvalue)") +
  scale_colour_manual(values=c(A="red", B="green",C="gray59",D="gray59",E="gray59")) +
  geom_vline(xintercept = 50, linetype="dashed", color = "blue", size=.25) +
  geom_vline(xintercept = -50, linetype="dashed", color = "blue", size=.25) +
  geom_hline(yintercept = -log10(0.01), linetype="dashed", color = "blue", size=.25)


pl.CHG.oF = ggplot(data=CHG.oF.10, aes(x=meth.diff, y=-log10(qvalue))) +
  geom_point(alpha=0.5, size=1,aes(colour=threshold)) + 
  theme(legend.position="none") +
  ylim(0,-log10(min(CHG.nocorr.10$qvalue))) +
  ggtitle("F-test") +
  xlab("Percent methylation difference") + ylab("-log10(qvalue)") +
  scale_colour_manual(values=c(A="red", B="green",C="gray59",D="gray59",E="gray59")) +
  geom_vline(xintercept = 50, linetype="dashed", color = "blue", size=.25) +
  geom_vline(xintercept = -50, linetype="dashed", color = "blue", size=.25) +
  geom_hline(yintercept = -log10(0.01), linetype="dashed", color = "blue", size=.25)

pl.CHG.numrec <- ggplot(data=CHG.numrec, aes(x=Correction,y=Tiles,fill=Correction)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=Tiles), vjust=2, size=3.5) +
  ggtitle("Number of DM tiles") +
  theme(legend.position="none") +
  scale_fill_npg()


pl.CHG.all <- plot_grid(pl.CHG.nocorr, pl.CHG.oChi, pl.CHG.oF, pl.CHG.numrec, nrow=2, labels="AUTO")
save_plot("/home/nia/project/nia/results/3.4-volcanoPlots/volcanoCHG.png", pl.CHG.all, ncol=2, nrow=2)

# CHH
CHH.nocorr.10 <- getData(diff.CHH.3x.100.nocorr)[getData(diff.CHH.3x.100.nocorr)$chr=="10",]
CHH.nocorr.10 <- CHH.nocorr.10 %>% mutate(threshold=ifelse(meth.diff <= -50 & qvalue <= 0.01 , "A",
                                                           ifelse(meth.diff >= 50 & qvalue <= 0.01 , "B",
                                                                  ifelse(abs(meth.diff) >= 50, "C",
                                                                         ifelse(qvalue <= 0.01,"D","E")))))
CHH.oChi.10 <- getData(diff.CHH.3x.100.oChi)[getData(diff.CHH.3x.100.oChi)$chr=="10",]
CHH.oChi.10 <- CHH.oChi.10 %>% mutate(threshold=ifelse(meth.diff <= -50 & qvalue <= 0.01 , "A",
                                                       ifelse(meth.diff >= 50 & qvalue <= 0.01 , "B",
                                                              ifelse(abs(meth.diff) >= 50, "C",
                                                                     ifelse(qvalue <= 0.01,"D","E")))))
CHH.oF.10 <- getData(diff.CHH.3x.100.oF)[getData(diff.CHH.3x.100.oF)$chr=="10",]
CHH.oF.10 <- CHH.oF.10 %>% mutate(threshold=ifelse(meth.diff <= -50 & qvalue <= 0.01 , "A",
                                                   ifelse(meth.diff >= 50 & qvalue <= 0.01 , "B",
                                                          ifelse(abs(meth.diff) >= 50, "C",
                                                                 ifelse(qvalue <= 0.01,"D","E")))))
CHH.numrec <- data.frame(Correction=c("Uncorrected", "Chi-squared", "F-test"),
                         Tiles=c(CHH.3x.100.nocorr.50p.all@num.records, CHH.3x.100.oChi.50p.all@num.records, CHH.3x.100.oF.50p.all@num.records))

pl.CHH.nocorr = ggplot(data=CHH.nocorr.10, aes(x=meth.diff, y=-log10(qvalue))) +
  geom_point(alpha=0.5, size=1,aes(colour=threshold)) + 
  theme(legend.position="none") +
  ylim(0,-log10(min(CHH.nocorr.10$qvalue))) +
  ggtitle("Uncorrected") +
  xlab("Percent methylation difference") + ylab("-log10(qvalue)") +
  scale_colour_manual(values=c(A="red", B="green",C="gray59",D="gray59",E="gray59")) +
  geom_vline(xintercept = 50, linetype="dashed", color = "blue", size=.25) +
  geom_vline(xintercept = -50, linetype="dashed", color = "blue", size=.25) +
  geom_hline(yintercept = -log10(0.01), linetype="dashed", color = "blue", size=.25)


pl.CHH.oChi = ggplot(data=CHH.oChi.10, aes(x=meth.diff, y=-log10(qvalue))) +
  geom_point(alpha=0.5, size=1,aes(colour=threshold)) + 
  theme(legend.position="none") +
  ylim(0,-log10(min(CHH.nocorr.10$qvalue))) +
  ggtitle("Chi-squared test") +
  xlab("Percent methylation difference") + ylab("-log10(qvalue)") +
  scale_colour_manual(values=c(A="red", B="green",C="gray59",D="gray59",E="gray59")) +
  geom_vline(xintercept = 50, linetype="dashed", color = "blue", size=.25) +
  geom_vline(xintercept = -50, linetype="dashed", color = "blue", size=.25) +
  geom_hline(yintercept = -log10(0.01), linetype="dashed", color = "blue", size=.25)


pl.CHH.oF = ggplot(data=CHH.oF.10, aes(x=meth.diff, y=-log10(qvalue))) +
  geom_point(alpha=0.5, size=1,aes(colour=threshold)) + 
  theme(legend.position="none") +
  ylim(0,-log10(min(CHH.nocorr.10$qvalue))) +
  ggtitle("F-test") +
  xlab("Percent methylation difference") + ylab("-log10(qvalue)") +
  scale_colour_manual(values=c(A="red", B="green",C="gray59",D="gray59",E="gray59")) +
  geom_vline(xintercept = 50, linetype="dashed", color = "blue", size=.25) +
  geom_vline(xintercept = -50, linetype="dashed", color = "blue", size=.25) +
  geom_hline(yintercept = -log10(0.01), linetype="dashed", color = "blue", size=.25)

pl.CHH.numrec <- ggplot(data=CHH.numrec, aes(x=Correction,y=Tiles,fill=Correction)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=Tiles), vjust=2, size=3.5) +
  ggtitle("Number of DM tiles") +
  theme(legend.position="none") +
  scale_fill_npg()


pl.CHH.all <- plot_grid(pl.CHH.nocorr, pl.CHH.oChi, pl.CHH.oF, pl.CHH.numrec, nrow=2, labels="AUTO")
save_plot("/home/nia/project/nia/results/3.4-volcanoPlots/volcanoCHH.png", pl.CHH.all, ncol=2, nrow=2)



# Per chromosome, all comps -----------------------------------------------
input <- list(diff.CpG.3x.100.nocorr, diff.CHG.3x.100.nocorr, diff.CHH.3x.100.nocorr,
              diff.CpG.3x.100.oChi, diff.CHG.3x.100.oChi, diff.CHH.3x.100.oChi,
              diff.CpG.3x.100.oF, diff.CHG.3x.100.oF, diff.CHH.3x.100.oF,
              diff.aCpG.3x.100, diff.aCHG.3x.100, diff.aCHH.3x.100,
              diff.bCpG.3x.100, diff.bCHG.3x.100, diff.bCHH.3x.100,
              diff.cCpG.3x.100, diff.cCHG.3x.100, diff.cCHH.3x.100)
names(input) <- c("full.CpG.nocorr","full.CHG.nocorr","full.CHH.nocorr",
                  "full.CpG.oChi","full.CHG.oChi","full.CHH.oChi",
                  "full.CpG.oF","full.CHG.oF","full.CHH.oF",
                  "aCpG", "aCHG", "aCHH",
                  "bCpG", "bCHG", "bCHH",
                  "cCpG", "cCHG", "cCHH")
titles <- list("Full model, no correction", "Full model, no correction", "Full model, no correction",
               "Full model, Chi-squared", "Full model, Chi-squared", "Full model, Chi-squared",
               "Full model, F-test", "Full model, F-test", "Full model, F-test",
               "A","A","A",
               "B","B","B",
               "C","C","C")

outV <- Map(function(x, title) {
  y <- getData(x)
  z <- y %>% mutate(threshold=ifelse(meth.diff <= -50 & qvalue <= 0.01 , "A",
                                    ifelse(meth.diff >= 50 & qvalue <= 0.01 , "B",
                                           ifelse(abs(meth.diff) >= 50, "C",
                                                  ifelse(qvalue <= 0.01,"D","E")))))
  a <- split(z, z$chr)
  
  chrV <- vector(mode = "list", length = 10)
  
  for (i in 1:10) {
    chrV[[i]] <- ggplot(data=a[[i]], aes(x=meth.diff, y=-log10(qvalue))) +
      geom_point(alpha=0.5, size=1,aes(colour=threshold)) + 
      theme(legend.position="none") +
      xlim(-100,100) +
      ggtitle(paste(title)) +
      xlab("Percent methylation difference") + ylab("-log10(qvalue)") +
      scale_colour_manual(values=c(A="red", B="green",C="gray59",D="gray59",E="gray59")) +
      geom_vline(xintercept = 50, linetype="dashed", color = "blue", size=.25) +
      geom_vline(xintercept = -50, linetype="dashed", color = "blue", size=.25) +
      geom_hline(yintercept = -log10(0.01), linetype="dashed", color = "blue", size=.25)
  }
  
  return(chrV)
  
}, input, titles)


# CpG
for (i in 1:10) {
  yVal <- max(layer_scales(outV$full.CpG.nocorr[[i]])$y$range$range,layer_scales(outV$full.CpG.oChi[[i]])$y$range$range,layer_scales(outV$full.CpG.oF[[i]])$y$range$range,
            layer_scales(outV$aCpG[[i]])$y$range$range,layer_scales(outV$bCpG[[i]])$y$range$range,layer_scales(outV$cCpG[[i]])$y$range$range)
  title.plot <- ggdraw() + 
    draw_label(paste("Chromosome", i, "CpG"), fontface = 'bold', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  plot1 <- plot_grid(outV$full.CpG.nocorr[[i]] + ylim(0,yVal), outV$full.CpG.oChi[[i]] + ylim(0,yVal), outV$full.CpG.oF[[i]] + ylim(0,yVal),
                     outV$aCpG[[i]] + ylim(0,yVal), outV$bCpG[[i]] + ylim(0,yVal), outV$cCpG[[i]] + ylim(0,yVal), 
                     labels="AUTO", ncol=2, nrow=3)
  plot2 <- plot_grid(title.plot, plot1, labels=NULL, rel_heights = c(0.1, 1), ncol=1)
  save_plot(paste0("/home/nia/project/nia/results/3.4-volcanoPlots/CpG-chr", i,"-volcano.png"),plot2, ncol=2, nrow=3)
}

# CHG
for (i in 1:10) {
  yVal <- max(layer_scales(outV$full.CHG.nocorr[[i]])$y$range$range,layer_scales(outV$full.CHG.oChi[[i]])$y$range$range,layer_scales(outV$full.CHG.oF[[i]])$y$range$range,
              layer_scales(outV$aCHG[[i]])$y$range$range,layer_scales(outV$bCHG[[i]])$y$range$range,layer_scales(outV$cCHG[[i]])$y$range$range)
  title.plot <- ggdraw() + 
    draw_label(paste("Chromosome", i, "CHG"), fontface = 'bold', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  plot1 <- plot_grid(outV$full.CHG.nocorr[[i]] + ylim(0,yVal), outV$full.CHG.oChi[[i]] + ylim(0,yVal), outV$full.CHG.oF[[i]] + ylim(0,yVal),
                     outV$aCHG[[i]] + ylim(0,yVal), outV$bCHG[[i]] + ylim(0,yVal), outV$cCHG[[i]] + ylim(0,yVal), 
                     labels="AUTO", ncol=2, nrow=3)
  plot2 <- plot_grid(title.plot, plot1, labels=NULL, rel_heights = c(0.1, 1), ncol=1)
  save_plot(paste0("/home/nia/project/nia/results/3.4-volcanoPlots/CHG-chr", i,"-volcano.png"),plot2, ncol=2, nrow=3)
}

# CHH
for (i in 1:10) {
  yVal <- max(layer_scales(outV$full.CHH.nocorr[[i]])$y$range$range,layer_scales(outV$full.CHH.oChi[[i]])$y$range$range,layer_scales(outV$full.CHH.oF[[i]])$y$range$range,
              layer_scales(outV$aCHH[[i]])$y$range$range,layer_scales(outV$bCHH[[i]])$y$range$range,layer_scales(outV$cCHH[[i]])$y$range$range)
  title.plot <- ggdraw() +
    draw_label(paste("Chromosome", i, "CHH"), fontface = 'bold', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  plot1 <- plot_grid(outV$full.CHH.nocorr[[i]] + ylim(0,yVal), outV$full.CHH.oChi[[i]] + ylim(0,yVal), outV$full.CHH.oF[[i]] + ylim(0,yVal),
                     outV$aCHH[[i]] + ylim(0,yVal), outV$bCHH[[i]] + ylim(0,yVal), outV$cCHH[[i]] + ylim(0,yVal),
                     labels="AUTO", ncol=2, nrow=3)
  plot2 <- plot_grid(title.plot, plot1, labels=NULL, rel_heights = c(0.1, 1), ncol=1)
  save_plot(paste0("/home/nia/project/nia/results/3.4-volcanoPlots/CHH-chr", i,"-volcano.png"),plot2, ncol=2, nrow=3)
}
