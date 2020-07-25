###################################################################
### PREPARING ANNOTATED TILES FOR DOWNSTREAM ANALYSIS
### Local R script 
###################################################################

# set up environment
setwd("D:/Nia/OneDrive - University of Guelph/Working/DMRanalysis")
library(dplyr)

# CpG ---------------------------------------------------------------------
# Load in files
complete.CpG<-read.csv(file="gene-full.CpG.csv", sep=",")
complete.aCpG<-read.csv(file="gene-A.CpG.csv", sep=",")
complete.bCpG<-read.csv(file="gene-B.CpG.csv", sep=",")
complete.cCpG<-read.csv(file="gene-C.CpG.csv", sep=",")

## Split on feature
tile.CpG<-split(complete.CpG, complete.CpG$feature)
names(tile.CpG) <- c("exons", "introns", "promoters")
tile.aCpG<-split(complete.aCpG, complete.aCpG$feature)
names(tile.aCpG) <- c("exons", "introns", "promoters")
tile.bCpG<-split(complete.bCpG, complete.bCpG$feature)
names(tile.bCpG) <- c("exons", "introns", "promoters")
tile.cCpG<-split(complete.cCpG, complete.cCpG$feature)
names(tile.cCpG) <- c("exons", "introns", "promoters")

## Split on promoter/genic
split1.CpG<-split(complete.CpG, !complete.CpG$feature=="promoters")
names(split1.CpG) <- c("promoters", "genic")
split1.aCpG<-split(complete.aCpG, !complete.aCpG$feature=="promoters")
names(split1.aCpG) <- c("promoters", "genic")
split1.bCpG<-split(complete.bCpG, !complete.bCpG$feature=="promoters")
names(split1.bCpG) <- c("promoters", "genic")
split1.cCpG<-split(complete.cCpG, !complete.cCpG$feature=="promoters")
names(split1.cCpG) <- c("promoters", "genic")
#Get unique tiles
feat.CpG <- lapply(split1.CpG, function(x) x[!duplicated(x[c("chr","start")]), ])
feat.aCpG <- lapply(split1.aCpG, function(x) x[!duplicated(x[c("chr","start")]), ])
feat.bCpG <- lapply(split1.bCpG, function(x) x[!duplicated(x[c("chr","start")]), ])
feat.cCpG <- lapply(split1.cCpG, function(x) x[!duplicated(x[c("chr","start")]), ])

## Split on feature type, only unique tiles (no tile that overlaps any features included)
u.CpG<-complete.CpG[!(duplicated(complete.CpG[,1:2]) | duplicated(complete.CpG[,1:2], fromLast = TRUE)), ]
u.aCpG<-complete.aCpG[!(duplicated(complete.aCpG[,1:2]) | duplicated(complete.aCpG[,1:2], fromLast = TRUE)), ]
u.bCpG<-complete.bCpG[!(duplicated(complete.bCpG[,1:2]) | duplicated(complete.bCpG[,1:2], fromLast = TRUE)), ]
u.cCpG<-complete.cCpG[!(duplicated(complete.cCpG[,1:2]) | duplicated(complete.cCpG[,1:2], fromLast = TRUE)), ]
unique.CpG<-split(u.CpG, interaction(!u.CpG$feature=="promoters",u.CpG$meth.diff<0))
names(unique.CpG) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
unique.aCpG<-split(u.aCpG, interaction(!u.aCpG$feature=="promoters",u.aCpG$meth.diff<0))
names(unique.aCpG) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
unique.bCpG<-split(u.bCpG, interaction(!u.bCpG$feature=="promoters",u.bCpG$meth.diff<0))
names(unique.bCpG) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
unique.cCpG<-split(u.cCpG, interaction(!u.cCpG$feature=="promoters",u.cCpG$meth.diff<0))
names(unique.cCpG) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")

## Split on methylation >0<
split2.CpG<-split(complete.CpG, complete.CpG$meth.diff<0)
names(split2.CpG) <- c("hyper", "hypo")
split2.aCpG<-split(complete.aCpG, complete.aCpG$meth.diff<0)
names(split2.aCpG) <- c("hyper", "hypo")
split2.bCpG<-split(complete.bCpG, complete.bCpG$meth.diff<0)
names(split2.bCpG) <- c("hyper", "hypo")
split2.cCpG<-split(complete.cCpG, complete.cCpG$meth.diff<0)
names(split2.cCpG) <- c("hyper", "hypo")
#Retrieve unique names (for easier access in later steps)
hyp.CpG <- lapply(split2.CpG, function(x) x[!duplicated(x[c("name")]), 'name'])
hyp.aCpG <- lapply(split2.aCpG, function(x) x[!duplicated(x[c("name")]), 'name'])
hyp.bCpG <- lapply(split2.bCpG, function(x) x[!duplicated(x[c("name")]), 'name'])
hyp.cCpG <- lapply(split2.cCpG, function(x) x[!duplicated(x[c("name")]), 'name'])

## Frequency of tile:gene mappings
freq.CpG <- complete.CpG[!duplicated(complete.CpG[c("chr","start","name")]),] %>% 
  group_by(name) %>% summarise(counts=n())
freq.aCpG <- complete.aCpG[!duplicated(complete.aCpG[c("chr","start","name")]),] %>% 
  group_by(name) %>% summarise(counts=n())
freq.bCpG <- complete.bCpG[!duplicated(complete.bCpG[c("chr","start","name")]),] %>% 
  group_by(name) %>% summarise(counts=n())
freq.cCpG <- complete.cCpG[!duplicated(complete.cCpG[c("chr","start","name")]),] %>% 
  group_by(name) %>% summarise(counts=n())

## Split on both categories
split3.CpG<-split(complete.CpG, interaction(!complete.CpG$feature=="promoters",complete.CpG$meth.diff<0))
names(split3.CpG) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
split3.aCpG<-split(complete.aCpG, interaction(!complete.aCpG$feature=="promoters",complete.aCpG$meth.diff<0))
names(split3.aCpG) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
split3.bCpG<-split(complete.bCpG, interaction(!complete.bCpG$feature=="promoters",complete.bCpG$meth.diff<0))
names(split3.bCpG) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
split3.cCpG<-split(complete.cCpG, interaction(!complete.cCpG$feature=="promoters",complete.cCpG$meth.diff<0))
names(split3.cCpG) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
#Retrieve unique names
gene.CpG <- lapply(split3.CpG, function(x) x[!duplicated(x[c("name")]), 'name'])
gene.aCpG <- lapply(split3.aCpG, function(x) x[!duplicated(x[c("name")]), 'name'])
gene.bCpG <- lapply(split3.bCpG, function(x) x[!duplicated(x[c("name")]), 'name'])
gene.cCpG <- lapply(split3.cCpG, function(x) x[!duplicated(x[c("name")]), 'name'])

# CHG ---------------------------------------------------------------------
# Load in files
complete.CHG<-read.csv(file="gene-full.CHG.csv", sep=",")
complete.aCHG<-read.csv(file="gene-A.CHG.csv", sep=",")
complete.bCHG<-read.csv(file="gene-B.CHG.csv", sep=",")
complete.cCHG<-read.csv(file="gene-C.CHG.csv", sep=",")

## Split on feature
tile.CHG<-split(complete.CHG, complete.CHG$feature)
names(tile.CHG) <- c("exons", "introns", "promoters")
tile.aCHG<-split(complete.aCHG, complete.aCHG$feature)
names(tile.aCHG) <- c("exons", "introns", "promoters")
tile.bCHG<-split(complete.bCHG, complete.bCHG$feature)
names(tile.bCHG) <- c("exons", "introns", "promoters")
tile.cCHG<-split(complete.cCHG, complete.cCHG$feature)
names(tile.cCHG) <- c("exons", "introns", "promoters")

## Split on promoter/genic
split1.CHG<-split(complete.CHG, !complete.CHG$feature=="promoters")
names(split1.CHG) <- c("promoters", "genic")
split1.aCHG<-split(complete.aCHG, !complete.aCHG$feature=="promoters")
names(split1.aCHG) <- c("promoters", "genic")
split1.bCHG<-split(complete.bCHG, !complete.bCHG$feature=="promoters")
names(split1.bCHG) <- c("promoters", "genic")
split1.cCHG<-split(complete.cCHG, !complete.cCHG$feature=="promoters")
names(split1.cCHG) <- c("promoters", "genic")
#Get unique tiles
feat.CHG <- lapply(split1.CHG, function(x) x[!duplicated(x[c("chr","start")]), ])
feat.aCHG <- lapply(split1.aCHG, function(x) x[!duplicated(x[c("chr","start")]), ])
feat.bCHG <- lapply(split1.bCHG, function(x) x[!duplicated(x[c("chr","start")]), ])
feat.cCHG <- lapply(split1.cCHG, function(x) x[!duplicated(x[c("chr","start")]), ])

## Split on feature type, only unique tiles (no tile that overlaps any features included)
u.CHG<-complete.CHG[!(duplicated(complete.CHG[,1:2]) | duplicated(complete.CHG[,1:2], fromLast = TRUE)), ]
u.aCHG<-complete.aCHG[!(duplicated(complete.aCHG[,1:2]) | duplicated(complete.aCHG[,1:2], fromLast = TRUE)), ]
u.bCHG<-complete.bCHG[!(duplicated(complete.bCHG[,1:2]) | duplicated(complete.bCHG[,1:2], fromLast = TRUE)), ]
u.cCHG<-complete.cCHG[!(duplicated(complete.cCHG[,1:2]) | duplicated(complete.cCHG[,1:2], fromLast = TRUE)), ]
unique.CHG<-split(u.CHG, interaction(!u.CHG$feature=="promoters",u.CHG$meth.diff<0))
names(unique.CHG) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
unique.aCHG<-split(u.aCHG, interaction(!u.aCHG$feature=="promoters",u.aCHG$meth.diff<0))
names(unique.aCHG) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
unique.bCHG<-split(u.bCHG, interaction(!u.bCHG$feature=="promoters",u.bCHG$meth.diff<0))
names(unique.bCHG) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
unique.cCHG<-split(u.cCHG, interaction(!u.cCHG$feature=="promoters",u.cCHG$meth.diff<0))
names(unique.cCHG) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")

## Split on methylation >0<
split2.CHG<-split(complete.CHG, complete.CHG$meth.diff<0)
names(split2.CHG) <- c("hyper", "hypo")
split2.aCHG<-split(complete.aCHG, complete.aCHG$meth.diff<0)
names(split2.aCHG) <- c("hyper", "hypo")
split2.bCHG<-split(complete.bCHG, complete.bCHG$meth.diff<0)
names(split2.bCHG) <- c("hyper", "hypo")
split2.cCHG<-split(complete.cCHG, complete.cCHG$meth.diff<0)
names(split2.cCHG) <- c("hyper", "hypo")
#Retrieve unique names (for easier access in later steps)
hyp.CHG <- lapply(split2.CHG, function(x) x[!duplicated(x[c("name")]), 'name'])
hyp.aCHG <- lapply(split2.aCHG, function(x) x[!duplicated(x[c("name")]), 'name'])
hyp.bCHG <- lapply(split2.bCHG, function(x) x[!duplicated(x[c("name")]), 'name'])
hyp.cCHG <- lapply(split2.cCHG, function(x) x[!duplicated(x[c("name")]), 'name'])

## Frequency of tile:gene mappings
freq.CHG <- complete.CHG[!duplicated(complete.CHG[c("chr","start","name")]),] %>% 
  group_by(name) %>% summarise(counts=n())
freq.aCHG <- complete.aCHG[!duplicated(complete.aCHG[c("chr","start","name")]),] %>% 
  group_by(name) %>% summarise(counts=n())
freq.bCHG <- complete.bCHG[!duplicated(complete.bCHG[c("chr","start","name")]),] %>% 
  group_by(name) %>% summarise(counts=n())
freq.cCHG <- complete.cCHG[!duplicated(complete.cCHG[c("chr","start","name")]),] %>% 
  group_by(name) %>% summarise(counts=n())

## Split on both categories
split3.CHG<-split(complete.CHG, interaction(!complete.CHG$feature=="promoters",complete.CHG$meth.diff<0))
names(split3.CHG) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
split3.aCHG<-split(complete.aCHG, interaction(!complete.aCHG$feature=="promoters",complete.aCHG$meth.diff<0))
names(split3.aCHG) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
split3.bCHG<-split(complete.bCHG, interaction(!complete.bCHG$feature=="promoters",complete.bCHG$meth.diff<0))
names(split3.bCHG) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
split3.cCHG<-split(complete.cCHG, interaction(!complete.cCHG$feature=="promoters",complete.cCHG$meth.diff<0))
names(split3.cCHG) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
#Retrieve unique names
gene.CHG <- lapply(split3.CHG, function(x) x[!duplicated(x[c("name")]), 'name'])
gene.aCHG <- lapply(split3.aCHG, function(x) x[!duplicated(x[c("name")]), 'name'])
gene.bCHG <- lapply(split3.bCHG, function(x) x[!duplicated(x[c("name")]), 'name'])
gene.cCHG <- lapply(split3.cCHG, function(x) x[!duplicated(x[c("name")]), 'name'])

# CHH ---------------------------------------------------------------------
# Load in files
complete.CHH<-read.csv(file="gene-full.CHH.csv", sep=",")
complete.aCHH<-read.csv(file="gene-A.CHH.csv", sep=",")
complete.bCHH<-read.csv(file="gene-B.CHH.csv", sep=",")
complete.cCHH<-read.csv(file="gene-C.CHH.csv", sep=",")

## Split on feature
tile.CHH<-split(complete.CHH, complete.CHH$feature)
names(tile.CHH) <- c("exons", "introns", "promoters")
tile.aCHH<-split(complete.aCHH, complete.aCHH$feature)
names(tile.aCHH) <- c("exons", "introns", "promoters")
tile.bCHH<-split(complete.bCHH, complete.bCHH$feature)
names(tile.bCHH) <- c("exons", "introns", "promoters")
tile.cCHH<-split(complete.cCHH, complete.cCHH$feature)
names(tile.cCHH) <- c("exons", "introns", "promoters")

## Split on promoter/genic
split1.CHH<-split(complete.CHH, !complete.CHH$feature=="promoters")
names(split1.CHH) <- c("promoters", "genic")
split1.aCHH<-split(complete.aCHH, !complete.aCHH$feature=="promoters")
names(split1.aCHH) <- c("promoters", "genic")
split1.bCHH<-split(complete.bCHH, !complete.bCHH$feature=="promoters")
names(split1.bCHH) <- c("promoters", "genic")
split1.cCHH<-split(complete.cCHH, !complete.cCHH$feature=="promoters")
names(split1.cCHH) <- c("promoters", "genic")
#Get unique tiles
feat.CHH <- lapply(split1.CHH, function(x) x[!duplicated(x[c("chr","start")]), ])
feat.aCHH <- lapply(split1.aCHH, function(x) x[!duplicated(x[c("chr","start")]), ])
feat.bCHH <- lapply(split1.bCHH, function(x) x[!duplicated(x[c("chr","start")]), ])
feat.cCHH <- lapply(split1.cCHH, function(x) x[!duplicated(x[c("chr","start")]), ])

## Split on feature type, only unique tiles (no tile that overlaps any features included)
u.CHH<-complete.CHH[!(duplicated(complete.CHH[,1:2]) | duplicated(complete.CHH[,1:2], fromLast = TRUE)), ]
u.aCHH<-complete.aCHH[!(duplicated(complete.aCHH[,1:2]) | duplicated(complete.aCHH[,1:2], fromLast = TRUE)), ]
u.bCHH<-complete.bCHH[!(duplicated(complete.bCHH[,1:2]) | duplicated(complete.bCHH[,1:2], fromLast = TRUE)), ]
u.cCHH<-complete.cCHH[!(duplicated(complete.cCHH[,1:2]) | duplicated(complete.cCHH[,1:2], fromLast = TRUE)), ]
unique.CHH<-split(u.CHH, interaction(!u.CHH$feature=="promoters",u.CHH$meth.diff<0))
unique.CHH$hyperPromoters<-data.frame(name=character(0))
unique.CHH$hyperGenic<-data.frame(name=character(0))
names(unique.CHH) <- c("hypoPromoters","hypoGenic","hyperPromoters","hyperGenic")
unique.aCHH<-split(u.aCHH, interaction(!u.aCHH$feature=="promoters",u.aCHH$meth.diff<0))
names(unique.aCHH) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
unique.bCHH<-split(u.bCHH, interaction(!u.bCHH$feature=="promoters",u.bCHH$meth.diff<0))
names(unique.bCHH) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
unique.cCHH<-split(u.cCHH, interaction(!u.cCHH$feature=="promoters",u.cCHH$meth.diff<0))
names(unique.cCHH) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")

## Split on methylation >0<
split2.CHH<-split(complete.CHH, complete.CHH$meth.diff<0)
split2.CHH$hyper<-data.frame(name=character(0))
names(split2.CHH) <- c("hypo","hyper")
split2.aCHH<-split(complete.aCHH, complete.aCHH$meth.diff<0)
names(split2.aCHH) <- c("hyper", "hypo")
split2.bCHH<-split(complete.bCHH, complete.bCHH$meth.diff<0)
names(split2.bCHH) <- c("hyper", "hypo")
split2.cCHH<-split(complete.cCHH, complete.cCHH$meth.diff<0)
names(split2.cCHH) <- c("hyper", "hypo")
#Retrieve unique names (for easier access in later steps)
hyp.CHH <- lapply(split2.CHH, function(x) x[!duplicated(x[c("name")]), 'name'])
hyp.aCHH <- lapply(split2.aCHH, function(x) x[!duplicated(x[c("name")]), 'name'])
hyp.bCHH <- lapply(split2.bCHH, function(x) x[!duplicated(x[c("name")]), 'name'])
hyp.cCHH <- lapply(split2.cCHH, function(x) x[!duplicated(x[c("name")]), 'name'])

## Frequency of tile:gene mappings
freq.CHH <- complete.CHH[!duplicated(complete.CHH[c("chr","start","name")]),] %>% 
  group_by(name) %>% summarise(counts=n())
freq.aCHH <- complete.aCHH[!duplicated(complete.aCHH[c("chr","start","name")]),] %>% 
  group_by(name) %>% summarise(counts=n())
freq.bCHH <- complete.bCHH[!duplicated(complete.bCHH[c("chr","start","name")]),] %>% 
  group_by(name) %>% summarise(counts=n())
freq.cCHH <- complete.cCHH[!duplicated(complete.cCHH[c("chr","start","name")]),] %>% 
  group_by(name) %>% summarise(counts=n())

## Split on both categories
split3.CHH<-split(complete.CHH, interaction(!complete.CHH$feature=="promoters",complete.CHH$meth.diff<0))
split3.CHH$hyperPromoters<-data.frame(name=character(0))
split3.CHH$hyperGenic<-data.frame(name=character(0))
names(split3.CHH) <- c("hypoPromoters","hypoGenic","hyperPromoters","hyperGenic")
split3.aCHH<-split(complete.aCHH, interaction(!complete.aCHH$feature=="promoters",complete.aCHH$meth.diff<0))
names(split3.aCHH) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
split3.bCHH<-split(complete.bCHH, interaction(!complete.bCHH$feature=="promoters",complete.bCHH$meth.diff<0))
names(split3.bCHH) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
split3.cCHH<-split(complete.cCHH, interaction(!complete.cCHH$feature=="promoters",complete.cCHH$meth.diff<0))
names(split3.cCHH) <- c("hyperPromoters","hyperGenic","hypoPromoters","hypoGenic")
#Retrieve unique names
gene.CHH <- lapply(split3.CHH, function(x) x[!duplicated(x[c("name")]), 'name'])
gene.aCHH <- lapply(split3.aCHH, function(x) x[!duplicated(x[c("name")]), 'name'])
gene.bCHH <- lapply(split3.bCHH, function(x) x[!duplicated(x[c("name")]), 'name'])
gene.cCHH <- lapply(split3.cCHH, function(x) x[!duplicated(x[c("name")]), 'name'])
save.image(file=".RData")

# get highly mapped genes
freqlist <- list(CpG=freq.CpG, aCpG=freq.aCpG, bCpG=freq.bCpG, cCpG=freq.cCpG,
                 CHG=freq.CHG, aCHG=freq.aCHG, bCHG=freq.bCHG, cCHG=freq.cCHG,
                 CHH=freq.CHH, aCHH=freq.aCHH, bCHH=freq.bCHH, cCHH=freq.cCHH)
out<-lapply(freqlist, function(x) x[x$counts >= 8,])