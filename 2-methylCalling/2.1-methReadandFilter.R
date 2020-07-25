###################################################################
### METHYLKIT METHREAD AND FILTERING
### R script to create tabix files of filtered cytosine methylation
### by sequence context
###################################################################

# set up environment
setwd("/scratch/nia/manuFinal")
library(methylKit)
load(file=".RData")

# read in methylation proportion files and filter
print("CpG ----------------------------------------------------------------")
CpG.raw=methRead(list("/scratch/nia/manuFinal/aMut_CpG.txt", "/scratch/nia/manuFinal/bMut_CpG.txt",
                      "/scratch/nia/manuFinal/cMut_CpG.txt", "/scratch/nia/manuFinal/aWT_CpG.txt",
                      "/scratch/nia/manuFinal/bWT_CpG.txt", "/scratch/nia/manuFinal/cWT_CpG.txt"),
                 sample.id=list("aMutCpG","bMutCpG","cMutCpG","aWTCpG", "bWTCpG", "cWTCpG"),
                 assembly="b73", treatment=c(1,1,1,0,0,0), context="CpG", dbtype = "tabix", dbdir = "methylDB", mincov=1)

CpG.3x=filterByCoverage(CpG.raw,lo.count=3,lo.perc=NULL,hi.count=NULL,hi.perc=99.9, suffix="3x", dbdir="methylDB")

print("CHG ----------------------------------------------------------------")
CHG.raw=methRead(list("/scratch/nia/manuFinal/aMut_CHG.txt", "/scratch/nia/manuFinal/bMut_CHG.txt",
                      "/scratch/nia/manuFinal/cMut_CHG.txt", "/scratch/nia/manuFinal/aWT_CHG.txt",
                      "/scratch/nia/manuFinal/bWT_CHG.txt", "/scratch/nia/manuFinal/cWT_CHG.txt"),
                 sample.id=list("aMutCHG","bMutCHG","cMutCHG","aWTCHG", "bWTCHG", "cWTCHG"),
                 assembly="b73", treatment=c(1,1,1,0,0,0), context="CHG", dbtype = "tabix", dbdir = "methylDB", mincov=1)

CHG.3x=filterByCoverage(CHG.raw,lo.count=3,lo.perc=NULL,hi.count=NULL,hi.perc=99.9, suffix="3x", dbdir="methylDB")

print("CHH ----------------------------------------------------------------")
CHH.raw=methRead(list("/scratch/nia/manuFinal/aMut_CHH.txt", "/scratch/nia/manuFinal/bMut_CHH.txt",
                      "/scratch/nia/manuFinal/cMut_CHH.txt", "/scratch/nia/manuFinal/aWT_CHH.txt",
                      "/scratch/nia/manuFinal/bWT_CHH.txt", "/scratch/nia/manuFinal/cWT_CHH.txt"),
                 sample.id=list("aMutCHH","bMutCHH","cMutCHH","aWTCHH", "bWTCHH", "cWTCHH"),
                 assembly="b73", treatment=c(1,1,1,0,0,0), context="CHH", dbtype = "tabix", dbdir = "methylDB", mincov=1)

CHH.3x=filterByCoverage(CHH.raw,lo.count=3,lo.perc=NULL,hi.count=NULL,hi.perc=99.9, suffix="3x", dbdir="methylDB")

# save workspace image for later loading
save.image(file=".RData")
save.image(file="backupRData/2.1-backup.RData")

q(save="yes")