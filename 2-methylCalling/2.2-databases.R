###################################################################
### METHYLKIT METHYLATION CALLING
### R script for methylation tiling and calling
###################################################################

# set up environment
setwd("/scratch/nia/manuFinal")
load(file=".RData")
library(methylKit)
library(GenomicRanges)
library(genomation)

# prepare GRanges object of regions not affected by introgression
chr2region <- GRanges(seqnames=c("1", "2", "2", "3", "4", "5", "6", "7", "8", "9", "10"), 
              ranges=IRanges(start=c(1,1,182000000,1,1,1,1,1,1,1,1),
                             end=c(307041717,22000000,244442276,235667834,246994605,223902240,174033170,182381542,181122637,159769782,150982314)),
              strand="*")
seqlengths(chr2region)=c(307041717,244442276,235667834,246994605,223902240,174033170,182381542,181122637,159769782,150982314)

# CpG ---------------------------------------------------------------------
print("CpG ----------------------------------------------------------------")
# Tile into 100bp windows
CpG.3x.100.preSubset=tileMethylCounts(CpG.3x,win.size=100,step.size=100,cov.bases=1,mc.cores=10,suffix="100")
# Remove introgression-affected region
CpG.3x.100.subset=selectByOverlap(CpG.3x.100.preSubset, chr2region)
CpG.3x.100=makeMethylDB(CpG.3x.100.subset, "subsetDB")
# Reorganize into pairwise genotype-specific comparisons
aCpG.3x.100=reorganize(CpG.3x.100, sample.ids=c("aMutCpG", "aWTCpG"), treatment=c(1,0), dbdir="subsetDB")
bCpG.3x.100=reorganize(CpG.3x.100, sample.ids=c("bMutCpG", "bWTCpG"), treatment=c(1,0), dbdir="subsetDB")
cCpG.3x.100=reorganize(CpG.3x.100, sample.ids=c("cMutCpG", "cWTCpG"), treatment=c(1,0), dbdir="subsetDB")
# Unite calls for each comparison
meth.CpG.3x.100=unite(CpG.3x.100, destrand=FALSE, min.per.group = 2L, mc.cores=10,suffix="CpG_3x_100")
meth.aCpG.3x.100=unite(aCpG.3x.100, destrand=FALSE, mc.cores=10,suffix="aCpG_3x_100")
meth.bCpG.3x.100=unite(bCpG.3x.100, destrand=FALSE, mc.cores=10,suffix="bCpG_3x_100")
meth.cCpG.3x.100=unite(cCpG.3x.100, destrand=FALSE, mc.cores=10,suffix="cCpG_3x_100")
# Calculate differential methylation
diff.CpG.3x.100.nocorr=calculateDiffMeth(meth.CpG.3x.100, mc.cores=10, suffix="nocorr") # no overdispersion correction, full model
diff.CpG.3x.100.oChi=calculateDiffMeth(meth.CpG.3x.100, mc.cores=10,overdispersion="MN", test="Chisq",suffix="oChi") # chi2 overdispersion correction, full model
diff.CpG.3x.100.oF=calculateDiffMeth(meth.CpG.3x.100, overdispersion="MN", mc.cores=10, suffix="oF") # F-test overdispersion correction, full model
diff.aCpG.3x.100=calculateDiffMeth(meth.aCpG.3x.100, mc.cores=10) # E genotype-specific pairwise comparison
diff.bCpG.3x.100=calculateDiffMeth(meth.bCpG.3x.100, mc.cores=10) # G genotype-specific pairwise comparison
diff.cCpG.3x.100=calculateDiffMeth(meth.cCpG.3x.100, mc.cores=10) # H genotype-specific pairwise comparison
# Call significant differential methylation
CpG.3x.100.nocorr.50p.all=getMethylDiff(diff.CpG.3x.100.nocorr,difference=50, qvalue=0.01, type="all", suffix="50p_all") # no overdispersion correction, full model
CpG.3x.100.oChi.50p.all=getMethylDiff(diff.CpG.3x.100.oChi,difference=50, qvalue=0.01, type="all", suffix="50p_all") # chi2 overdispersion correction, full model
CpG.3x.100.oF.50p.all=getMethylDiff(diff.CpG.3x.100.oF,difference=50, qvalue=0.01, type="all", suffix="50p_all") # F-test overdispersion correction, full model
aCpG.3x.100.50p.all=getMethylDiff(diff.aCpG.3x.100,difference=50, qvalue=0.01, type="all", suffix="50p_all") # E genotype-specific pairwise comparison
bCpG.3x.100.50p.all=getMethylDiff(diff.bCpG.3x.100,difference=50, qvalue=0.01, type="all", suffix="50p_all") # G genotype-specific pairwise comparison
cCpG.3x.100.50p.all=getMethylDiff(diff.cCpG.3x.100,difference=50, qvalue=0.01, type="all", suffix="50p_all") # H genotype-specific pairwise comparison

# CHG ---------------------------------------------------------------------
print("CHG ----------------------------------------------------------------")
# Tile into 100bp windows
CHG.3x.100.preSubset=tileMethylCounts(CHG.3x,win.size=100,step.size=100,cov.bases=1,mc.cores=10,suffix="100")
# Remove introgression-affected region
CHG.3x.100.subset=selectByOverlap(CHG.3x.100.preSubset, chr2region)
CHG.3x.100=makeMethylDB(CHG.3x.100.subset, "subsetDB")
# Reorganize into pairwise genotype-specific comparisons
aCHG.3x.100=reorganize(CHG.3x.100, sample.ids=c("aMutCHG", "aWTCHG"), treatment=c(1,0), dbdir="subsetDB")
bCHG.3x.100=reorganize(CHG.3x.100, sample.ids=c("bMutCHG", "bWTCHG"), treatment=c(1,0), dbdir="subsetDB")
cCHG.3x.100=reorganize(CHG.3x.100, sample.ids=c("cMutCHG", "cWTCHG"), treatment=c(1,0), dbdir="subsetDB")
# Unite calls for each comparison
meth.CHG.3x.100=unite(CHG.3x.100, destrand=FALSE, min.per.group = 2L, mc.cores=10,suffix="CHG_3x_100")
meth.aCHG.3x.100=unite(aCHG.3x.100, destrand=FALSE, mc.cores=10,suffix="aCHG_3x_100")
meth.bCHG.3x.100=unite(bCHG.3x.100, destrand=FALSE, mc.cores=10,suffix="bCHG_3x_100")
meth.cCHG.3x.100=unite(cCHG.3x.100, destrand=FALSE, mc.cores=10,suffix="cCHG_3x_100")
# Calculate differential methylation
diff.CHG.3x.100.nocorr=calculateDiffMeth(meth.CHG.3x.100, mc.cores=10, suffix="nocorr") # no overdispersion correction, full model
diff.CHG.3x.100.oChi=calculateDiffMeth(meth.CHG.3x.100, mc.cores=10,overdispersion="MN", test="Chisq",suffix="oChi") # chi2 overdispersion correction, full model
diff.CHG.3x.100.oF=calculateDiffMeth(meth.CHG.3x.100, overdispersion="MN", mc.cores=10, suffix="oF") # F-test overdispersion correction, full model
diff.aCHG.3x.100=calculateDiffMeth(meth.aCHG.3x.100, mc.cores=10) # E genotype-specific pairwise comparison
diff.bCHG.3x.100=calculateDiffMeth(meth.bCHG.3x.100, mc.cores=10) # G genotype-specific pairwise comparison
diff.cCHG.3x.100=calculateDiffMeth(meth.cCHG.3x.100, mc.cores=10) # H genotype-specific pairwise comparison
# Call significant differential methylation
CHG.3x.100.nocorr.50p.all=getMethylDiff(diff.CHG.3x.100.nocorr,difference=50, qvalue=0.01, type="all", suffix="50p_all") # no overdispersion correction, full model
CHG.3x.100.oChi.50p.all=getMethylDiff(diff.CHG.3x.100.oChi,difference=50, qvalue=0.01, type="all", suffix="50p_all") # chi2 overdispersion correction, full model
CHG.3x.100.oF.50p.all=getMethylDiff(diff.CHG.3x.100.oF,difference=50, qvalue=0.01, type="all", suffix="50p_all") # F-test overdispersion correction, full model
aCHG.3x.100.50p.all=getMethylDiff(diff.aCHG.3x.100,difference=50, qvalue=0.01, type="all", suffix="50p_all") # E genotype-specific pairwise comparison
bCHG.3x.100.50p.all=getMethylDiff(diff.bCHG.3x.100,difference=50, qvalue=0.01, type="all", suffix="50p_all") # G genotype-specific pairwise comparison
cCHG.3x.100.50p.all=getMethylDiff(diff.cCHG.3x.100,difference=50, qvalue=0.01, type="all", suffix="50p_all") # H genotype-specific pairwise comparison

# CHH ---------------------------------------------------------------------
print("CHH ----------------------------------------------------------------")
# Tile into 100bp windows
CHH.3x.100.preSubset=tileMethylCounts(CHH.3x,win.size=100,step.size=100,cov.bases=1,mc.cores=10,suffix="100")
# Remove introgression-affected region
CHH.3x.100.subset=selectByOverlap(CHH.3x.100.preSubset, chr2region)
CHH.3x.100=makeMethylDB(CHH.3x.100.subset, "subsetDB")
# Reorganize into pairwise genotype-specific comparisons
aCHH.3x.100=reorganize(CHH.3x.100, sample.ids=c("aMutCHH", "aWTCHH"), treatment=c(1,0), dbdir="subsetDB")
bCHH.3x.100=reorganize(CHH.3x.100, sample.ids=c("bMutCHH", "bWTCHH"), treatment=c(1,0), dbdir="subsetDB")
cCHH.3x.100=reorganize(CHH.3x.100, sample.ids=c("cMutCHH", "cWTCHH"), treatment=c(1,0), dbdir="subsetDB")
# Unite calls for each comparison
meth.CHH.3x.100=unite(CHH.3x.100, destrand=FALSE, min.per.group = 2L, mc.cores=10,suffix="CHH_3x_100")
meth.aCHH.3x.100=unite(aCHH.3x.100, destrand=FALSE, mc.cores=10,suffix="aCHH_3x_100")
meth.bCHH.3x.100=unite(bCHH.3x.100, destrand=FALSE, mc.cores=10,suffix="bCHH_3x_100")
meth.cCHH.3x.100=unite(cCHH.3x.100, destrand=FALSE, mc.cores=10,suffix="cCHH_3x_100")
# Calculate differential methylation
diff.CHH.3x.100.nocorr=calculateDiffMeth(meth.CHH.3x.100, mc.cores=10, suffix="nocorr") # no overdispersion correction, full model
diff.CHH.3x.100.oChi=calculateDiffMeth(meth.CHH.3x.100, mc.cores=10,overdispersion="MN", test="Chisq",suffix="oChi") # chi2 overdispersion correction, full model
diff.CHH.3x.100.oF=calculateDiffMeth(meth.CHH.3x.100, overdispersion="MN", mc.cores=10, suffix="oF") # F-test overdispersion correction, full model
diff.aCHH.3x.100=calculateDiffMeth(meth.aCHH.3x.100, mc.cores=10) # E genotype-specific pairwise comparison
diff.bCHH.3x.100=calculateDiffMeth(meth.bCHH.3x.100, mc.cores=10) # G genotype-specific pairwise comparison
diff.cCHH.3x.100=calculateDiffMeth(meth.cCHH.3x.100, mc.cores=10) # H genotype-specific pairwise comparison
# Call significant differential methylation
CHH.3x.100.nocorr.50p.all=getMethylDiff(diff.CHH.3x.100.nocorr,difference=50, qvalue=0.01, type="all", suffix="50p_all") # no overdispersion correction, full model
CHH.3x.100.oChi.50p.all=getMethylDiff(diff.CHH.3x.100.oChi,difference=50, qvalue=0.01, type="all", suffix="50p_all") # chi2 overdispersion correction, full model
CHH.3x.100.oF.50p.all=getMethylDiff(diff.CHH.3x.100.oF,difference=50, qvalue=0.01, type="all", suffix="50p_all") # F-test overdispersion correction, full model
aCHH.3x.100.50p.all=getMethylDiff(diff.aCHH.3x.100,difference=50, qvalue=0.01, type="all", suffix="50p_all") # E genotype-specific pairwise comparison
bCHH.3x.100.50p.all=getMethylDiff(diff.bCHH.3x.100,difference=50, qvalue=0.01, type="all", suffix="50p_all") # G genotype-specific pairwise comparison
cCHH.3x.100.50p.all=getMethylDiff(diff.cCHH.3x.100,difference=50, qvalue=0.01, type="all", suffix="50p_all") # H genotype-specific pairwise comparison

# save workspace for later loading
save.image(file=".RData")
save.image(file="backupRData/2.2-databases-backup.RData")