###################################################################
### G WILD-TYPE PROCESSBISMARKALN
###################################################################

setwd("/scratch/nia/manuFinal")
library(methylKit)

processBismarkAln(location = "/home/nia/project/nia/alignments/bWT/B_10.bam",
                  sample.id="bWT", assembly="b73", read.context="none",
                  save.context=c("CpG","CHH","CHG"), save.folder=getwd(),
                  nolap=FALSE, mincov=1, minqual=20)