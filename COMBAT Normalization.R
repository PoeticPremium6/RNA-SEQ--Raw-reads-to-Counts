install.packages("COMBAT")
library("COMBAT")

# read SNP P values
file1 <- paste(path.package("COMBAT"),"extdata","SNP_info.txt.gz",sep="/")
snp.info  <- read.table(file1, header = TRUE, as.is=TRUE)
snp.pvals <- as.matrix(snp.info[,2])

# read reference genotype
file2 <- paste(path.package("COMBAT"),"extdata","SNP_ref.txt.gz",sep="/")
snp.ref   <- read.table(file2, header = TRUE)
snp.ref   <- as.matrix(snp.ref)
#call COMBAT

COMBAT(snp.pvals, snp.ref, nperm=100, ncores=2)

install.packages("deSeq2")

countMtx <- read.table("gene.matrix")

setwd("~/Documents/Comp_Micro/Tuxedo/htseq_counts")
count = cbind(UTI44.4A_htseq_counts.out, UTI44.4B_htseq_counts_pos_sorted.out, UTI44.4A_htseq_counts_pos_sorted.out)
