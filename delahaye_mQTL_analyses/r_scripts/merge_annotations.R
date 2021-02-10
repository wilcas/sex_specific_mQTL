################################################################################
# Merge two QTL based annotations and add a FDR < 0.05 annotation
################################################################################
library(data.table)
library(parallel)

argv <- commandArgs(trailingOnly = TRUE)
print(argv)
annot1 <- argv[[1]]
suffix1 <- argv[[2]]
mqtl1 <- fread(argv[[3]], key = "SNP")[FDR < 0.05]
print("read mqtl1")
annot2 <- argv[[4]]
suffix2 <- argv[[5]]
mqtl2 <- fread(argv[[6]], key = "SNP")[FDR < 0.05]
print("read mqtl2")
out <- mclapply(
  1:22,
  function(i) {
    a1 <- fread(sprintf(annot1, i))
    a1$all_cis <- as.numeric(a1$SNP %chin% mqtl1$SNP)
    colnames(a1)[5:ncol(a1)] <- paste0(colnames(a1)[5:ncol(a1)], suffix1)
    a2 <- fread(sprintf(annot2, i))
    a2$all_cis <- as.numeric(a2$SNP %chin% mqtl2$SNP)
    colnames(a2)[5:ncol(a2)] <- paste0(colnames(a2)[5:ncol(a2)], suffix2)
    fwrite(
      merge(a1, a2, by = c("CHR", "BP", "SNP", "CM")),
      sprintf("merged%s%s.%d.annot.gz", suffix1, suffix2, i),
      row.names = F,
      sep = "\t",
      quote = F
    )
  },
  mc.cores = 8
)
