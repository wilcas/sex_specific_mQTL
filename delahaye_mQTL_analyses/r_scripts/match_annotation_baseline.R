library(data.table)
library(parallel)


argv <- commandArgs(trailingOnly = TRUE)
annot_str <- argv[[1]]
out <- mclapply(
  1:22,
  function(i) {
    annot <- fread(sprintf(annot_str, i))
    annot <- annot[!duplicated(SNP)]
    bim <- fread(sprintf("/arc/project/st-dennisjk-1/shared/data/1000G_phase3_by_chr/plink/mapped_maf01_chr%d.bim", i))[, .(V1, V4, V2, V3)]
    colnames(bim) <- c("CHR", "BP", "SNP", "CM")
    merged <- merge(bim, annot, by = colnames(bim), all.x = TRUE)
    for (j in 5:ncol(annot)) {
      to_select <- c("CHR", "BP", "SNP", "CM", colnames(annot)[j])
      fwrite(
        merged[, ..to_select],
        sprintf(argv[[2]], colnames(annot)[j], i),
        row.names = F,
        sep = "\t",
        quote = F
      )
    }
    return(c(nrow(bim), nrow(merged), nrow(annot)))
  },
  mc.cores = 8
)
print(out)
