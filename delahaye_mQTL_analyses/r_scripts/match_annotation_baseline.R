library(data.table)
library(parallel)


argv <- commandArgs(trailingOnly=TRUE)
annot_str <- argv[[1]]
out <- mclapply(
  1:22,
  function(i) {
    baseline <- fread(sprintf("/scratch/st-dennisjk-1/wcasazza/1000G_phase3_ldsc/baseline_annot_files/baselineLD.%d.annot.gz",i))
    annot <- fread(sprintf(annot_str, i))
    merged <- merge(baseline,annot, by = c("CHR", "BP", "SNP", "CM"))
    for(j in 5:ncol(annot)){
      to_select <- c("CHR", "BP", "SNP", "CM", colnames(annot)[j])
      fwrite(
        merged[,..to_select],
        sprintf(argv[[2]],colnames(annot)[j], i),
        row.names = F,
        sep = "\t",
        quote = F
      )
    }
  },
  mc.cores = 8
)
