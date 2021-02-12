library(data.table)
library(parallel)


argv <- commandArgs(trailingOnly=TRUE)
annot_str <- argv[[1]]
out <- mclapply(
  1:22,
  function(i) {
    baseline <- fread(sprintf("/scratch/st-dennisjk-1/wcasazza/1000G_phase3_ldsc/baseline_annot_files/baselineLD.%d.annot.gz",i))
    annot <- fread(sprintf(annot_str, i))
    for(j in 5:ncol(annot)){
      fwrite(
        merge(baseline,annot, by = c("CHR", "BP", "SNP", "CM"))[,c("CHR", "BP", "SNP", "CM",colnames(annot)[j])],
        sprintf(paste0(colnames(annot)[j],"_",argv[[2]]), i),
        row.names = F,
        sep = "\t",
        quote = F
      )
    }
  },
  mc.cores = 8
)
