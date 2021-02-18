library(data.table)
library(parallel)

argv <- commandArgs(trailingOnly = TRUE)
mqtl <- fread(argv[[1]], key = "SNP")[FDR < 0.05]
annot_str <- argv[[2]]
out_suffix <- argv[[3]]

out <- mclapply(
  1:22,
  function(i){
    new_annot <- fread(sprintf(annot_str,i))[,c(1,2,3,4)]
    new_annot$all_cis <- as.numeric(new_annot$SNP %chin% mqtl$SNP)
    fwrite(new_annot,sprintf("%s.%d.annot.gz",out_suffix,i),sep='\t',quote=F,row.names=F)
  },
  mc.cores=8
)
