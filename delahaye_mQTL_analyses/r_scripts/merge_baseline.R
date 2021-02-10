library(data.table)
library(parallel)

setwd("/scratch/st-dennisjk-1/wcasazza/1000G_phase3_ldsc/baseline_annot_files/")

out <- mclapply(
  1:22,
  function(i){
    files <- dir(sprintf("chr%d",i),full.names=T)
    names <- gsub(".*chr.*/(.*)\\.annot", "\\1",files)
    dt <- do.call(cbind, lapply(files, function(f)fread(f,header=T)))
    colnames(dt) <- names
    dt$base <- 1
    bim <- fread(sprintf("/project/st-dennisjk-1/shared/data/1000G_phase3_by_chr/plink/mapped_maf05_chr%d.bim", i),header=F)
    bim <- bim[,.(CHR=i,BP=V4,SNP=V2,CM=V3)]
    fwrite(cbind(bim,dt), sprintf("baselineLD.%d.annot.gz",i),sep="\t",col.names=T,row.names=F,quote=F)
  },
  mc.cores = 16
)
