library(data.table)
library(parallel)
argv <- commandArgs(trailingOnly=TRUE)
probe_pos <- fread("matrix_eqtl_data/probe_pos.txt")
kg_snps <- fread(
  "/arc/project/st-dennisjk-1/shared/data/1000G_phase3_info/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz",
  skip="#CHROM",
  drop=c("QUAL","FILTER","INFO")
)
all_mQTLs <- fread(argv[[2]])
mQTL_in_1kg <- all_mQTLs[SNP %chin% kg_snps[`#CHROM` == argv[[1]]]$ID]
out <- mclapply(
  unique(mQTL_in_1kg$gene),
  function(cpg){
    dt <- mQTL_in_1kg[gene == cpg & `p-value` < 0.05,.(SNP,Z =`t-stat`)]
    if(nrow(dt) > 0){
      fwrite(
        dt,
        sprintf("caviar_files/%s/chr%s/%s.z",argv[[3]],argv[[1]],cpg),
        sep="\t",
        col.names=FALSE,
        quote=FALSE
      )
      return(data.table())
    }else{
      return(data.table())
    }
  },
  mc.cores=16
)
