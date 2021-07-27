library(data.table)
library(parallel)
argv <- commandArgs(trailingOnly = TRUE)

all_annotations <- rbindlist(
  mclapply(
    dir(
      sprintf(
        "%s/chr%s",
        argv[[2]],
        argv[[1]]
      ),
      pattern = ".*.out_post",
      full.names = TRUE
    ),
    function(fname) {
      cpg <- gsub(".*/(.*)\\.out_post", "\\1", fname)
      df <- fread(fname)
      df$probe <- cpg
      return(df)
    },
    mc.cores = 32
  )
)
snp_data <- fread("/scratch/st-dennisjk-1/wcasazza/delahaye_QC/matrix_eqtl_data/snp_pos.txt")[CHR == sprintf("chr%s", argv[[1]])]
maxCPP <- all_annotations[Causal_Post._Prob. > 0, .(SNP = SNP_ID, probe, CPP = Causal_Post._Prob.)]
matched <- match(maxCPP$SNP, snp_data$SNP)

maxCPP <- maxCPP[, .(CHR = gsub("chr", "", snp_data$CHR[matched]), BP = snp_data$POS[matched], SNP, probe, CPP)]
fwrite(maxCPP, file = sprintf("%s.%s.txt.gz", argv[[3]], argv[[1]]), row.names = F, quote = F, sep = "\t")
