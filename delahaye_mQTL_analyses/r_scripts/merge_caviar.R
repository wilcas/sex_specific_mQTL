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
    fread,
    mc.cores = 32
  )
)
snp_data <- fread("/scratch/st-dennisjk-1/wcasazza/delahaye_QC/matrix_eqtl_data/snp_pos.txt")[CHR == sprintf("chr%s", argv[[1]])]
maxCPP <- all_annotations[, .(SNP = SNP_ID, maxCPP = max(Causal_Post._Prob.)), by = "SNP_ID"]
matched <- match(maxCPP$SNP, snp_data$SNP)

maxCPP <- maxCPP[, .(CHR = gsub("chr", "", snp_data$CHR[matched]), BP = snp_data$POS[matched], SNP = SNP, CM = 0, maxCPP, bin_CPP = as.numeric(maxCPP != 0))]
fwrite(maxCPP, file = sprintf("%sdelahaye_caviar_cpp_annotations_chr%s.annot.gz", argv[[3]], argv[[1]]), row.names = F, quote = F, sep = "\t")
