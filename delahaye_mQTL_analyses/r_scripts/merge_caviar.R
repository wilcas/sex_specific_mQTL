library(data.table)
library(parallel)
argv <- commandArgs(trailingOnly = TRUE)

if (file.exists(sprintf("%sdelahaye_caviar_cpp_annotations_chr%s.annot.gz", argv[[2]], argv[[1]]))) {
  print("fixing annotation")
  tmp_annot <- fread(sprintf("%sdelahaye_caviar_cpp_annotations_chr%s.annot.gz", argv[[2]], argv[[1]]))
  tmp <- fread(sprintf("/arc/project/st-dennisjk-1/shared/data/1000G_phase3_by_chr/plink/chr%s_mapped.bim", argv[[1]]))
  merged <- merge(tmp_annot, tmp, by.x = "SNP", by.y = "V2", all = T)
  merged <- merged[, .(CHR = argv[[1]], BP = V4, SNP = SNP, CM = V3, maxCPP, bin_CPP)]
  merged <- merged[match(tmp$V2, SNP)]
  setnafill(merged, fill = 0, cols = c("maxCPP", "bin_CPP"))
  fwrite(merged, sprintf("%s.%s.annot.gz", argv[[2]], argv[[1]]), row.names = F, quote = F, sep = "\t")
} else {
  all_annotations <- rbindlist(
    mclapply(
      dir(
        sprintf(
          "caviar_files/%s/chr%s",
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
  snp_data <- fread("matrix_eqtl_data/snp_pos.txt")[CHR == sprintf("chr%s", argv[[1]])]
  maxCPP <- all_annotations[, .(SNP = SNP_ID, maxCPP = max(Causal_Post._Prob.)), by = "SNP_ID"]
  matched <- match(maxCPP$SNP, snp_data$SNP)

  maxCPP <- maxCPP[, .(CHR = gsub("chr", "", snp_data$CHR[matched]), BP = snp_data$POS[matched], SNP = SNP, CM = 0, maxCPP, bin_CPP = as.numeric(maxCPP != 0))]
  fwrite(maxCPP, file = sprintf("%sdelahaye_caviar_cpp_annotations_chr%s.annot.gz", argv[[2]], argv[[1]]), row.names = F, quote = F, sep = "\t")
}
