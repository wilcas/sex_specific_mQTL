library(data.table)
library(parallel)


argv <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(argv[[3]])
annot_str <- argv[[1]]
annot <- fread(sprintf(annot_str, i))
annot$CM <- NULL
annot <- annot[!duplicated(SNP)]
bim <- fread(sprintf("/arc/project/st-dennisjk-1/shared/data/1000G_EUR_ldsc_data/1000G_EUR_Phase3_plink/1000G.EUR.QC.%d.bim", i))[, .(V1, V4, V2, V3)]
#bim <- fread(sprintf("/arc/project/st-dennisjk-1/shared/data/1000G_phase3_by_chr/plink/mapped_maf01_chr%d.bim", i))[, .(V1, V4, V2, V3)]
colnames(bim) <- c("CHR", "BP", "SNP", "CM")
merged <- merge(bim, annot, by = colnames(bim)[-c(4)], all.x = TRUE)
setnafill(merged, fill = 0, cols = 5:ncol(merged))
for (j in 4:ncol(annot)) {
  print(j)
  to_select <- c("CHR", "BP", "SNP", "CM", colnames(annot)[j])
  fwrite(
    merged[, ..to_select],
    sprintf(argv[[2]], colnames(annot)[j], i),
    row.names = F,
    sep = "\t",
    quote = F
  )
}
c(nrow(bim), nrow(merged), nrow(annot))
