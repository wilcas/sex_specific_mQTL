library(data.table)
library(foreach)
library(doParallel)
library(mediation)
registerDoParallel(16)
get_celltypes <- function(celltype) {
  df <- fread(sprintf("matrix_eqtl_data/mQTL_covar_6_methy_PC_%s.txt", celltype))
  dat <- as.numeric(df[id == celltype])[-c(1)]
  return(dat)
}
argv <- commandArgs(trailingOnly = TRUE)
genotypes <- fread("matrix_eqtl_data/all_imputed_snps_matrixeQTL.txt")
setindex(genotypes, SNP)
methylation <- fread("matrix_eqtl_data/methylation_data_matrixeQTL.txt")
setindex(methylation, V1)
sig_sex <- fread(argv[[1]])
celltypes <- data.frame(
  Stromal = get_celltypes("Stromal"),
  Endothelial = get_celltypes("Endothelial"),
  Hofbauer = get_celltypes("Hofbauer"),
  Trophoblasts = get_celltypes("Trophoblasts"),
  nRBC = get_celltypes("nRBC"),
  Syncytiotrophoblast = get_celltypes("Syncytiotrophoblast")
)


covar <- fread("matrix_eqtl_data/mQTL_covar_6_methy_PC.txt")
cov_mat <- t(covar[, -c(1)])
colnames(cov_mat) <- covar$id
for (celltype in colnames(celltypes)) {
  results_mediation <- foreach(i = 1:nrow(sig_sex)) %dopar% {
    cpg <- sig_sex$gene[i]
    snp <- sig_sex$SNP[i]
    suppressWarnings(l <- as.numeric(genotypes[snp, -c(1), on = "SNP"]))
    suppressWarnings(t <- as.numeric(methylation[cpg, -c(1), on = "V1"]))
    g <- as.numeric(celltypes[, celltype])
    df <- cbind(as.data.frame(cov_mat), data.frame(l = l, g = g, t = t))
    df$cell_int <- l * g
    df$sex_int <- l * df$Sex
    med.fit <- lm(cell_int ~ . - t - cell_int, data = df)
    out.fit <- lm(t ~ . - t, data = df)
    med.out <- try(mediate(med.fit, out.fit, treat = "sex_int", mediator = "cell_int"))
    tmp <- try(data.frame(
      SNP = snp,
      probe = cpg,
      ACME_eff = med.out$d0,
      ACME_low = med.out$d0.ci[1],
      ACME_high = med.out$d0.ci[2],
      ADE_eff = med.out$d1,
      ADE_low = med.out$d1.ci[1],
      ADE_high = med.out$d1.ci[2],
      ACME_p = med.out$d0.p,
      ADE_p = med.out$d1.p
    ))
    if (class(tmp) == "try-error") {
      return(data.frame())
    } else {
      return(tmp)
    }
  }
  fwrite(rbindlist(results_mediation), sprintf("cell_type_mediation_%s_%s_sig.txt", celltype, argv[[2]]), sep = "\t", row.names = F, quote = F)
}
