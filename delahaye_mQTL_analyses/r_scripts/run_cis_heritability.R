library(data.table)
library(parallel)
library(batchtools)

run_reml <- function(i) {
  probes <- fread("/scratch/st-dennisjk-1/wcasazza/delahaye_QC/matrix_eqtl_data/probe_pos.txt")
  snps <- fread("/scratch/st-dennisjk-1/wcasazza/delahaye_QC/matrix_eqtl_data/snp_pos.txt")

  all_cpgs <- colnames(fread("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_methy_mvalues_noheader_IDs.phen", nrows = 1))[-c(1, 2)]
  grm_cmd <- paste(
    "/arc/project/st-dennisjk-1/software/plink2.0/plink2",
    "--bfile", "/scratch/st-dennisjk-1/wcasazza/delahaye_QC/placenta_regulatory_landscape/RootStudyConsentSet_phs001717.PlacentalRegulation.v1.p1.c1.HMB-IRB-PUB-COL-MDS/genotype_qc/imputed_geno_for_heritability_maf01",
    "--extract", "/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_snps/%s_snps.txt",
    "--thread-num", "16",
    "--make-grm-list",
    "--out", "/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_grm/%s",
    "&& gzip /scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_grm/%s.grm"
  )
  cmd <- paste(
    "/arc/project/st-dennisjk-1/software/reacta",
    "--grm", "~/analysis_of_pd_dnam/gcta_analysis/delahaye_grm/%s",
    "--reml",
    "--pheno", "/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_methy_mvalues_noheader_IDs.phen",
    "--qcovar", "/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_covariates.qcovar",
    "--mpheno", "%d",
    "--thread-num", "16",
    "--out", "~/analysis_of_pd_dnam/gcta_analysis/delahaye_output/%s"
  )
  cpg <- probes$geneid[i]
  s1 <- probes$s1[i]
  s2 <- probes$s2[i]
  chr <- probes$chr[i]
  mpheno <- which(all_cpgs == cpg)
  if (!file.exists(sprintf("gcta_analysis/delahaye_snps/%s_snps.txt", cpg))) {
    fwrite(
      snps[CHR == chr & (POS >= s1 - 75000 & POS <= s2 + 75000), .(SNP)],
      sprintf("gcta_analysis/delahaye_snps/%s_snps.txt", cpg),
      quote = F,
      col.names = F,
      row.names = F
    )
  }
  if (!file.exists(sprintf("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_grm/%s.grm.gz", cpg))) {
    tmp <- sprintf(grm_cmd, cpg, cpg, cpg)
    system(tmp, intern = FALSE)
  }
  if (!file.exists(sprintf("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_output/%s.hsq", cpg))) {
    tmp <- sprintf(cmd, cpg, mpheno, cpg)
    system(tmp, intern = FALSE)
  }
  system(sprintf("rm /scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_grm/%s*", cpg))
  system(sprintf("rm /scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_snps/%s_snps.txt", cpg))
}
i <- splitIndices(1:nrow(probes))
# @TODO batch submit to cluster code
out <- mclapply(1:nrow(probes), run_reml, mc.cores = 32)
