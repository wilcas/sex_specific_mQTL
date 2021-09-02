library(data.table)
library(batchtools)

probes <- fread("/scratch/st-dennisjk-1/wcasazza/delahaye_QC/matrix_eqtl_data/probe_pos.txt")[chr != "chrX"]
run_reml <- function(j) {
  library(data.table)
  probes <- fread("/scratch/st-dennisjk-1/wcasazza/delahaye_QC/matrix_eqtl_data/probe_pos.txt")[chr != "chrX"]
  snps <- fread("/scratch/st-dennisjk-1/wcasazza/delahaye_QC/matrix_eqtl_data/snp_pos.txt")

  all_cpgs <- colnames(fread("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_methy_mvalues.phen", nrows = 1))[-c(1, 2)]
  grm_cmd <- paste(
    "/arc/project/st-dennisjk-1/software/plink2.0/plink2",
    "--bfile", "/scratch/st-dennisjk-1/wcasazza/delahaye_QC/placenta_regulatory_landscape/RootStudyConsentSet_phs001717.PlacentalRegulation.v1.p1.c1.HMB-IRB-PUB-COL-MDS/genotype_qc/imputed_geno_for_heritability_maf01",
    "--extract", "/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_snps/%s_snps.txt",
    "--thread-num", "16",
    "--make-grm-list",
    "--out", "/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_grm/%s",
    "&& gzip /scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_grm/%s.grm"
  )
  reacta_cmd <- paste(
    "/arc/project/st-dennisjk-1/software/reacta",
    "--grm", "/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_grm/%s",
    "--reml",
    "--pheno", "/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_methy_mvalues_noheader_IDs.phen",
    "--qcovar", "/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_covariates.qcovar",
    "--mpheno", "%d",
    "--thread-num", "16",
    "--out", "/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_output/%s"
  )
  parallelMap::parallelMap(
    function(i) {
      cpg <- probes$geneid[i]
      s1 <- probes$s1[i]
      s2 <- probes$s2[i]
      chr <- probes$chr[i]
      mpheno <- which(all_cpgs == cpg)
      if (!file.exists(sprintf("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_snps/%s_snps.txt", cpg))) {
        fwrite(
          snps[CHR == chr & (POS >= s1 - 75000 & POS <= s2 + 75000), .(SNP)],
          sprintf("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_snps/%s_snps.txt", cpg),
          quote = F,
          col.names = F,
          row.names = F
        )
      }
      grm_file <- sprintf("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_grm/%s.grm.gz", cpg)
      hsq_file <- sprintf("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_output/%s.hsq", cpg)
      if (!file.exists(grm_file) & !file.exists(hsq_file)) {
        tmp <- sprintf(grm_cmd, cpg, cpg, cpg)
        system(tmp, intern = FALSE)

        tmp <- sprintf(reacta_cmd, cpg, mpheno, cpg)
        system(tmp, intern = FALSE)
      }
      system(sprintf("rm /scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_grm/%s*", cpg))
      system(sprintf("rm /scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/gcta_analysis/delahaye_snps/%s_snps.txt", cpg))
    },
    j
  )
}

unlink("/scratch/st-dennisjk-1/wcasazza/registry_default/", recursive = TRUE)
tmp <- makeRegistry("/scratch/st-dennisjk-1/wcasazza/registry_default", make.default = TRUE)
tmp$cluster.functions <- makeClusterFunctionsTORQUE("~/.config/batchtools/torque-lido.tmpl")
idx <- 1:nrow(probes)
job_ids <- batchMap(fun = run_reml, split(idx, cut(idx, 12)), reg = tmp)
res <- list(walltime = 60 * 60 * 24, memory = 128, pm.backend = "multicore", ncpus = 32, allocation = "st-dennisjk-1", conda.env = "~/miniconda3/envs/misc_bio/bin/activate")
submitJobs(ids = job_ids, res = res)
