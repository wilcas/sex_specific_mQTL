library(tidyverse)
library(data.table)
library(coloc)
library(bigsnpr)
library(parallel)
library(glue)
argv <- list("neonatal", 4) # commandArgs(trailingOnly = TRUE)

marginal_bonf <- fread("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/marginal_mcpg_bonf.txt.gz", key = "SNP")[p < 0.05]
rds <- snp_readBed2("/arc/project/st-dennisjk-1/shared/data/1000G_EUR_ldsc_data/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL.bed", backingfile = tempfile())

reference <- snp_attach(rds)

compute_coloc <- function(SNP, mqtl, gwas, method = "coloc", type = "quant", s = NULL) { # SNP must be in SNP column of mqtl and gwas
  D1 <- list(
    beta = mqtl$b,
    varbeta = mqtl$SE^2,
    snp = mqtl$SNP,
    position = mqtl$BP,
    N = 400,
    MAF = mqtl$Freq,
    type = "quant"
  )
  if (is.null(s)) {
    D2 <- list(
      pvalues = pnorm(-abs(gwas$Z)) * 2,
      z = gwas$Z,
      snp = gwas$SNP,
      MAF = snp_MAF(reference$genotypes, ind.col = match(gwas$SNP, reference$map$marker.ID)),
      N = min(gwas$N),
      type = type
    )
  } else {
    D2 <- list(
      pvalues = pnorm(-abs(gwas$Z)) * 2,
      z = gwas$Z,
      snp = gwas$SNP,
      MAF = snp_MAF(reference$genotypes, ind.col = match(gwas$SNP, reference$map$marker.ID)),
      N = min(gwas$N),
      type = type,
      s = s
    )
  }
  if (method == "susie") {
    LD <- snp_cor(reference$genotypes, ind.col = match(SNP, reference$map$marker.ID))^2
    LD <- as.matrix(LD)
    colnames(LD) <- SNP
    rownames(LD) <- SNP
    D1$LD <- LD
    D2$LD <- LD
    S1 <- runsusie(D1)
    S2 <- runsusie(D2)
    return(coloc.susie(S1, S2))
  } else if (method == "coloc") {
    return(coloc.abf(D1, D2))
  } else {
    return(NULL)
  }
}

if (argv[[1]] == "PGC") {
  sumstat_files <- c(
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/adhd_jul2017.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/anxiety.meta.full.fs.tbl.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/anxiety.meta.full.cc.tbl.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/pgcAN2.2019-07.vcf.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/AUDIT_UKB_2018_AJP.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/Cannabis_ICC_23andmetop_UKB_het.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/iPSYCH-PGC_ASD_Nov2017.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/ocd_aug2017.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/pgc_alcdep.eur_unrelated.aug2018_release.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/pgc_bip_2018.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/PGC_UKB_depression_genome-wide.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/MDD2018_ex23andMe.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/PGC3_SCZ_wave3_public.v2.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/tag.cpd.tbl.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/TS_Oct2018.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/pts_all_freeze2_overall.results.sumstats.gz"
  )
  trait_names <- c(
    "ADHD",
    "ANXFS",
    "ANX",
    "AN",
    "AUDIT",
    "CUD",
    "ASD",
    "OCD",
    "ALC",
    "BIP",
    "MDD",
    "MDD2018",
    "SCZ",
    "TAG_CPD",
    "TS",
    "PTSD"
  )
  sample_prev <- c(
    0.36,
    NULL,
    0.33,
    0.23,
    NULL,
    0.30,
    0.40,
    0.28,
    0.26,
    0.15,
    0.34,
    0.35,
    0.42,
    NULL,
    0.34,
    0.15
  )
} else if (argv[[1]] == "neonatal") {
  sumstat_files <- c(
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/Pubertal_growth_PGF_PGM_combined.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/EGG_Obesity_Meta_Analysis_1.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_hay_fever_irnt.gwas.imputed_v3.both_sexes.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/EGG-TotalGWG-Offspring.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/Pubertal_growth_PTF_PTM_combined.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_hay_fever_raw.gwas.imputed_v3.both_sexes.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_asthma_raw.gwas.imputed_v3.both_sexes.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_asthma_irnt.gwas.imputed_v3.both_sexes.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_diabetes_raw.gwas.imputed_v3.both_sexes.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_diabetes_irnt.gwas.imputed_v3.both_sexes.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/EGG_BMI_HapMap_DISCOVERY.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/Pubertal_growth_10F_12M_combined.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/EGG_HC_DISCOVERY.v2.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/CHILD_ONSET_ASTHMA.20180501.allchr.assoc.GC.sumstats.gz"
  )
  trait_names <- c(
    "PGF_PGM",
    "EGG_OBESITY",
    "HAY_FEVER(irnt)",
    "GWG",
    "PTF_PTM",
    "HAY_FEVER(raw)",
    "AGE_ASTHMA(raw)",
    "AGE_ASTHMA(irnt)",
    "AGE_DIABETES(raw)",
    "AGE_DIABETES(irnt)",
    "EGG_BMI",
    "10F_12M",
    "EGG_GC",
    "CHILD_ASTHMA"
  )
  sample_prev <- c(
    NULL,
    0.399,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    0.0288
  )
}
result <- list()
i <- as.numeric(argv[[2]])
if (file.exists(glue("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/{trait_names[i]}_colocalization.txt"))) {
  quit()
}

gwas <- fread(sumstat_files[i])
tmp_marginal_bonf <- marginal_bonf[intersect(gwas$SNP, reference$map$marker.ID), on = "SNP", nomatch = 0]
print("Reading in elligble CpG sites")
eligible_cpg <- unlist(mclapply(
  unique(tmp_marginal_bonf$Probe)[1:100],
  function(probe) {
    mqtl <- tmp_marginal_bonf[Probe == probe]
    gwas_tmp <- gwas[mqtl$SNP, on = "SNP"]
    return(min(mqtl$p) < 5e-8 & max(abs(gwas_tmp$Z)) > 5.45)
  },
  mc.cores = 32
))
if (sum(eligible_cpg) > 0) {
  print("starting colocalization")
  test <- mclapply(
    unique(tmp_marginal_bonf$Probe)[eligible_cpg],
    function(probe) {
      mqtl <- tmp_marginal_bonf[Probe == probe]
      gwas_tmp <- gwas[SNP %in% mqtl$SNP]
      invisible(
        capture.output(res <- compute_coloc(
          mqtl$SNP,
          mqtl,
          gwas_tmp,
          method = "coloc",
          s = sample_prev[i],
          type = ifelse(is.null(sample_prev[i]), "quant", "cc")
        )$summary)
      )
      return(res)
    },
    mc.cores = 32
  )
  names(test) <- unique(tmp_marginal_bonf$Probe)[eligible_cpg]
  result <- rbindlist(lapply(test, function(x) data.table(t(x))), idcol = "probe")

  fwrite(
    result,
    glue("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/{trait_names[i]}_colocalization.txt"),
    quote = F,
    sep = "\t",
    row.names = F
  )
}
