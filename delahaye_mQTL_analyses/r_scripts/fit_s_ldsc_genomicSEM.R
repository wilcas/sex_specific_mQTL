library(GenomicSEM)
library(data.table)
library(glue)
sumstat_files <- c(
  "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_sumstats/adhd_jul2017.gz",
  "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_sumstats/anxiety.meta.full.fs.tbl.gz",
  "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_sumstats/anxiety.meta.full.cc.tbl.gz",
  "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_sumstats/pgcAN2.2019-07_refmt.vcf.tsv.gz",
  "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_sumstats/AUDIT_UKB_2018_AJP_fixed.txt.gz",
  "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_sumstats/Cannabis_ICC_23andmetop_UKB_het_fixed.txt.gz",
  "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_sumstats/iPSYCH-PGC_ASD_Nov2017.gz",
  "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_sumstats/ocd_aug2017.gz",
  "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_sumstats/pgc_alcdep.eur_unrelated.aug2018_release_refmt.txt.gz",
  "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_sumstats/pgc_bip_2018.gz",
  "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_sumstats/PGC_UKB_depression_genome-wide_fixed.txt.gz",
  "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_sumstats/MDD2018_ex23andMe.gz",
  "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_sumstats/PGC3_SCZ_wave3_public.v2.tsv.gz",
  "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_sumstats/tag.cpd.tbl.gz",
  "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_sumstats/TS_Oct2018.gz",
  "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_sumstats/pts_all_freeze2_overall.results.gz"
)
sample_sizes <- c(55374, 17310, 17310, 72517, 141932, 184765, 46351, 9725, 38686, 198882, 500199, 173005, 161405, 74053, 14307, 64357)
trait_names <- c("ADHD", "ANXFS", "ANX", "AN", "AUDIT", "CUD", "ASD", "OCD", "ALC", "BIP", "MDD", "MDD2018", "SCZ", "TAG_CPD", "TS", "PTSD")
sample_prev <- c(0.36, NA, 0.33, 0.23, NA, 0.30, 0.40, 0.28, 0.26, 0.15, 0.34, 0.35, 0.42, NA, 0.34, 0.15) # Sample prevalence
population_prev <- c(0.094, NA, 0.19, 0.006, NA, 0.16, 0.019, 0.012, 0.053, 0.028, 0.071, 0.071, 0.0075, NA, 0.0089, 0.068) # Population prevalence

traits <- paste0("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/delahaye_mQTL_analyses/", trait_names, ".sumstats.gz")
LDSCoutput <- readRDS("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/LDSC_covariances_EUR.RDS")
s_ld_dir <- "/scratch/st-dennisjk-1/wcasazza/1000G_phase3_ldsc/single_delahaye_annotations"
baseline <- "/scratch/st-dennisjk-1/wcasazza/1000G_v2.2_baseline/baselineLD"
frq <- "/scratch/st-dennisjk-1/wcasazza/1000G_Phase3_frq/1000G.EUR.QC"
wld <- "/scratch/st-dennisjk-1/wcasazza/weights_hm3_no_hla/weights"

annotation <- commandArgs(trailingOnly = TRUE)[[1]]
s_ld <- glue("{s_ld_dir}/{annotation}")
result <- s_ldsc(
  traits = traits,
  sample.prev = sample_prev,
  population.prev = population_prev,
  ld = c(baseline, s_ld),
  wld = wld,
  frq = frq,
  trait.names = trait_names
)
saveRDS(result, glue("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/SLDSCoutput_{annotation}.RDS"))
