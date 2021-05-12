library(GenomicSEM)
library(data.table)
LDSCoutput <- readRDS("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/LDSC_covariances.RDS")
p_sumstats <- fread("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/delahaye_mQTL_analyses/common_factor_sumstats_15_traits.txt.gz")
pfactor <- commonfactorGWAS(covstruc = LDSCoutput, SNPs = p_sumstats)

fwrite(pfactor,"/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/common_factor_gwas_result.txt.gz",quote=F,row.names=F,sep='\t')
