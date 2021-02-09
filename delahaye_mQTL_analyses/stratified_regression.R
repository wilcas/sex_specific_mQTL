library(MASS)
library(sfsmisc)
library(data.table)
library(parallel)
argv <- commandArgs(trailingOnly=T)
run_rlm <- function(i, male_selected, female_selected){
  row <- sig_sex_fdr[i,]
  X <-  as.numeric(genotype[row$SNP,-c(1),on="SNP"])
  y <- as.numeric(methylation[row$gene,-c(1),on="cpg"])
  is_male <- pheno[nrow(pheno),-c(1)] == 1
  male_snp <- X[male_selected]
  female_snp <- X[female_selected]
  covar <- t(pheno[-nrow(pheno),-c(1)])
  X <- cbind(covar[,-c(nrow(pheno))],X)
  if(length(unique(male_snp)) != 1){
    male_coef<- summary(lm(y[male_selected & is.finite(y) & !is.na(y)] ~X[male_selected & is.finite(y) & !is.na(y),]))$coef
    simple_fit <- lm(y[male_selected & is.finite(y) & !is.na(y)] ~X[male_selected & is.finite(y) & !is.na(y),ncol(X)])
    beta <- summary(simple_fit)$coef[2,1]
    int <- summary(simple_fit)$coef[1,1]
    aFC <- sign(beta)*log2(abs(2*beta / int) + 1)
    res_m <- data.frame(res_m=male_coef[nrow(male_coef),4],beta_m=male_coef[nrow(male_coef),1],se_m=male_coef[nrow(male_coef),2],aFC_m=aFC)
  }else{
    res_m <- data.frame(res_m=NA,beta_m=NA,se_m=NA,aFC_m=NA)
  }
  if(length(unique(female_snp)) != 1){
    female_coef <- summary(lm(y[female_selected & is.finite(y)& !is.na(y)] ~X[female_selected & is.finite(y)& !is.na(y),]))$coef
    simple_fit <- lm(y[female_selected & is.finite(y)& !is.na(y)] ~X[female_selected & is.finite(y)& !is.na(y),ncol(X)])
    beta <- summary(simple_fit)$coef[2,1]
    int <- summary(simple_fit)$coef[1,1]
    aFC <- sign(beta)*log2(abs(2*beta / int) + 1)
    res_f<- data.frame(res_f=female_coef[nrow(female_coef),4],beta_f=female_coef[nrow(female_coef),1],se_f=female_coef[nrow(female_coef),2],aFC_f=aFC)
  }else{
    res_f<- data.frame(res_f=NA,beta_f=NA,se_f=NA,aFC_f=NA)
  }
  return(cbind(row,res_m,res_f))
}
sex_spec <-fread("matrix_eqtl_data/cis_mQTL_9_methy_PC_all_sex_interaction.txt")
genotype <- fread("matrix_eqtl_data/all_imputed_matrixeQTL.txt")
setindex(genotype, SNP)
methylation <- fread("matrix_eqtl_data/methylation_matrixeQTL.txt")
setindex(methylation,cpg)
pheno <- fread("matrix_eqtl_data/mQTL_covar_9_methy_PC.txt")
is_male <- pheno[nrow(pheno),-c(1)] == 1
male_idx <- which(is_male)[sample(1:sum(is_male),100)]
male_selected <- rep(FALSE,length(is_male))
male_selected[male_idx] <- TRUE
is_female <- pheno[nrow(pheno),-c(1)] == 0
female_idx <- which(!is_male)[sample(1:sum(is_female),100)]
female_selected <- rep(FALSE,length(is_female))
female_selected[female_idx] <- TRUE
sig_sex_fdr <- sex_spec[FDR < 0.05]
res <- rbindlist(lapply(1:nrow(sig_sex_fdr),function(i) run_rlm(i, male_selected,female_selected)))
fwrite(res,sprintf("delahaye_stratified_regression_fdr_hits_rep_%s.txt",argv[[1]]),quote=F,sep="\t",row.names=F)
