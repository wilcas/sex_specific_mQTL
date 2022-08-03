## ----setup, echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(minfi)
library(wateRmelon)
library(data.table)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
metadata <- read.delim("/scratch/st-dennisjk-1/wcasazza/delahaye_QC/placenta_regulatory_landscape/RootStudyConsentSet_phs001717.PlacentalRegulation.v1.p1.c1.HMB-IRB-PUB-COL-MDS/PhenotypeFiles/Placental_Regulation_Phenotypes_clean.txt")
metadata$sample <- as.character(metadata$SUBJECT_ID)


## ----read GEO data and separate meth unmeth objects---------------------------------------------------------------------------------------------------------------------------------------------------

# Magda's code

# read in the 1st 5 lines of the file
methy_data1 <- fread("/scratch/st-dennisjk-1/wcasazza/delahaye_QC/placenta_regulatory_landscape/RootStudyConsentSet_phs001717.PlacentalRegulation.v1.p1.c1.HMB-IRB-PUB-COL-MDS/GenotypeFiles/phg001182.v1.PlacentalRegulation.genotype-original-submission.divers_arrays.c1.HMB-IRB-PUB-COL-MDS.update/corrected_run44_NICHD_ctrl-bkg.txt")
methy_data2 <- fread("/scratch/st-dennisjk-1/wcasazza/delahaye_QC/placenta_regulatory_landscape/RootStudyConsentSet_phs001717.PlacentalRegulation.v1.p1.c1.HMB-IRB-PUB-COL-MDS/GenotypeFiles/phg001182.v1.PlacentalRegulation.genotype-original-submission.divers_arrays.c1.HMB-IRB-PUB-COL-MDS.update/corrected_run47_NICHD_ctrl_bkg.txt")[, -c(1)]
methy_data <- cbind(methy_data1, methy_data2)
run_44 <- unique(gsub("\\..*", "", colnames(methy_data1)))
run_47 <- unique(gsub("\\..*", "", colnames(methy_data2)))
metadata$batch <- ifelse(metadata$SUBJECT_ID %in% run_44, "run_44", ifelse(metadata$SUBJECT_ID %in% run_47, "run_47", NA))
colnames(methy_data) # figure out naming convention, should contain columns for methylated, unmethylated, and detection p values for all samples.
methy_data <- as.data.frame(methy_data)
# separate the methy_data data into 3 matrices: Detection.Pval, Unmethylated.Signal and Methylated.Signal (these can sometimes be named A/B/Green/Red)
rownames(methy_data) <- methy_data[, 1] # first column contains probe IDs, make  rownames
methy_data <- methy_data[, -1]

Detection.1 <- methy_data[, grep(".Detection Pval", colnames(methy_data))]
colnames(Detection.1) <- gsub(".Detection Pval", "", colnames(Detection.1))
colnames(Detection.1) # 36

unmeth.1 <- methy_data[, grep(".Signal_A", colnames(methy_data))]
colnames(unmeth.1) <- gsub(".Signal_A", "", colnames(unmeth.1))
colnames(unmeth.1) # 36

meth.1 <- methy_data[, grep(".Signal_B", colnames(methy_data))]
colnames(meth.1) <- gsub(".Signal_B", "", colnames(meth.1))
colnames(meth.1) # 36

# sanity check that everything matches
identical(rownames(meth.1), rownames(unmeth.1))
identical(rownames(meth.1), rownames(Detection.1))

identical(colnames(meth.1), colnames(unmeth.1))
identical(colnames(meth.1), colnames(Detection.1))



## ----create mset--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

mset <- new("MethylSet",
  SummarizedExperiment(assays = SimpleList(
    Meth = as.matrix(meth.1),
    Unmeth = as.matrix(unmeth.1)
  )),
  annotation = c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"),
  preprocessMethod = c(
    rg.norm = "Raw (no normalization or bg correction)",
    minfi = "1.24.0", manifest = "0.4.0"
  )
)
# pData(mset) <- DataFrame(pData1) # add in pData
mset # this object should be 485577 rows long (if it doesn't contain SNP probes will be 485512)



## ----check mset---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# code from Victor

# this code checks that the methylset has been constructed from the correct meth/unmeth columns by pulling the beta values and making sure that the beta value density distribution makes sense. the hypomethylated peak should be bigger if the mset has been constructed properly.

grset <- mapToGenome(ratioConvert(mset))

# plot density, if the meth/unmeth intensities are flipped, then the hypermethylated peak will be bigger than the hypomethylated peak.
betas <- as.data.frame(getBeta(mset))
head(betas)
ggplot(betas, aes(x = `8608`)) +
  geom_density()



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
betas <- getBeta(mset)
dim(betas) # 485512     36

head(colnames(betas)) # sentrix ID and position concatenated
colnames(betas) <- make.unique(metadata$sample)

head(colnames(mset))
colnames(mset) <- make.unique(metadata$sample)



## ----check sex----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# extract the probes from the Illumina manifest corresponding to the X/Y chromosome
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Locations)
sexProbeIDs <- Locations[which(Locations$chr == "chrY" | Locations$chr == "chrX"), ]
yProbes <- Locations[which(Locations$chr == "chrY"), ]
# pull XY probes
betasXY <- betas[which(rownames(betas) %in% rownames(sexProbeIDs)), ] # 11648
dim(betasXY) # 11648 by 46 good
sum(is.na(betasXY)) # 0

# pull xist probes
betasXIST <- betas[rownames(betas) %in% c("cg03554089", "cg12653510", "cg05533223", "cg11717280", "cg20698282"), ]
dim(betasXIST)

# color by metadata sex
sexLabel <- as.vector(metadata$Gender)
table(sexLabel) # 23 males and 23 females
sexLabel <- gsub("Female", "#CC79A7", sexLabel)
sexLabel <- gsub("Male", "#0072B2", sexLabel)

batchLabel <- palette(rainbow(10))[as.factor(metadata$batch)]
# cluster on XY probes - some mislabeled samples
plot(ape::as.phylo(hclust(dist(t(betasXY)))), lab4ut = "axial", type = "unrooted", edge.width = 0.5, cex = 0.8, tip.color = sexLabel)
plot(ape::as.phylo(hclust(dist(t(betasXY)))), lab4ut = "axial", type = "unrooted", edge.width = 0.5, cex = 0.8, tip.color = batchLabel)
# cluster on XIST probes - some mislabelled samples
plot(ape::as.phylo(hclust(dist(t(betasXIST)))), lab4ut = "axial", type = "unrooted", edge.width = 0.5, cex = 0.8, tip.color = sexLabel)
plot(ape::as.phylo(hclust(dist(t(betasXIST)))), lab4ut = "axial", type = "unrooted", edge.width = 0.5, cex = 0.8, tip.color = batchLabel)

# just looking at batch
# plot(ape::as.phylo(hclust(dist(t(betas)))), lab4ut="axial", type = "unrooted", edge.width=0.5, cex=0.8, tip.color=batchLabel)



## ----check sample quality-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

qc <- getQC(mset)
head(qc)

plotQC(qc) # all samples are good quality (yay) none need to be removed

# add QC to your metadata (can join pData and qc) if any samples look odd at this step, decide on a threshold for m,edian intensity to retain the majority of the samples.

mset_filt <- mset
betas <- getBeta(mset_filt)



## ----poor quality probes------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

detP <- methy_data[, grepl("Detection Pval", colnames(methy_data))]
colnames(detP) <- gsub("\\..*", "", colnames(detP))
head(detP)

failed <- detP > 0.01
colMeans(failed) # average # failed probes per sample
sum(rowMeans(failed) > 0.2) # 323 failed probes (detP>0.01) in >20% samps
failed <- rownames(detP[rowMeans(failed) > 0.2, ]) # vector of 350 failed probe IDs
head(failed)

# remove detP failing probes
dim(betas) # 485512
betas <- betas[!(rownames(betas) %in% failed), ]
dim(betas) # 485189


naBeta <- is.na(betas)
naBeta <- print(sum(rowMeans(naBeta) > 0.2)) # 79 missing betas in >20% samps

# remove NA  probes
dim(betas) # 485189
betas <- betas[!(rownames(betas) %in% naBeta), ]
dim(betas) # 485110



## ----SNP and XY probes--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# SNP probes - already gone (by preprocessNoob)
rs_probes <- grep(pattern = "^rs", x = rownames(betas)) # vector of rs rownames
length(rs_probes) # 0
# XY probes
dim(betas) # 342762
betasXY <- betas[(rownames(betas) %in% rownames(sexProbeIDs)), ]
dim(betasXY) # 9813 preprocessed X probes

betas <- betas[!(rownames(betas) %in% rownames(yProbes)) & !(rownames(betas) %in% rs_probes), ]
dim(betas) # 332949 preprocessed autosomal probes



## ----cross hybridizing--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
priceAnno <- read.csv("/arc/project/st-dennisjk-1/wcasazza/methylation_QC_resources/price_et_al_450k_annot.csv")
priceXY <- priceAnno[grepl("YES", priceAnno$XY_Hits), ] # 12388
priceAuto <- priceAnno[grepl("YES", priceAnno$Autosomal_Hits), ] # 40650

# remove these CH probes
dim(betas) # 485110
betas <- betas[!(rownames(betas) %in% priceXY$ID), ]
betas <- betas[!(rownames(betas) %in% priceAuto$ID), ]
dim(betas) # 443324


## ----polymorphic probes, eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------
##
## # read in chen polymoprhic probes
## Chen <- readxl::read_xlsx("/arc/project/st-dennisjk-1/wcasazza/methylation_QC_resources/48640-polymorphic-CpGs-Illumina450k.xlsx")
##
## # Polymorphic CpGs (including polymorphic SBE sites) correspond to entries where BASE_FROM_SBE is less than or equal to 1 for infinium II probes and less than or equal to 2 for infinium I probes.
## str(Chen)
## Chen$BASE_FROM_SBE <- as.factor(Chen$BASE_FROM_SBE)
## levels(Chen$BASE_FROM_SBE) # 0, 1, and 2 probes are either polymorphic or have poly SBE
##
## # check that none of the "BASE_FROM_SBE" values are >1 for tye II probes
## ProbeII <- Chen[which(Chen$INFINIUM == "II"),]
## dim(ProbeII) # 57536 rows
## ProbeII$BASE_FROM_SBE <- factor(ProbeII$BASE_FROM_SBE)
## levels(ProbeII$BASE_FROM_SBE) # 0 and 1
##
## # select only those probes that have a MAF > 1%
## dim(Chen) # 70889
## Chen_SNP <- Chen[which(Chen$AF > 0.01),]
## dim(Chen_SNP) # 14140
##
## # remove all of the probes identified as polymorphic by Chen et al.
## dim(betas) # 443324
## betas <- betas[!(rownames(betas) %in% Chen_SNP$PROBE),]
## dim(betas) # 431409
##


## ----subset mset--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# define list of filtered probes
probes_keep <- as.character(rownames(betas))

mset_filt <- mset[probes_keep, -c(76, 89)]
dim(mset_filt)

all.equal(rownames(betas), rownames(mset_filt)) # should be TRUE



## ----preprocessquantile-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

set.seed(42)
grset.quantile <- preprocessQuantile(mset_filt, sex = na.omit(recode(as.character(metadata$Gender[match(colnames(mset_filt), metadata$SUBJECT_ID)]), "Male" = "M", "Female" = "F")))
betas.quantile <- getBeta(grset.quantile)
probeDesign <- data.frame(Type = getProbeType(grset.quantile, withColor = FALSE))
probeDesign$Name <- rownames(grset.quantile)
plotBetasByType(betas.quantile[, 1], probeTypes = probeDesign, main = "Quantile")




## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
write.table(metadata, "DNAm_w_paired_genotype_metadata.txt", row.names = F, quote = F, sep = "\t")
fwrite(as.data.frame(betas.quantile) %>% rownames_to_column("cpg"), "processed_DNAm_delahaye_quantile_norm.txt", quote = F, row.names = F, sep = "\t")
