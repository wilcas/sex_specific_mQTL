library(data.table)
library(glue)
library(parallel)

argv <- commandArgs(trailingOnly = TRUE)
# 1: CpGs (data.table(probe))
# 2: Caviar results (data.table(
# 3: fmt for output

cpgs <- fread(argv[[1]])
caviar <- fread(argv[[2]])
fmt <- argv[[3]]

caviar <- caviar[cpgs, .(maxCPP = max(CPP)), by = .(SNP), on = .(probe)]
# Below assumes baseline in hard-coded location
mclapply(
  1:22,
  function(i) {
    annot <- fread(glue("/scratch/st-dennisjk-1/wcasazza/1000G_v2.2_baseline/baselineLD.{i}.annot.gz"))
    out <- caviar[annot, .(CHR, BP, SNP, CM, maxCPP = nafill(maxCPP, fill = 0)), on = .(SNP)]
    fwrite(out, glue(fmt), sep = "\t", row.names = FALSE, quote = FALSE)
  },
  mc.cores = 8
)
