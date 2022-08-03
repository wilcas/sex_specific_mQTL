library(data.table)
library(parallel)
annot_str <- commandArgs(trailingOnly = TRUE)[[1]]
annot <- fread(annot_str)
fwrite(annot[SNP != "."], annot_str, sep = "\t", quote = F, row.names = F)
