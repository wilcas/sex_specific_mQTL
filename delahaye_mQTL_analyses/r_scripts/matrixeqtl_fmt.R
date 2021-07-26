library(data.table)

argv <- commandArgs(trailingOnly = TRUE)

df <- fread(argv[[1]])

fwrite(df[, .(SNP, gene = Gene, beta = b, `t-stat` = b / SE, `p-value` = p, FDR = p.adjust(p, method = "BH"))], argv[[2]], sep = "\t", quote = F, row.names = F)
