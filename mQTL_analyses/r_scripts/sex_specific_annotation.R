library(data.table)
library(glue)
library(parallel)

merge_chromosome <- function(f1, f2, out_name){
  t1 <- fread(f1)
  t2 <- fread(f2)
  col_val <- colnames(t1)[5]
  if(col_val == "maxCPP"){
    t1[, c(col_val) := {
        tmp <- pmax(t1[[col_val]], t2[[col_val]])
        tmp[!(t1[[..col_val]] & t2[[..col_val]])] <- 0
        tmp
      }
    ]
    res <- t1
  }else{
    t1[, c(col_val) := as.numeric(t1[[..col_val]] & t2[[..col_val]])]
    res <- t1
  }
  fwrite(res, out_name, sep='\t', row.names=F, quote=F)
}

for(sex in c("male", "female")){
  for(annotation in c("maxCPP","bin_CPP","all_cis")){
    mclapply(
      1:22,
      function(i){
        merge_chromosome(
          glue("/scratch/st-dennisjk-1/wcasazza/1000G_phase3_ldsc/single_delahaye_annotations/sex_interaction_meta_{annotation}.{i}.annot.gz"),
          glue("/scratch/st-dennisjk-1/wcasazza/1000G_phase3_ldsc/single_delahaye_annotations/{sex}_meta_{annotation}.{i}.annot.gz"),
          glue("/scratch/st-dennisjk-1/wcasazza/1000G_phase3_ldsc/single_delahaye_annotations/{sex}_specific_meta_{annotation}.{i}.annot.gz")
        )
      },
      mc.cores=16
    )
  }
}
