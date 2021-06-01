#!/bin/bash
ANNOT=( "male_meta" "female_meta" "marginal_meta" "sex_interaction_meta" )

for cur_annot in ${ANNOT[@]}; do 
  cur_mqtl="/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/delahaye_richs_${cur_annot}_matrixeqtl.txt" 
  mqtl_job=$(qsub -v MQTLS="${cur_mqtl}",Z_DIR="/tmp/${cur_annot}",FILE_PREFIX="/scratch/st-dennisjk-1/wcasazza/1000G_phase3_ldsc/single_delahaye_annotations/${cur_annot}" run_caviar_per_chromosome.pbs)

  qsub -W depend=afterok:${mqtl_job} -v ANNOTATIONS="${cur_annot}" get_ldsc.pbs

done
