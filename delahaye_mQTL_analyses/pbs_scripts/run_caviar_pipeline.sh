#!/bin/bash
ANNOT=("male_meta" "female_meta" "marginal_meta" "sex_interaction_meta")

# "male_specific_meta" "female_specific_meta")
for cur_annot in ${ANNOT[@]}; do
	if [ $cur_annot = "cord_mqtl" ]; then
		cur_mqtl="/scratch/st-dennisjk-1/wcasazza/ariesmqtl/cord.ALL.M.tab"
	else
		cur_mqtl="/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/delahaye_richs_${cur_annot}_matrixeqtl_maf01.txt.gz"
	fi
	mqtl_job=$(qsub -v ANNOT_DIR="${cur_annot}",MQTLS="${cur_mqtl}",Z_DIR="/tmp/${cur_annot}",FILE_PREFIX="/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/${cur_annot}" run_caviar_detection_per_chromosome.pbs)

	#qsub -W depend=afterok:${mqtl_job} -v ANNOTATIONS="${cur_annot}" get_ldsc.pbs
	#qsub -v ANNOTATIONS="${cur_annot}" get_ldsc.pbs

done
