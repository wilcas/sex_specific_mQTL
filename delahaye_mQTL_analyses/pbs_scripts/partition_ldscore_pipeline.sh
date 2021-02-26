#!/bin/bash

base_annot_job=$(qsub make_annot.pbs)
base_annot_merge=$(qsub -W depend=afterok:${base_annot_job} run_baseline_merge.pbs)
qsub -W depend=afterok:${base_annot_merge} get_ldsc_baseline.pbs

interaction_mqtl="/scratch/st-dennisjk-1/wcasazza/delahaye_QC/matrix_eqtl_data/cis_mQTL_9_methy_PC_all_sex_interaction.txt"
marginal_mqtl="/scratch/st-dennisjk-1/wcasazza/delahaye_QC/matrix_eqtl_data/cis_mQTL_9_methy_PC_all.txt"

interaction_job=$(qsub -v MQTL="${interaction_mqtl}",Z_DIR="/scratch/st-dennisjk-1/wcasazza/delahaye_QC/caviar_files/sex_interaction",FILE_PREFIX="/scratch/st-dennisjk-1/wcasazza/1000G_phase3_ldsc/single_delahaye_annotations/sex_interaction" run_caviar_merge.pbs)
marginal_job=$(qsub -v MQTL="${marginal_mqtl}",Z_DIR="/scratch/st-dennisjk-1/wcasazza/delahaye_QC/caviar_files/marginal",FILE_PREFIX="/scratch/st-dennisjk-1/wcasazza/1000G_phase3_ldsc/single_delahaye_annotations/marginal" run_caviar_merge.pbs)

qsub -W depend=afterok:${interaction_job}:${marginal_job} get_ldsc.pbs
