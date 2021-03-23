#!/bin/bash
z_score_job=$(qsub -v MQTLS="/scratch/st-dennisjk-1/wcasazza/ariesmqtl/cord.ALL.M.tab",Z_DIR="cord_mqtl" run_caviar_per_chromosome.pbs)
