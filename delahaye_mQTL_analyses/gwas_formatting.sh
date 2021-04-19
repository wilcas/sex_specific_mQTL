#!/bin/bash
source /home/wcasazza/miniconda3/bin/activate
conda activate ldsc
cd /scratch/st-dennisjk-1/wcasazza/tmp_GWAS/
echo "" > missed_files.txt
for f in pgc_sumstats/*.gz; do
    out_f=${f%.*}
    skip=""
    if [[ ! -f "pgc_formatted_sumstats/${out_f##*/}.sumstats.gz" ]]; then
        cmd="/arc/project/st-dennisjk-1/software/ldsc/munge_sumstats.py"
        cmd+=" --sumstats ${f}"
        cmd+=" --out pgc_formatted_sumstats/${out_f##*/}"
        for arg in $(zgrep -v "^#" ${f} | head -n 1 | tr -d '\r\n'); do
            case ${arg} in
                "TotalN")
                    cmd+=" --N-col TotalN"
                ;;
                "totalN")
                    cmd+=" --N-col totalN"
                ;;
                "Neff")
                    cmd+=" --N-col Neff"
                ;;
                "Nca")
                    cmd+=" --N-cas-col Nca"
                ;;
                "NCAS")
                    cmd+=" --N-cas-col NCAS"
                ;;
                "NCON")
                    cmd+=" --N-con-col NCON"
                ;;
                "Nco")
                    cmd+=" --N-con-col Nco"
                ;;
                "Z")
                    cmd+=" --signed-sumstats Z,0"
                    break
                ;;
                "BETA")
                    cmd+=" --signed-sumstats BETA,0"
                ;;
                "LogOR")
                    cmd+=" --signed-sumstats LogOR,0"
                ;;
                "Log10BF")
                    cmd+=" --signed-sumstats Log10BF,0"
                ;;
                "REF")
                    cmd+=" --a1 REF"
                ;;
                "ALT")
                    cmd+=" --a2 ALT"
                ;;
                "ID")
                    cmd+=" --snp ID"
                ;;
                *)
                    echo "${arg} ignored."
                ;;
            esac
        done
            # Datasets without N column and other special cases
            case "${out_f##*/}" in
                "adhd_jul2017")
                    cmd+=" --N-cas 20183"
                    cmd+=" --N-con 35191"
                ;;
                "AUDIT_UKB_2018_AJP.txt") #Focus on TOTAL Audit score only
                    cmd+=" --signed-sumstats beta_T,0"
                    cmd+=" --p p_T"
                    cmd+=" --a1 a_0"
                    cmd+=" --a2 a_1"
                ;;
                "iPSYCH-PGC_ASD_Nov2017")
                    cmd+=" --N-cas 18381"
                    cmd+=" --N-con 27969"
                ;;
                "ocd_aug2017")
                    cmd+=" --N-cas 2688"
                    cmd+=" --N-con 7037"
                ;;
                "pgc_adhd_females")
                    cmd+=" --N-cas 4945"
                    cmd+=" --N-con 16246"
                ;;
                "pgc_adhd_males")
                    cmd+=" --N-cas 14154"
                    cmd+=" --N-con 17948"
                ;;
                "pgc_alcdep.trans_re2_unrel_geno.aug2018_release.txt")
                    skip=1
                ;;
                "pgc_alcdep.trans_mantra_unrel_geno.aug2018_release.txt")
                    skip=1
                ;;
                "pgc.cross.full.2013-03.txt")
                    cmd+=" --N-cas 33332"
                    cmd+=" --N-con 27888"
                ;;
                "PGC_UKB_depression_genome-wide.txt")
                    cmd+=" --N-cas 246363"
                    cmd+=" --N-con 561190"
                ;;
                "tag.cpd.tbl")
                    cmd+=" --N 74053"
                    cmd+=" --signed-sumstats OR,0"
                ;;
                "TS_Oct2018")
                    cmd+=" --N-cas 4819"
                    cmd+=" --N-con 9788"
                ;;
                *)
                    echo "Header of ${f} was parsed."
                ;;
            esac
        if [[ ! ${skip} ]]; then
            eval $cmd
        fi
    fi
done
gwas_out=""
for f in /scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/*.gz; do
    tmp="${f/.sumstats.gz/}"
    tmp="${tmp##*/}"
    gwas_out+=" $tmp"
done
echo $gwas_out
ls /scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/*.gz | wc -l
