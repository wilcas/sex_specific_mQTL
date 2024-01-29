[![DOI](https://zenodo.org/badge/337494687.svg)](https://zenodo.org/badge/latestdoi/337494687)

# Sex-dependent placental mQTL provide insight into the prenatal origins of childhood-onset traits and conditions

Code and notebooks corresponding to analysis done in the above titled publication [^1]. mQTL summary statistics and related annotations can be found on this [Open Science Framework site](https://osf.io/9r4wf/). 

## Directory structure
- [NICHD quality control for DNA methylation and genotyping](./NICHD_all_QC)
- [RICHS quality control for DNA methylation and genotyping](./richs_all_QC/)
- [Analysis code for publication](./mQTL_analyses/)
    - [main results figures and data processing](./mQTL_analyses/main_analyses/)
    - [Supporting R scripts](./mQTL_analyses/r_scripts/)
    - [Supporting TORQUE/MOAB scripts, written in bash](./mQTL_analyses/pbs_scripts/)

## Environment file
All packages used in this project are included in the `environment.yml` file. They can be installed [using anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)



[^1]: William Casazza, Amy M. Inkster, Giulia F. Del Gobbo, Victor Yuan, Fabien Delahaye, Carmen Marsit, Yongjin P. Park, Wendy P. Robinson, Sara Mostafavi, Jessica K. Dennis. Sex-dependent placental mQTL provide insight into the prenatal origins of childhood onset traits and conditions. *iScience*, 2024. https://doi.org/10.1016/j.isci.2024.109047
