# Quality control (QC) for NICHD DNA methylation and genotyping data
## Executing this code

QC code for genotyping data was written using [Snakemake](https://snakemake.readthedocs.io/en/stable/executing/cli.html), and used according to the pipeline outlined in [`delahaye_genotyping_qc.nb.html`](./delahaye_genotyping_qc.nb.html). QC code can be found in [`Snakefile`](./Snakefile), and additional scripts called by this program are available in [`scripts`](./scripts/).

Code for DNAm QC was written in R, primarily making use of [`minfi`](https://doi.org/doi:10.18129/B9.bioc.minfi). It is stored in [`delahaye_dnam_qc.nb.html`](./delahaye_dnam_qc.nb.html).

## Viewing analysis

Analysis is present in [R-notebook files](https://bookdown.org/yihui/rmarkdown/notebook.html) (with the extension `.nb.html`) and can be downloaded, opened, and viewed in a web browser.
