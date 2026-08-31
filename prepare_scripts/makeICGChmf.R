#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)
library(ggplot2)

alleles = c("A","C","G","T")

all_icgc = list.files(path = "/groups/umcg-fg/tmp01/projects/non-coding-somatic/cancer/icgc_results/frequency_matrices/downloads/ICGC/", 
                         pattern = "*.purple.somatic.vcf.gz$", 
                         recursive = TRUE,
                         full.names = T)

for (my_file in all_icgc){
  print(my_file)
my_file = "/groups/umcg-fg/tmp01/projects/non-coding-somatic/cancer/icgc_results/frequency_matrices/downloads/ICGC//PACA-AU/DO49178/purple/DO49178T.purple.somatic.vcf.gz"
icgc_donor_id = strsplit(my_file, "/")[[1]][14]
project_code = strsplit(my_file, "/")[[1]][13]
vcf_file = fread(cmd = glue("zcat {my_file}")) %>% .[FILTER == "PASS" & REF %in% alleles & ALT %in% alleles]

vcf_file[, mutation_id := paste0(`#CHROM`, "_", POS, "_", REF, "_", ALT)]
vcf_file[, icgc_donor_id := icgc_donor_id]
vcf_file[, project_code := project_code]


for (chromosome in vcf_file[, unique(`#CHROM`)]){
 fwrite(vcf_file[`#CHROM` == chromosome, .(mutation_id, icgc_donor_id, project_code, `#CHROM`, POS, REF, ALT)], 
        glue("/groups/umcg-fg/tmp01/projects/non-coding-somatic/cancer/icgc_results/frequency_matrices/artifacts/icgcHMF/{chromosome}_ssm.tsv"),
        append = my_file != all_icgc[1]) 
}
}
}

all_ssm = list.files(path = "/groups/umcg-fg/tmp01/projects/non-coding-somatic/cancer/icgc_results/frequency_matrices/artifacts/icgcHMF/", 
                      pattern = "*_ssm.tsv$", 
                      full.names = T)
icgc_subgroups = fread("/groups/umcg-fg/tmp01/projects/non-coding-somatic/ICGC_per_tumor_type/icgc_subgroups.tsv")

for (my_file in all_ssm){
  my_ssm = fread(my_file)
  
  my_ssm[icgc_subgroups, on=.(icgc_donor_id=icgc_donor_id), gel_cancer := my_cancer_type]
  setnames(my_ssm, "icgc_donor_id", "donor_id")
  
  fwrite(my_ssm[, .(`#CHROM`,	POS,	REF,	ALT,	gel_cancer,	donor_id)],
         "/groups/umcg-fg/tmp01/projects/non-coding-somatic/cancer/icgc_results/frequency_matrices/artifacts/icgcHMF/icgcHMF_all_donorPurplePASS.tsv", 
         sep = "\t",
         append = my_file != all_ssm[1])
}
