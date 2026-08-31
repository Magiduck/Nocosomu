#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(magrittr)
library(ggplot2)
library(splitstackshape)

# make plots showcasing the frequencies of ICGC mutations per tumor type
# How many mutations per type

parquet_file_mutations <- "/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/6000_upstream_300_downstream_regions_of_protein_coding_genes/6000_upstream_300_downstream_regions_of_protein_coding_genes.parquet"
parquet_file_mutations <- "/groups/umcg-fg/tmp01/projects/non-coding-somatic/ICGC_per_tumor_type/data_views/columnar_view/all_simple_somatic_mutations_ICGC_release28_unsorted.parquet"

icgc_subgroups = fread("/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/scratch_cgut/nocosomu_iDriver/icgc_subgroups.tsv.gz")

#get real mutation data
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
project_gene_mutations <- dbGetQuery(con, glue("SELECT DISTINCT icgc_donor_id, mutation_type, project_code, COUNT(DISTINCT icgc_mutation_id)
                                                 FROM '{parquet_file_mutations}'
                                                 WHERE sequencing_strategy = 'WGS'
                                                 GROUP BY icgc_donor_id, mutation_type, project_code")) %>% as.data.table
setnames(project_gene_mutations, "count(DISTINCT icgc_mutation_id)", "nevents")
project_gene_mutations[icgc_subgroups, on=.(icgc_donor_id == icgc_donor_id), subgroup := my_cancer_type]
project_gene_mutations[, subgroup_f := factor(subgroup, levels = project_gene_mutations[, median(nevents), by = subgroup][order(V1), subgroup])]

xlabs <- paste(levels(project_gene_mutations$subgroup_f)," (n=",table(unique(project_gene_mutations[,.(subgroup_f, icgc_donor_id)]) %>% .[, subgroup_f]),")",sep="")

ggplot(project_gene_mutations, aes(subgroup_f, log10(nevents))) +
  geom_boxplot(aes(fill = mutation_type, col = "#FFFFFF"), position = "dodge", outlier.shape = 4, outlier.size = 1, outlier.alpha = 0.1, linewidth = 0.1) +
  scale_x_discrete(labels = xlabs) +
  scale_color_manual(values = "#555555") +
  theme_bw() +
  coord_flip()
