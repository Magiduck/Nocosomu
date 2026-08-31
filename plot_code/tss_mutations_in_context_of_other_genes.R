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

parquet_dir = "/groups/umcg-fg/tmp01/projects/non-coding-somatic/ICGC_per_tumor_type/data_views/columnar_view/"

parquet_consequence <- glue("{parquet_dir}/all_simple_somatic_mutations_ICGC_release28_unsorted.parquet")
parquet_regions <- "/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/6000_upstream_300_downstream_regions_of_protein_coding_genes/6000_upstream_300_downstream_regions_of_protein_coding_genes_OLD.parquet"

significant_tss_mutations = unique(fread("/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/6000_upstream_300_downstream_regions_of_protein_coding_genes/all_simple_somatic_mutations_ICGC_release28_chrsorted_WGS/SuRE_ResNet_Attention_SuRE4n_tss_selection_m300_p100_stranded_EnhA_intersection_intersection_K562_TSS_EnhA_strand_False_trial_20_model_CNN/plots/ISM_abs/donors_with_tss_mutations.tsv"))
setkey(significant_tss_mutations, chromosome, chromosome_start, chromosome_end)
significant_tss_mutations[, chromosome_str := as.character(chromosome)]
significant_tss_donors = paste0("('", significant_tss_mutations[,paste0(unique(icgc_donor_id), collapse = "','")], "')")

#Are other genes recurrently KO? (Support for)
non_coding <- "('intergenic_region','intron_variant','upstream_gene_variant','exon_variant','downstream_gene_variant','synonymous_variant')"
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
gene_consequence_donors <- current_mutations <- dbGetQuery(con, glue("SELECT DISTINCT icgc_donor_id,
                                                 project_code, 
                                                 consequence_type,
                                                 COUNT(DISTINCT icgc_mutation_id),
                                                 aa_mutation,
                                                 cds_mutation,
                                                 gene_affected,
                                                 transcript_affected
                                                 FROM '{parquet_consequence}'
                                                 WHERE consequence_type NOT IN {non_coding}
                                                 GROUP BY icgc_donor_id, project_code, consequence_type, gene_affected, transcript_affected, aa_mutation, cds_mutation")) %>% as.data.table
dbDisconnect(con)


ptv <- c("stop_gained",
        "frameshift_variant",
        "disruptive_inframe_insertion",
        "disruptive_inframe_deletion",
        "start_lost")

significant_tss_mutations[,.N,by=.(associated_gene)][order(N)]

gene_consequence_donors[consequence_type %in% ptv, .N, by = icgc_donor_id]

gene_ptv_donors <- gene_consequence_donors[consequence_type %in% ptv,]

patients_with_tert <- significant_tss_mutations[associated_gene == "ENSG00000164362", unique(icgc_donor_id)]
projects_with_tert <- gene_consequence_donors[icgc_donor_id %in% patients_with_tert, unique(project_code)]
patients_in_same_projects <- gene_consequence_donors[project_code %in% projects_with_tert & !(icgc_donor_id %in% patients_with_tert),unique(icgc_donor_id)]


perm_results <- lapply(1:1000, function(x){
  perm_measurements = gene_ptv_donors[icgc_donor_id %in% sample(patients_in_same_projects, 87, replace = F), .N, by = gene_affected]
  perm_measurements[, n_fold := x]
  return(perm_measurements)
}) %>% do.call("rbind", .)

real_results = gene_ptv_donors[icgc_donor_id %in% patients_with_tert, .N, by = gene_affected]
real_results[, n_fold := "real"]

max_perm_results <- perm_results[,.SD[order(-N)][1],by=.(gene_affected)]

merged_real_perm = max_perm_results[real_results, on=.(gene_affected == gene_affected)]
merged_real_perm[is.na(N), N:= 0]
merged_real_perm[i.N > N][order(i.N)]

perm_plot_data <- rbind(real_results, perm_results)


#Many patients with multiple significant genes
significant_tss_mutations[,.N,by=icgc_donor_id][,hist(N)]

#tert then other genes with many patients
#significant_tss_mutations[,.N,by=associated_gene][,hist(N)]


