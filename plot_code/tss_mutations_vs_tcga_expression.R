# How is the expression of these genes

tcga <- fread("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/TCGA_AfterQC/output/TCGA_ALL_TPM_log2_QNorm_QCed_CovCorrected_AllCovariates_new.txt")
tcga_annotation <- fread("/groups/umcg-fg/tmp01/projects/genenetwork/recount3/tissuePredictions/samplesWithPrediction_16_09_22_noOutliers.txt")

significant_tss_mutations = fread("/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/6000_upstream_300_downstream_regions_of_protein_coding_genes/all_simple_somatic_mutations_ICGC_release28_chrsorted_WGS/SuRE_ResNet_Attention_SuRE4n_tss_selection_m300_p100_stranded_EnhA_intersection_intersection_K562_TSS_EnhA_strand_False_trial_20_model_CNN/plots/ISM_abs/donors_with_tss_mutations.tsv")
setkey(significant_tss_mutations, chromosome, chromosome_start, chromosome_end)
significant_tss_mutations[, chromosome_str := as.character(chromosome)]
significant_tss_donors = paste0("('", significant_tss_mutations[,paste0(unique(icgc_donor_id), collapse = "','")], "')")


significant_gene_mask <- tcga$V1 %in% unique(significant_tss_mutations$associated_gene)
significant_gene_mask <- tcga$V1 %in% sample(tcga$V1, 66, replace = F)

tcga_result <- tcga[,1:20][ , sapply(.SD, function(x){wilcox_result <- wilcox.test(x[significant_gene_mask],x[!significant_gene_mask], alternative = "two.sided"); return(wilcox_result$p.value)}), .SDcols = !"V1"]

reshape2::melt(tcga_result)
