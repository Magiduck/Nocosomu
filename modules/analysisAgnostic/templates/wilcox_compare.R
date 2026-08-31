#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)
library(ggplot2)
library(dplyr)

#Prepare allele changes

##################

#borzoi prepare (harmonize prediction column)
#sure prepare (dummy gene)

current_chromosome = "$chromosome"
# 
# #get preliminary files
minimum_donors = $minimum_donors
initial_region_selection = fread("region_selection.tsv")

donorMatches = fread("no_background.tsv")
background_donors = fread("background.tsv")
qc_passing_background = fread("background_mutations.tsv")

#second unique performed since only real mutations once epr occurrence
raw_real_mutations = fread("real_mutations.tsv")
real_mutations = raw_real_mutations %>% unique %>% .[unique(donorMatches[,.(donor_id,subgroup, associated_region,real_mutation_id)]), on=.(ID==real_mutation_id), allow.cartesian=T]
real_mutations[, label:="observed_mutation"]

#no unique performed since backgrounds can clash among mutations
background_mutations = donorMatches[,.(donor_id,subgroup, associated_region,matching_background_mutation_id)] %>% unique %>% .[qc_passing_background, on=.(matching_background_mutation_id=ID), nomatch=0, allow.cartesian = TRUE]
background_mutations[, label:="comparison_group"]
setnames(background_mutations, "matching_background_mutation_id", "ID")

#some background mutations were observed in some patients so they were calculated in realmutations
sometimes_observed_background_mutations = donorMatches[,.(donor_id,subgroup, associated_region,matching_background_mutation_id)] %>% unique %>% .[fread("real_mutations.tsv"), on=.(matching_background_mutation_id=ID), nomatch=0, allow.cartesian = TRUE]
sometimes_observed_background_mutations[, label:="comparison_group"]
setnames(sometimes_observed_background_mutations, "matching_background_mutation_id", "ID")

final_background_mutations = rbind(background_mutations, sometimes_observed_background_mutations, fill=TRUE) %>% unique

# --- DEBUGGING: Explicit Distance Check ---
dm_debug = copy(donorMatches)
dm_debug[, c("r_chr", "r_start", "r_end", "r_ref", "r_alt") := tstrsplit(real_mutation_id, "_")]
dm_debug[, c("b_chr", "b_start", "b_end", "b_ref", "b_alt") := tstrsplit(matching_background_mutation_id, "_")]
dm_debug[, distance := abs(as.numeric(r_start) - as.numeric(b_start))]

print("DEBUG - Coordinate Distance Check:")
print(paste("Max distance between matched real and background IDs:", max(dm_debug[["distance"]], na.rm=TRUE)))
print(paste("Min distance between matched real and background IDs:", min(dm_debug[["distance"]], na.rm=TRUE)))
# ----------------------------------------

# Restore the stump variants filter to maintain matched statistical integrity.
# We export the missing backgrounds to diagnose why AlphaGenome dropped them.
mutations_with_background = donorMatches[matching_background_mutation_id %in% qc_passing_background[,unique(ID)], real_mutation_id]
mutations_without_background = donorMatches[!matching_background_mutation_id %in% qc_passing_background[,unique(ID)], real_mutation_id]
stump_variants = setdiff(mutations_without_background, mutations_with_background)

missing_background_ids = donorMatches[!matching_background_mutation_id %in% qc_passing_background[,unique(ID)], unique(matching_background_mutation_id)]
fwrite(data.table(missing_ID = missing_background_ids), "DEBUG_missing_backgrounds.txt", sep="\t")

#Patient count per region
sufficient_patient_counts = real_mutations[,.(donor_id, associated_region,subgroup)] %>% 
  .[,.N,by=.(associated_region,subgroup)] %>% .[N>=minimum_donors,]

region_selection = initial_region_selection[sufficient_patient_counts, on=.(associated_region==associated_region,cancer_type==subgroup), nomatch=0]

# Combine, drop missing predictions, and enforce the stump variants filter
combined_mutations = rbind(real_mutations, final_background_mutations, fill=TRUE)

print("DEBUG - V1 Tracking:")
print(paste("Total real mutations before filters:", nrow(combined_mutations[label == "observed_mutation"])))

mutations_no_na = combined_mutations[!is.na(model_predictions)]
surviving_real = mutations_no_na[label == "observed_mutation"]
print(paste("Real mutations surviving NA drop:", nrow(surviving_real)))

# Diagnose exactly why surviving real mutations fail the stump filter
surviving_real_ids = surviving_real[, unique(ID)]
expected_backgrounds = donorMatches[real_mutation_id %in% surviving_real_ids, unique(matching_background_mutation_id)]
missing_expected = setdiff(expected_backgrounds, qc_passing_background[, unique(ID)])

print(paste("For the", length(surviving_real_ids), "surviving real mutations, we expect", length(expected_backgrounds), "backgrounds."))
print(paste("Of those expected backgrounds,", length(missing_expected), "are missing from qc_passing_background."))
if (length(missing_expected) > 0) {
    print("First 5 unexpectedly missing background IDs:")
    print(head(missing_expected, 5))
}

all_mutations_w_pancancer = mutations_no_na[!ID %in% stump_variants]
print(paste("Real mutations surviving stump filter:", nrow(all_mutations_w_pancancer[label == "observed_mutation"])))

# Keep only mutations in regions with sufficient patients and valid genes
enough_mutations = all_mutations_w_pancancer[region_selection, on=.(associated_region == associated_region, subgroup == cancer_type), nomatch=0]
enough_mutations = enough_mutations[!is.na(gene)]
# -----------------------------------------------------------------

#diagnose why there are not enough x or y observations for some regions
diagnostic = enough_mutations[, .(nrow(.SD[label == "observed_mutation"]), nrow(.SD[label == "comparison_group"])), by=.(subgroup,associated_region,gene)]

# --- DEBUGGING: Export diagnostic counts and filter invalid groups ---
fwrite(diagnostic, "DEBUG_diagnostic_counts.tsv", sep="\t")

valid_groups = diagnostic[V1 >= 1 & V2 >= 1, .(subgroup, associated_region, gene)]
if (nrow(valid_groups) == 0) {
    stop("DEBUG HALT: All groups have 0 observations in either the observed (V1) or comparison (V2) group. See DEBUG_diagnostic_counts.tsv.")
}
enough_mutations = enough_mutations[valid_groups, on=.(subgroup, associated_region, gene), nomatch=0]
# ---------------------------------------------------------------------

auc_denominator = diagnostic[, n0n1 := V1*V2][, c('subgroup', 'associated_region', 'gene', 'n0n1')]
n_mutations = diagnostic[, n0 := V1][, c('subgroup', 'associated_region', 'gene', 'n0')]

#jitter results to break ties
jittered_enrichment <- function(enough_mutations, number_of_tests, test, iter, auc_denominator, n_mutations){
  all_results = lapply(1:iter, function(x){
    enough_mutations[, noise := runif(nrow(enough_mutations), -0.00001, 0.00001)]
    enough_mutations[, jittered_predictions := model_predictions + noise]
    calculate_enrichments(enough_mutations, test, auc_denominator, n_mutations)
  }) %>% do.call("rbind", .) %>% .[, .(w_pvalue = median(w_pvalue), auc = median(auc), n_observed_mutations = median(n_observed_mutations)), by=.(subgroup, associated_region, gene)]
  
  significance_threshold = 0.05/number_of_tests
  
  all_results[, is_bf_significant := FALSE]
  all_results[w_pvalue < significance_threshold, is_bf_significant := TRUE]
  
  return(all_results)
}

calculate_enrichments <- function(enough_mutations, test, auc_denominator, n_mutations){
    
  if(test == "greater"){
    enough_mutations[, use_predictions:= abs(jittered_predictions)]
  } else{
    enough_mutations[, use_predictions:= jittered_predictions]
  }
  
  all_wilcox_against_null = enough_mutations[, wilcox.test(.SD[label == "observed_mutation",use_predictions], 
                                                     .SD[label == "comparison_group",use_predictions], 
                                                     alternative = test)[c(1,3)], 
                                                     by=.(subgroup,associated_region,gene)]

  setnames(all_wilcox_against_null, "p.value", "w_pvalue")
  setnames(all_wilcox_against_null, "statistic", "w_statistic")
  all_wilcox_against_null[auc_denominator, on=.NATURAL, auc := w_statistic / n0n1]
  all_wilcox_against_null[n_mutations, on=.NATURAL, n_observed_mutations := n0]
      
  return(all_wilcox_against_null[order(w_pvalue)])
}

number_of_tests = region_selection[,length(unique(associated_region))]

#the comparisons are pre filtered to the regions that had enough mutations
all_wilcox_against_null <- jittered_enrichment(enough_mutations, number_of_tests, "two.sided", 5, auc_denominator, n_mutations)
all_wilcox_against_null_abs <- jittered_enrichment(enough_mutations, number_of_tests, "greater", 5, auc_denominator, n_mutations)

fwrite(all_wilcox_against_null, glue("regions_summary.tsv"), sep = "\t")
fwrite(all_wilcox_against_null_abs, glue("regions_summary_abs.tsv"), sep = "\t")
fwrite(enough_mutations, glue("tested_mutations.tsv.gz"), sep = "\t")
