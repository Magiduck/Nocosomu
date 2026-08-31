#!/usr/bin/env Rscript
library(data.table)

dt = fread("raw_predictions.tsv")

cat("DEBUG: Total rows loaded from merged prediction file:", nrow(dt), "\n")

if ("variant_name" %in% colnames(dt)) {
    setnames(dt, "variant_name", "ID")
} else if ("variant_id" %in% colnames(dt)) {
    setnames(dt, "variant_id", "ID")
}

if ("raw_score" %in% colnames(dt) && !"model_predictions" %in% colnames(dt)) {
    setnames(dt, "raw_score", "model_predictions")
}

if ("gene_name" %in% colnames(dt)) {
    setnames(dt, "gene_name", "gene")
} else if ("gene_id" %in% colnames(dt)) {
    setnames(dt, "gene_id", "gene")
}

if (nrow(dt) > 0 && grepl("^chr[0-9a-zA-Z]+:[0-9]+:[a-zA-Z]+>[a-zA-Z]+", dt[["ID"]][1])) {
    dt[, c("ag_chr", "ag_pos", "ag_refalt") := tstrsplit(ID, ":")]
    dt[, ag_chr := sub("^chr", "", ag_chr)]
    dt[, ag_pos := as.integer(ag_pos)]
    dt[, c("ag_ref", "ag_alt") := tstrsplit(ag_refalt, ">")]
    
    dt[, ag_start := ag_pos - 1] 
    dt[, ag_end := ag_pos]
    
    dt[, ID := paste(ag_chr, ag_start, ag_end, ag_ref, ag_alt, sep = "_")]
    dt[, c("ag_chr", "ag_pos", "ag_end", "ag_refalt", "ag_ref", "ag_alt", "ag_start") := NULL]
}

if (!"ID" %in% colnames(dt) || !"model_predictions" %in% colnames(dt)) {
    stop(paste("Missing required columns. Available columns are:", paste(colnames(dt), collapse=", ")))
}

cat("DEBUG: Total unique variant IDs successfully parsed:", length(unique(dt[["ID"]])), "\n")

# fwrite(data.table(ID = dt[["ID"]]), "DEBUG_loaded_ids.txt", sep="\t", col.names=FALSE)

group_cols <- c("ID")
if ("tri_context_ref" %in% colnames(dt)) group_cols <- c(group_cols, "tri_context_ref")
if ("tri_context_alt" %in% colnames(dt)) group_cols <- c(group_cols, "tri_context_alt")

if ("gene" %in% colnames(dt)) {
    bed_path <- "$gene_annotation"
    bed_dt <- fread(cmd = paste("zcat -f", bed_path))
    setnames(bed_dt, c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "bed_gene_versioned"))
    
    bed_dt[, chrom := sub("^chr", "", chrom)]
    bed_dt[, bed_gene := sub("[.][0-9]+", "", bed_gene_versioned)]
    setkey(bed_dt, chrom, chromStart, chromEnd)
    
    dt[, c("v_chr", "v_start", "v_end") := tstrsplit(ID, "_", keep = 1:3)]
    dt[, v_start := as.integer(v_start)]
    dt[, v_end := as.integer(v_end)]
    setkey(dt, v_chr, v_start, v_end)
    
    # Execute physical overlap but retain all unmatched variants
    dt_mapped <- foverlaps(dt, bed_dt, by.x = c("v_chr", "v_start", "v_end"), by.y = c("chrom", "chromStart", "chromEnd"), nomatch = NA)
    
    # Calculate absolute score for secondary priority sorting
    dt_mapped[, abs_score := abs(model_predictions)]
    
    # Flag priorities
    dt_mapped[, physical_overlap := !is.na(bed_gene)]
    dt_mapped[, exact_match := physical_overlap & (gene == bed_gene)]
    
    # Sort dataset: Exact matches first, then physical overlaps (highest score), then no overlaps (highest score)
    setorder(dt_mapped, ID, -exact_match, -physical_overlap, -abs_score)
    
    # Retain the single highest priority row per unique variant-gene interaction
    dt_out <- dt_mapped[, .SD[1], by = .(ID, bed_gene)]
    
    # Force the gene identifier to match the BED file if a physical overlap exists, 
    # bypassing the GENCODE version mismatch.
    dt_out[physical_overlap == TRUE & exact_match == FALSE, gene := bed_gene]
    
    # Diagnostic counts
    exact_count <- nrow(dt_out[exact_match == TRUE])
    overlap_mismatch_count <- nrow(dt_out[physical_overlap == TRUE & exact_match == FALSE])
    no_overlap_count <- nrow(dt_out[physical_overlap == FALSE])
    total_retained <- exact_count + overlap_mismatch_count + no_overlap_count
    expected_variants <- length(unique(dt[["ID"]]))
    
    cat("\n=== MAPPING DIAGNOSTICS ===\n")
    cat("1. Exact BED overlap and gene match:             ", exact_count, "\n")
    cat("2. Physical overlap only (gene ID harmonized):   ", overlap_mismatch_count, "\n")
    cat("3. No physical overlap (max score fallback):     ", no_overlap_count, "\n")
    cat("------------------------------------------------------\n")
    cat("TOTAL variant-gene interactions retained:        ", total_retained, "\n")
    cat("UNIQUE variants based on input:                  ", expected_variants, "\n")
    cat("(Note: Retained count may exceed unique variants due to multi-gene overlaps)\n")
    cat("===========================\n\n")
    
    cols_to_keep <- c(group_cols, "gene", "model_predictions")
    dt_out <- dt_out[, ..cols_to_keep]
} else {
    dt_out = dt[, .(model_predictions = mean(model_predictions, na.rm = TRUE)), by = group_cols]
    dt_out[, gene := "unknown"]
}

fwrite(dt_out, "processed_predictions.tsv", sep=",", quote=FALSE)
