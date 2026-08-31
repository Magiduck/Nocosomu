#!/usr/bin/env Rscript
library(data.table)

# Nextflow injects the variables directly into this template
input_file = "${vcflikeFile}"
output_file = "${vcflikeFile.baseName}_alphagenome_format.tsv"

dt = fread(input_file)

if ("ID" %in% colnames(dt)) {
    setnames(dt, "ID", "variant_id")
} else if ("id" %in% colnames(dt)) {
    setnames(dt, "id", "variant_id")
} else if (!"variant_id" %in% colnames(dt)) {
    chr_col = if("#CHROM" %in% colnames(dt)) "#CHROM" else if("CHROM" %in% colnames(dt)) "CHROM" else "seqnames"
    pos_col = if("POS" %in% colnames(dt)) "POS" else "start"
    dt[, variant_id := paste0(get(chr_col), ":", get(pos_col), "_", REF, "/", ALT)]
}

if ("#CHROM" %in% colnames(dt)) {
    setnames(dt, "#CHROM", "CHROM")
} else if ("seqnames" %in% colnames(dt)) {
    setnames(dt, "seqnames", "CHROM")
}

# FIX: Convert 0-based BED coordinates to 1-based VCF 'POS'
if (!"POS" %in% colnames(dt)) {
    if ("end" %in% colnames(dt)) {
        # In BED format, the 'end' coordinate for a SNP aligns with the 1-based POS
        dt[, POS := end]
    } else if ("start" %in% colnames(dt)) {
        dt[, POS := start + 1]
    } else if ("STA" %in% colnames(dt)) {
        dt[, POS := STA + 1]
    }
}

dt_out = dt[, .(variant_id, CHROM = as.character(CHROM), POS, REF, ALT)]
dt_out[!grepl("^chr", CHROM, ignore.case=TRUE), CHROM := paste0("chr", CHROM)]

fwrite(dt_out, output_file, sep="\t", quote=FALSE)
