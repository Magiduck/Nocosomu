library(data.table)
library(glue)
library(magrittr)
library(tidyr)

what <- "v45"
what <- "v45lift37"

all_gencode <- fread(glue("/groups/umcg-fg/tmp02/projects/non-coding-somatic/shared/gencode/downloads/gencode.{what}.annotation.gtf.gz")) %>% .[V3=='exon' & grepl("protein_coding", V9) & grepl("MANE", V9)]

all_gencode_separated <- separate(all_gencode, V9, paste0("A", 1:8), sep = ';') %>% .[,c(1:8,16)] %>% unique

all_gencode_separated[, dummy_id := gsub("exon_id", "gene_id", A8)]
all_gencode_separated[, dummy_id := gsub("exon_id", "gene_id", dummy_id)]

all_gencode_separated[, ks := paste0(dummy_id,"; ",A8,";")]


fwrite(all_gencode_separated[,c(1:8,11)], glue("/groups/umcg-fg/tmp02/projects/non-coding-somatic/shared/gencode/artifacts/gencode.{what}.exonsAsGenes.gtf.gz"), col.names = F, sep = "\t", quote = F)

#introns
all_gencode <- fread(glue("/groups/umcg-fg/tmp02/projects/non-coding-somatic/shared/gencode/downloads/gencode.{what}.annotation.gtf.gz")) %>% .[V3=='exon' & grepl("protein_coding", V9) & grepl("MANE", V9)]
all_gencode_separated <- separate(all_gencode, V9, paste0("A", 1:8), sep = ';') %>% .[,c(1:8,10,15,16)] %>% unique

all_gencode_separated[,exon_number := as.integer(gsub("exon_number ", "", A7))]
all_gencode_separated[,next_exon_number := exon_number + 1]

all_gencode_introns = all_gencode_separated[all_gencode_separated, on=.(A2==A2, next_exon_number == exon_number), nomatch=0]
all_gencode_introns[V7=="+",V4:=V5]
all_gencode_introns[V7=="+",V5:=i.V4]
all_gencode_introns[V7=="-",V5:=V4]
all_gencode_introns[V7=="-",V4:=i.V5]


all_gencode_introns[, dummy_id := gsub("exon_id", "gene_id", A8)]

all_gencode_introns[, ks := paste0(dummy_id,"; ",A8,";")]
all_gencode_introns[, ks := gsub("ENSE","ENSI", ks)]


fwrite(all_gencode_introns[,c(1:8,26)], glue("/groups/umcg-fg/tmp02/projects/non-coding-somatic/shared/gencode/artifacts/gencode.{what}.intronsAsGenes.gtf.gz"), 
       col.names = F, sep = "\t", quote = F)


#make mappingfile
all_gencode_introns[,transcript_id := gsub("transcript_id ", "", A2)]
all_gencode_introns[,exon_5end := gsub("exon_id ", "", A8)]
all_gencode_introns[,exon_3end := gsub("exon_id ", "", i.A8)]
all_gencode_introns[,intron_id := paste0("intron_", exon_number, "_", next_exon_number)]

final_table = all_gencode_introns[,c(1,4,5,7,27:30)]
setnames(final_table, c("chromosome",
                        "start",
                        "end",
                        "strand",colnames(final_table)[5:8]))

fwrite(final_table, glue("/groups/umcg-fg/tmp02/projects/non-coding-somatic/shared/gencode/artifacts/gencode.{what}.intronsToExonsAndTranscript.tsv.gz"), sep = "\t", quote = F)

both_featues = rbind(fread(glue("/groups/umcg-fg/tmp02/projects/non-coding-somatic/shared/gencode/artifacts/gencode.{what}.exonsAsGenes.gtf.gz")),
                     fread(glue("/groups/umcg-fg/tmp02/projects/non-coding-somatic/shared/gencode/artifacts/gencode.{what}.intronsAsGenes.gtf.gz")))

fwrite(both_featues, glue("/groups/umcg-fg/tmp02/projects/non-coding-somatic/shared/gencode/artifacts/gencode.{what}.exonsAndIntronsAsGenes.gtf.gz"), col.names = F, sep = "\t", quote = F)
