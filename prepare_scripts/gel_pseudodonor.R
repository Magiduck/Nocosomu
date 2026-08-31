library(data.table)
library(splitstackshape)
library(magrittr)

fold_result = list.files(path = "./", 
                         pattern = "*frequency_table.tsv",
                         full.names = T)

all_gel = lapply(fold_result, function(x){
  m <- fread(x)
  m[, gel_cancer := basename(x)]
}) %>% do.call("rbind", .)

all_hmf_exp <- expandRows(all_gel, "count", count.is.col = TRUE, drop = TRUE) %>% as.data.table

#all_hmf_exp[, pre_donor_id := paste0("gel_pseudo_donor_", seq_len(.N)), by=.(V1,V2,V3,V4,V5)]
all_hmf_exp[, pre_donor_id := paste0("gel_pseudo_donor_", seq_len(.N)), by=.(`#CHROM`, POS, REF, gel_cancer)]

all_hmf_exp[, donor_id := paste0(pre_donor_id, "_",.GRP), by=gel_cancer]
all_hmf_exp[, pre_donor_id := NULL]

fwrite(all_hmf_exp, "./gel_all_pseudodonorQCed.tsv.gz", 
       sep = '\t')

