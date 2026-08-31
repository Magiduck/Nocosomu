library(data.table)
library(magrittr)
library(glue)
library(stringr)

what = '19'


#low_complexity

m <- fread(cmd=glue("zcat /groups/umcg-fg/tmp01/projects/non-coding-somatic/shared/anySTR/downloads/hg{what}.fa.out.gz | grep complexity"))
all_low_complexity = m[ , list(chromosome = V5, POS = seq(V6, V7)), by = 1:nrow(m)]
all_low_complexity[,fwrite(.SD[,.(POS)],glue("artifacts/hg{what}/{.BY}_repeatMaskerLowComplexity.txt.gz"), col.names = F), by=chromosome]
