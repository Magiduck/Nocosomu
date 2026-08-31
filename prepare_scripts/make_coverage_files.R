library(data.table)
library(magrittr)
library(glue)

#https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
what = '19'
m <- fread(glue("downloads/hg{what}/hg{what}.chrom.sizes"))
m[,V3:=0]
fwrite(m[,.(V1,V3,V2)], glue("downloads/hg{what}/sizes_hg{what}.bed"), sep = "\t", col.names = F)

#module load BEDTools
#zcat k36.umap.bedgraph.gz | grep -E "1\\.|0\\.9" > k36.umap.bedgraph_whitelist.bed
#bedtools subtract -a sizes_hg19.bed -b k36.umap.bedgraph_whitelist.bed > k36.umap.bedgraph_blacklist_updated.bed

m <- fread(glue("downloads/hg{what}/k36.umap.bedgraph_blacklist_updated.bed", skip=1))
all_blacklist = m[ , list(chromosome = V1, POS = seq(V2, V3)), by = 1:nrow(m)]

all_blacklist[,fwrite(.SD[,.(POS)],glue("artifacts/hg{what}/{.BY}_k36.umap.multiread_blacklist.txt.gz"), col.names = F), by=chromosome]
