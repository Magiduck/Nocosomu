include { import_ICGC ; import_GE ; elegible ; remove_germline ; filter_on_blacklist } from './processes.nf'
include { filter_basepairs as filter_repeats; filter_basepairs as filter_lowComplexity; filter_basepairs as filter_on_coverage; filter_basepairs as filter_on_mappability; filter_basepairs as filter_coding} from './processes.nf'
include { CollectFileTuple } from './../general'
include { publish_file2 as publish_import } from './../general'
include { publish_file2 as publish_elegible } from './../general'

workflow {
  chromosomes = Channel.fromList("$params.chromosomes".tokenize( ',' ))
  germline_variants_ch = chromosomes.map {tuple(it, "${params.germline_variants}/${it}_TrueSNV_TruePASS.csv")}
  filtered_bed = filter_on_blacklist("$params.blacklist", "$params.bed")

  if ("$params.data_source" =~ /ICGC/) {
    import_ch = Import_ICGC_parquet(Channel.fromPath("$params.dataset"), filtered_bed)
  } else if ("$params.data_source" =~ /GE/) {
    import_ch = Import_GE_tsv(Channel.fromPath("$params.dataset"), filtered_bed)
  }
  
  qc_filtered_variants = import_ch
  | map { tuple(it[0], it[1], it[2], it[3], "${params.below_coverage}/chr${it[2]}_positions_below_30.csv.gz", "$params.skipQC" ==~ /COV30/) }
  | filter_on_coverage
  | map { tuple(it[0], it[1], it[2], it[3], "${params.within_repeat}/chr${it[2]}_harmonizedDietlein.txt.gz", "$params.skipQC" ==~ /REPEAT/) }
  | filter_repeats
  | map { tuple(it[0], it[1], it[2], it[3], "${params.within_repeat}/chr${it[2]}_repeatMaskerLowComplexity.txt.gz", "$params.skipQC" ==~ /COMPLEXITY/ ) }
  | filter_lowComplexity
  | map { tuple(it[0], it[1], it[2], it[3], "${params.low_mappable}/chr${it[2]}_hg38.k36.umap.multiread_blacklist.txt.gz", "$params.skipQC" ==~ /UMAPK36/ ) }
  | filter_on_mappability
  | map { tuple(it[0], it[1], it[2], it[3], "${params.coding}/chr${it[2]}_CDS.txt.gz", "$params.skipQC" ==~ /CDS/ ) }
  | filter_coding
  | map { tuple(it[0], it[1], it[2], it[3], "${params.germline_variants}/chr${it[2]}_TrueSNV_TruePASS.csv.gz", params.germlineAlleleFrequency) }
  | map { tuple(it[0], it[1], it[2], it[3], "${params.germline_variants}/gnomad.joint.v4.1.sites.chr${it[2]}_FAF.tsv.gz", params.germlineAlleleFrequency) }
  | remove_germline
  
  qc_filtered_variants
  | map {tuple("${it[0]}/${it[1]}/${it[2]}", it[3])}
  | publish_import

  input_ch = qc_filtered_variants.map { tuple(it[2], it[3])}.filter { it[1].countLines() > 1 } //If no elegible mutations are found then drop the file

  elegible_regions(Channel.fromPath("$params.dataset"), filtered_bed, input_ch)
  | map { tuple("${it[0]}/${it[1]}/${it[2]}", it[3]) }
  | publish_elegible

}

workflow Import_ICGC_parquet {
  take:
    icgc_parquet
    bed_file
  
  main:
  chromosomes = Channel.fromList("$params.chromosomes".tokenize( ',' ))
  
  //subgroups = Channel.fromPath("$params.subgroup_file")
  //  .splitCsv(header: true, sep: '\t')
  //  .map { it.my_cancer_type }
  //  .unique()
  //  .view()
  
  import_ICGC(icgc_parquet.combine(bed_file)
                          .combine(chromosomes),
              "$params.subgroup_file")
  
  emit:
  import_ICGC.out
  
  
}

workflow Import_GE_tsv {
  take:
    ge_tsv
    bed_file
  
  main:
  chromosomes = Channel.fromList("$params.chromosomes".tokenize( ',' ))
  
  import_GE(ge_tsv.combine(bed_file)
                          .combine(chromosomes),
              "$params.subgroup_file")
  
  emit:
  import_GE.out
  
  
}

workflow elegible_regions {
  take:
    icgc_parquet
    bed_file
    mutations_file_ch
  
  main:
  chromosomes = Channel.fromList("$params.chromosomes".tokenize( ',' ))
  
  input_ch = chromosomes.combine(icgc_parquet).combine(bed_file).join(mutations_file_ch)
  .map {tuple(it[1], it[2], it[0], it[3])}

  elegible(input_ch,
           "$params.subgroup_file",
           "$params.minimum_donors")

  emit:
  elegible.out
  
}
