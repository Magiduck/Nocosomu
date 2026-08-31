process createISM {
  label 'rscript'
  
  input:
    tuple val(ID), path(bed, stageAs: 'regions.bed')
  
  output:
    tuple val(ID), path("ism.tsv")
  
  script:
    template 'ism.R'
  
}