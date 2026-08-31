process alternativeAlleles_background {
  label 'rscript'
  
  input:
  tuple val(ID), path(mutations_file), path(ism_file)
  
  output:
  tuple val(ID), path("alternativeAlleles_background.tsv"), emit: background
  tuple val(ID), path("alternativeAlleles_donorMapping.tsv"), emit: donorMapping
  script:
  template 'alternative_allele_background.R'
  
}

process samePromoter_background {
  label 'rscript'
  
  input:
  tuple val(ID), path(mutations_file), path(ism_file)
  
  output:
  tuple val(ID), path("samePromoter_background.tsv"), emit: background
  tuple val(ID), path("samePromoter_donorMapping.tsv"), emit: donorMapping
  
  script:
  template 'samePromoter_background.R'
  
}

process otherTumors_background {
  label 'rscript'
  
  input:
  tuple val(ID), path(mutations_file), path(ism_file)
  
  output:
  tuple val(ID), path("otherTumors_background.tsv"), emit: background
  tuple val(ID), path("otherTumors_donorMapping.tsv"), emit: donorMapping
  
  script:
  template 'otherTumors_background.R'
  
}

process no_background {
  label 'rscript'
  
  input:
  tuple val(ID), path(mutations_file), path(ism_file)
  
  output:
  tuple val(ID), path("sbs_no_background.tsv")
  
  script:
  template 'no_background.R'
  
}

process region_metadata {
  label 'rscript'
  
  input:
  tuple val(ID), path(mutations_file), path(elegible_regions)
  
  output:
  tuple val(ID), path("snvs_in_region.tsv")
  
  script:
  template 'region_metadata.R'
  
}
