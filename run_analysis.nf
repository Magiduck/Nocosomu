
workflow {
  elegible_regions = ICGC_elegible_regions(Channel.fromPath("$params.dataset"), Channel.fromPath("$params.bed"))
  
  #readn null model
  chromosomes = Channel.fromList("$params.chromosomes".tokenize( ',' ))
  mutation_files = chromosomes
  | map {tuple(it,path("${params.bedID}/${params.datasetID}/$it/merged_mutations_in_region.tsv"))}
  | sbs_mutations
     
     
  #run bakground model 
  all_backgrounds = alternativeAlleles_background(elegible_regions.join(mutation_files, by:0))
  
  
  
  #make predictions for background and null
  
  #call Wicoxon
  
  #merge results by steps
  
  #publish final file

}