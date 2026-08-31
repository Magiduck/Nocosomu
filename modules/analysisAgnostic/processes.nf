process transform_to_vcflike {
  label 'rscript'
  
  input:
  tuple val(ID), path(mutations_file)
  
  output:
  tuple val(ID), path("vcf_like.tsv")
  
  script:
  template 'transform_to_vcflike.R'
  
}

process region_subselection {
  label 'rscript'
  
  input:
  tuple val(ID), path(elegible_regions, stageAs: "elegible_regions.csv")
  each path(region_subselection, stageAs: "region_subselection.txt") //Weird
  
  output:
  tuple val(ID), path("elegible_regions_to_run.csv")
  
  script:
  template 'region_subselection.R'
  
}

process mutation_subselection {
  label 'rscript'
  
  input:
  tuple val(ID), path(mutations, stageAs: "mutations.csv"), path(region_subselection, stageAs: "elegible_regions.csv")
  
  output:
  tuple val(ID), path("mutations_to_run.csv")
  
  script:
  template 'mutation_subselection.R'
  
}

process add_trinucleotide {
  label 'pyscript'
  
  input:
  tuple val(ID), path(mutations_file)
  path reference_genome
  
  output:
  tuple val(ID), path("all_snvs_with_trinucleotide.tsv")
  
  script:
  """
  add_trinucleotide_context.py --inputPath $mutations_file --genomePath $reference_genome -o all_snvs_with_trinucleotide.tsv
  """
  
}

process agreggate_analysis_files {
  label 'rscript'
  
  input:
  tuple val(chromosome), val(analysis_id), val(model_id), path(background_mutations, stageAs: "background_mutations?.tsv"), val(dummy), path(real_mutations, stageAs: "real_mutations?.tsv"), path(background, stageAs: "background?.tsv"), path(background, stageAs: "no_background?.tsv"), path(region_selection, stageAs: "region_selection?.tsv")
  each path(region_aggregation, stageAs: "region_aggregation.tsv")
  
  output:
  tuple val('aggregation'), val(analysis_id), val(model_id), path("background_mutations.tsv"), val(dummy), path("real_mutations.tsv"), path("background.tsv"), path("no_background.tsv"), path("region_selection.tsv")
  
  when:
  region_aggregation.name != 'NO_FILE'
  
  script:
  template 'region_aggregation.R'
}


process wilcoxon_compare {
  label 'rscript'
  
  input:
  tuple val(chromosome), val(analysis_id), val(model_id), path(background_mutations, stageAs: "background_mutations.tsv"), val(dummy), path(real_mutations, stageAs: "real_mutations.tsv"), path(background, stageAs: "background.tsv"), path(background, stageAs: "no_background.tsv"), path(region_selection, stageAs: "region_selection.tsv")
  each minimum_donors
  
  output:
  tuple val(chromosome), val(analysis_id), val(model_id), path("regions_summary.tsv"), emit: summary
  tuple val(chromosome), val(analysis_id), val(model_id), path("regions_summary_abs.tsv"), emit: summary_abs
  tuple val(chromosome), val(analysis_id), val(model_id), path("tested_mutations.tsv.gz"), emit: tested_mutations
  
  script:
  template 'wilcox_compare.R'
  
}

process prediction_post_process_borzoi {
    label 'rscript'
    
    input:
    tuple val(ID), val(groupID), path(prediction_file, stageAs: "raw_predictions.tsv")

    output:
    tuple val(ID), val(groupID), path('processed_predictions.tsv')
    
    script:
    template 'prediction_post_process_borzoi.R'
    

}

process prediction_post_process_alphagenome {
    label 'rscript'
    
    input:
    tuple val(ID), val(groupID), path(prediction_file, stageAs: "raw_predictions.tsv")
    each gene_annotation

    output:
    tuple val(ID), val(groupID), path('processed_predictions.tsv')
    
    script:
    template 'prediction_post_process_alphagenome.R'
    

}

process prediction_post_process_sure {
    label 'rscript'
    
    input:
    tuple val(ID), val(groupID), path(prediction_file, stageAs: "raw_predictions.tsv")

    output:
    tuple val(ID), val(groupID), path('processed_predictions.tsv')
    
    script:
    template 'prediction_post_process_SuRE.R'
    

}

process prediction_post_process_cadd {
    label 'rscript'
    
    input:
    tuple val(ID), val(groupID), path(prediction_file, stageAs: "raw_predictions.tsv")

    output:
    tuple val(ID), val(groupID), path('processed_predictions.tsv')
    
    script:
    template 'prediction_post_process_cadd.R'
    

}
