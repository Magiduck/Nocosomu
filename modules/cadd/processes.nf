
process cadd_tsv_like {
  label 'rscript'
  
  input:
  tuple val(ID), path(tsv_file, stageAs: 'piece.csv')
  each path(cadd_scores)

  
  output:
  tuple val(ID), val('CADD'), path('cadd_prediction.tsv')
  
  script:
    """
    ${moduleDir}/resources/usr/bin/cadd_scores.R
    """
  
}

process cadd_merge_folds {
  label 'rscript'
  
  input:
  tuple val(ID), val(groupID), path('?_prediction.tsv')
  
  output:
  tuple val(ID), val(groupID), path("cadd_merged.tsv", arity: '1')
  
  script:
    template 'merge_cadd_folds.R'
  
}
