
process borzoi_tsv_like {
  label 'gpu'
  
  input:
  tuple val(ID), path(tsv_file)
  each path(borzoi_model)
  each path(reference_genome)
  each path(targets_borzoi)
  each path(gene_annotation)
  each path(borzoi_params)
  each borzoi_flags
  
  output:
  tuple val(ID), path('sed.h5')
  
  script:
  """
  ${params.borzoi_dir}/borzoi_sed.py \
    --stats SED,logSED,D1,D2,nD2,nDi,JS,REF,ALT \
    -f $reference_genome \
    -g $gene_annotation \
    -t $targets_borzoi \
    -o . \
    $borzoi_params \
    $borzoi_model $borzoi_flags\
    <(awk '{print \$1,\$3,\$6,\$4,\$5}' OFS='\t' $tsv_file)
  """
  
}

process borzoi_extract_sed {
  
  input:
  tuple val(ID), path(sed_file)
  each path(targets_borzoi)
  
  output:
    tuple val(ID), path('*_predictions.tsv', arity: '32')
  
  script:
    """
    ${params.borzoi_dir}/extract_sed_from_h5_gnomAD2.py \
        -i ./sed.h5 \
        -t $targets_borzoi \
        -o _predictions.tsv
    """
  
}

process borzoi_merge_folds {
  label 'rscript'
  
  input:
  tuple val(ID), val(groupID), path('?_predictions.tsv')
  val borzoi_folds
  
  output:
  tuple val(ID), val(groupID), path("borzoi_merged_*.tsv", arity: '1')
  
  script:
    template 'merge_borzoi_folds.R'
  
}

process borzoi_splice_effect {
  label 'rscript'
  
  input:
  tuple val(ID), val(groupID), path(prediction_file)
  each path(conversion_table)
  
  output:
  tuple val(ID), val(groupID), path("spliceEffect_*.tsv", arity: '1')
  
  script:
    template 'borzoi_splice_effect.R'
  
}

process borzoi_format {
  label 'rscript'
  
  input:
  tuple val(ID), path(big_file)
  val borzoi_columns
  
  output:
  tuple val(ID), path('output_file.tsv')
  
  script:
  template 'format_borzoi.R'
  
}
