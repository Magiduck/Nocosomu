process sureResNet_tsv_like {
  label 'gpu'
  
  input:
  tuple val(ID), path(tsv_file)
  each path(sure_model)
  each path(reference_genome)

  output:
  tuple val(ID), val(sure_model.simpleName), path("prediction.tsv")
  
  script:
    """
        ${params.sure_dir}/predict_sure_from_tsv.py \
        --model $sure_model \
        --input $tsv_file \
        --genome $reference_genome \
        --output prediction.tsv \
        --batchsize 1000 \
        --L_max 600
    """
  
}

process parm_merge {
  label 'rscript'
  
  input:
  tuple val(ID), val(groupID), path('?_predictions.tsv')
  
  output:
  tuple val(ID), val(groupID), path("*Strand.tsv", arity: '3')
  
  script:
    template 'merge_sure_data.R'
  
}
