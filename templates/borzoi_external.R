
Channel.fromPath("$params.vcflike")
.splitText(by: 10, file: true, keepHeader: true)
.map{ it -> tuple("$params.vcflike", it)}
.into {borzoi_external_in}

Channel.fromList("$params.borzoi_folds".tokenize( ',' )) //'f0,f1,f2,f3'
.into {borzoi_folds; borzoi_folds_external}

process borzoi_tsv_like {
  
  input:
  tuple path(og_tsv), tsv_file, borzoi_fold from borzoi_external_in.combine(borzoi_folds_external)
  file reference_genome from file("$params.reference_genome")
  file id_gene from file("$params.vcflikeID")
  file targets_borzoi from file("${params.borzoi_dir}/targets_gtex.txt")
  file gene_annotation from file("${params.borzoi_dir}/gencode41_basic_nort.gtf")
  file borzoi_params from file("${params.borzoi_dir}/params_pred.json")
  
  output:
    tuple path('*_predictions.tsv') into borzoi_out_external
  
  when:
    params.stage == "run_predictions_external"
  params.models =~ /borzoi/ && params.from_cache == 'NOT_PROVIDED'
  
  script:
    def track
  if (params.models =~ /gtexBlood/){
    track = 'RNA:blood'
  }
  """
        ${params.borzoi_dir}/borzoi_sed.py \
        --stats SED,logSED,D1,D2,nD2,nDi,JS,REF,ALT \
        -f $reference_genome \
        -g $gene_annotation \
        -t $track \
        -t $targets_borzoi \
        -o . \
        $borzoi_params \
        ${params.borzoi_dir}/saved_models/${borzoi_fold}/model0_best.h5 \
        $tsv_file

        ${params.borzoi_dir}/extract_sed_from_h5.py \
        -i ./sed.h5 \
        -a $id_gene \
        -f $tsv_file \
        -t $targets_borzoi \
        -x $track \
        -o ${og_tsv.baseName}_${borzoi_fold}_predictions.tsv

        
    """
  
}

borzoi_out_external
.collectFile(keepHeader:true)
.set{borzoi_collected}

process borzoi_post_processing {
  
  publishDir "${params.borzoi_outdir}", mode: 'copy', saveAs: { filename -> "$filename" }
  
  input:
    tuple tsv_file from borzoi_collected
  
  output:
    path '*.tsv'
  
  script:
    """
    cp $tsv_file ${tsv_file.baseName}.tsv
    """
  
}