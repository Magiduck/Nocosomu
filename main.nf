Channel.fromList("$params.chromosomes".tokenize( ',' ))
.into {chr_ism; chr_selection; chr_alternative_alleles; chr_regions; chr_mutatedBasePairs}

Channel.fromPath("$params.bed")
.into{bed_ism; bed_init; bed_compress; bed_ID; bed_mutatedBasePairs}

bed_ID
.combine(Channel.fromPath("$params.dataset"))
.into{bed_parquet; bed_dataset_for_analyses; bedDataset_models; bedDataset_Analyses; bedDataset_tissue}

process add_region_to_parquet {
    input:
    tuple bed_file, dataset, chromosome from bed_parquet.combine(chr_regions)
    
    output:
    path '*.parquet' into parquet_with_region_out
    
    when:
    params.stage == "add_region"
    
	script:
	template 'add_region_to_parquet.R'
}

process compress_parquet {
    publishDir "${params.project_folder}/${bed.simpleName}/", mode: 'copy', saveAs: { filename -> "${bed.simpleName}.parquet"}
    
    input:
    path parquet_files from parquet_with_region_out.flatten().collect(sort: true)
    path bed from bed_compress
    file reference_genome from file("$params.reference_genome")
    
    output:
    path 'unified.parquet'
    
	script:
        template 'unify_parquet.py'
}

process icgc_cancer_types {
    input:
    path bed_file from bed_init
    
    output:
    path 'subtypes_genes.csv' into subtypes_gene_selection
    path 'subtypes_chromosome.csv' into subtypes_chromosome_selection
    path 'subtypes.txt' into cancer_subtypes
    path 'genes.txt' into gene_selection_ch
    path 'genes.txt' into gene_selection_plot
    
    when:
    params.stage =~ /run/ && (params.mode =~ /ISM/ || params.mode =~ /mutatedBasePairs/) //hotfix
    
	script:
	template 'ICGC_cancer_types.R'
}

gene_selection_ch
.into {gene_selection_ch_ism; gene_selection_ch_mutatedBasePairs}

subtypes_gene_selection
.into {subtypes_gene_selection_alternativeAlleles; subtypes_gene_selection_samePromoter; subtypes_gene_selection_otherPromoters; subtypes_gene_selection_nearbyInPromoter}

cancer_subtypes
.splitCsv()
.into {cancer_types; cancer_types_wgene}

subtypes_chromosome_selection
.splitCsv(header: true)
.map {row -> tuple(row.chromosome, row.cancer_type)}
.filter { it[0] in "$params.chromosomes".tokenize( ',' ) }
.into {subtypes_chromosome_selection_alternative_alleles; subtypes_chromosome_selection_samePromoter; subtypes_chromosome_selection_otherPromoters; subtypes_chromosome_selection_nearbyInPromoter}

Channel.fromPath("$params.models")
.into{ models_single; models_id}

models_id
.combine(bed_dataset_for_analyses)
.into {models_predict_in; models_alternative_alleles_in; models_nearbyInPromoter_in; models_samePromoter_in; models_otherPromoters_in; models_plotISM_in}

models_single
.into {models_wgene; models_wAnalysis; models_wModel}

Channel.fromList("$params.analyses".tokenize( ',' ))
.into {analyses_wgene; analyses_wAnalysis}


process ism_inputs {
	input:
	tuple bed_file, chromosome from bed_ism.combine(chr_ism)
	each gene_selection from gene_selection_ch_ism

	output:
	path "*_ism.tsv" into ism_out
	

	when:
	params.mode == "ISM" && params.stage == "run_predictions"

	script:
	template 'ism.R'

	
}

process mutatedBasePairs_inputs {
	input:
	tuple bed_file, chromosome from bed_mutatedBasePairs.combine(chr_mutatedBasePairs)
	each gene_selection from gene_selection_ch_mutatedBasePairs

	output:
	path "*_ism.tsv" into mutatedBasePairs_out

	when:
	params.mode == "mutatedBasePairs" && params.stage == "run_predictions"

	script:
	template 'mutatedBasePair.R'
	
}

ism_out
.mix(mutatedBasePairs_out)
.flatten()
.filter{ !(it.simpleName =~ /NO_GENES/) }
.combine(models_predict_in)
.into {sure_in; sureResnet_in; enformer_in; vep_in; cache_in; borzoi_in}

process sure_tsv_like {

    publishDir "${params.project_folder}/${bed.simpleName}/${dataset.simpleName}/${model_pickle.simpleName}/artifacts/${params.mode}/", mode: 'copy', saveAs: { filename -> "$filename" }

    input:
    tuple tsv_file, model_pickle, bed, dataset from sure_in
    file reference_genome from file("$params.reference_genome")

    output:
    tuple path('*_predictions.tsv.gz'), model_pickle, bed, dataset into sure_out
    
    when:
    model_pickle =~ /SuRE4n.*pickle/ && params.from_cache == 'NOT_PROVIDED'

    script:
    """
        ${params.sure_dir}/predict_sure_from_tsv.py \
        --model $model_pickle \
        --input $tsv_file \
        --genome $reference_genome \
        --output ${tsv_file.baseName}_predictions.tsv \
        --batchsize 1000 \
        --L_max 600
        gzip ${tsv_file.baseName}_predictions.tsv
    """

}

process sureResnet_tsv_like {

    publishDir "${params.project_folder}/${bed.simpleName}/${dataset.simpleName}/${model_pickle.simpleName}/artifacts/${params.mode}/", mode: 'copy', saveAs: { filename -> "$filename" }

    input:
    tuple tsv_file, model_pickle, bed, dataset from sureResnet_in
    file reference_genome from file("$params.reference_genome")

    output:
    tuple path('*_predictions.tsv.gz'), model_pickle, bed, dataset into sureResnet_out
    
    when:
    model_pickle =~ /ResNet_Attention/ && params.from_cache == 'NOT_PROVIDED'

    script:
    """
        ${params.sure_dir}/predict_sure_from_tsv.py \
        --model $model_pickle \
        --input $tsv_file \
        --genome $reference_genome \
        --output ${tsv_file.baseName}_predictions.tsv \
        --batchsize 1000 \
        --L_max 600
        gzip ${tsv_file.baseName}_predictions.tsv
    """

}

process enformer_tsv_like {

    publishDir "${params.project_folder}/${bed.simpleName}/${dataset.simpleName}/${model_pickle.simpleName}/artifacts/${params.mode}/", mode: 'copy', saveAs: { filename -> "$filename" }

    input:
    tuple tsv_file, model_pickle, bed, dataset from enformer_in
    file reference_genome from file("$params.reference_genome")

    output:
    tuple path('*_predictions.tsv.gz'), model_pickle, bed, dataset into enformer_out
    
    when:
    
    model_pickle =~ /EleutherAI_enformer_official_rough/ && params.from_cache == 'NOT_PROVIDED'

    script:
    def cell_line
    def track
    if (model_pickle =~ /K562/){
      cell_line = 'K562'
    } else if (model_pickle =~ /HEPG2/){
      cell_line = 'HEPG2'
    }
    if (model_pickle =~ /CAGEtrack/){
      cell_line = 'CAGE'
    } else if (model_pickle =~ /DNASEtrack/){
      cell_line = 'DNASE'
    }
    """
        ${params.sure_dir_inhouse}/predict_enformer_from_tsv.py \
        -i $tsv_file \
        -f $reference_genome \
        -c $cell_line \
        -t $track \
        -o ${tsv_file.baseName}_predictions.tsv
        gzip ${tsv_file.baseName}_predictions.tsv
    """

}


Channel.fromList("$params.borzoi_folds".tokenize( ',' )) //'f0,f1,f2,f3'
.into {borzoi_folds; borzoi_folds_external}

process borzoi_tsv_like {

    input:
    tuple tsv_file, model_pickle, bed, dataset, borzoi_fold from borzoi_in.combine(borzoi_folds)
    file reference_genome from file("$params.reference_genome")
    file targets_borzoi from file("${params.borzoi_dir}/targets_gtex.txt")
    file gene_annotation from file("${params.gencode_annotation}")
    file borzoi_params from file("${params.borzoi_dir}/params_pred.json")

    output:
    tuple path('*_predictions.tsv'), model_pickle, bed, dataset into borzoi_out
    
    when:
    params.stage == "run_predictions"
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
        --ignore_allele_checks \
        ${params.borzoi_dir}/saved_models/${borzoi_fold}/model0_best.h5 \
        <(awk '{print \$1,\$3,\$6,\$4,\$5}' OFS='\t' $tsv_file)

        ${params.borzoi_dir}/extract_sed_from_h5_CGUT.py \
        -i ./sed.h5 \
        -t $targets_borzoi \
        -x $track \
        -o ${tsv_file.baseName}_${borzoi_fold}_predictions.tsv
        
    """

}

process from_cache {

    publishDir "${params.project_folder}/${bed.simpleName}/${dataset.simpleName}/${model_pickle.simpleName}/artifacts/${params.mode}/", mode: 'copy', saveAs: { filename -> "$filename" }

    input:
    tuple tsv_file, model_pickle, bed, dataset from cache_in

    output:
    tuple path('*_predictions.tsv.gz'), model_pickle, bed, dataset into cache_out
    
    when:
    params.from_cache != 'NOT_PROVIDED'

    script:
    template 'from_cache.R'

}

sure_out
.mix(sureResnet_out)
.mix(enformer_out)
.mix(cache_out)
.set {all_model_predictions_sure}

process prediction_post_process_sure {

    input:
    tuple prediction_file, model_pickle, bed, dataset from all_model_predictions_sure

    output:
    tuple path('*_processed_predictions.tsv'), model_pickle, bed, dataset into pp_sure_out
    
    script:
    template 'prediction_post_process_SuRE.R'
    

}

borzoi_out
.set {all_model_predictions_borzoi}

process prediction_post_process_borzoi {

    input:
    tuple prediction_file, model_pickle, bed, dataset from all_model_predictions_borzoi

    output:
    tuple path('*_processed_predictions.tsv'), model_pickle, bed, dataset into pp_borzoi_out
    
    script:
    template 'prediction_post_process_borzoi.R'
    

}

pp_sure_out
.mix(pp_borzoi_out)
.set {all_model_predictions}

process prediction_post_process {
    publishDir "${params.project_folder}/${bed.simpleName}/${dataset.simpleName}/${model_pickle.simpleName}/artifacts/${params.mode}/", mode: 'copy', saveAs: { filename -> "$filename" }

    input:
    tuple prediction_file, model_pickle, bed, dataset from all_model_predictions

    output:
    path '*.feather'
    
    script:
    template 'prediction_post_process.R'
    

}

process vep_tsv_like {

    publishDir "${params.project_folder}/${bed.simpleName}/${dataset.simpleName}/${model_pickle.simpleName}/artifacts/${params.mode}/", mode: 'copy', saveAs: { filename -> "$filename" }

    input:
    tuple tsv_file, model_pickle, bed, dataset from vep_in
    file reference_genome from file("$params.reference_genome")
    file vep_dir_cache from file("$params.vep_dir_cache")

    output:
    tuple path('*_vep_predictions.tsv.gz'), model_pickle, bed, dataset into vep_out
    
    when:
    model_pickle =~ /AlphaMissense/ && params.from_cache == 'NOT_PROVIDED'

    script:
    """
    vep --no_stats --offline --coding_only --protein --symbol --ccds --uniprot --canonical --merged \
    --dir_cache $vep_dir_cache \
    --fasta $reference_genome \
    --tab --fields 'Uploaded_variation,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,STRAND,SYMBOL,HGNC_ID,CANONICAL,CCDS,ENSP' \
    --compress_output gzip \
    --use_given_ref --gencode_basic \
    -i <(awk 'NR>1 {print \$1","\$3","\$3","\$4"/"\$5",+,"\$6}' $tsv_file | tr ',' '\t') \
    -o ${tsv_file.baseName}_vep_predictions.tsv.gz
    """

}

process alpha_missense {

    publishDir "${params.project_folder}/${bed.simpleName}/${dataset.simpleName}/${model_pickle.simpleName}/artifacts/${params.mode}/", mode: 'copy', saveAs: { filename -> "$filename" }

    input:
    tuple tsv_file, model_pickle, bed, dataset from vep_out

    output:
    path '*.feather'
    
    when:
    model_pickle =~ /AlphaMissense/

    script:
    template 'alpha_missense.R'

}

process alternative_allele_analysis {
publishDir "${params.project_folder}/${bed_file.simpleName}/${dataset.simpleName}/${model_pickle.simpleName}/plots/${params.mode}/alternativeAlleles/${cancer_type}/genes/", mode: 'copy', pattern: '*.pdf'
publishDir "${params.project_folder}/${bed_file.simpleName}/${dataset.simpleName}/${model_pickle.simpleName}/plots/${params.mode}/alternativeAlleles/${cancer_type}/", mode: 'copy', pattern: 'plotData*.tsv'

    input:
    tuple model_pickle, bed_file, dataset, chromosome, cancer_type from models_alternative_alleles_in.combine(subtypes_chromosome_selection_alternative_alleles)
    each gene_selection from subtypes_gene_selection_alternativeAlleles

    output:
    path '*.pdf'
    path '*_alternativeAlleles_out.tsv'
    path '*_alternativeAlleles.tsv' into alternative_alleles_out
    
    when:
    params.stage == 'run_analyses' && params.analyses =~ /alternativeAlleles/
    
    script:
    template 'alternative_allele_analysis.R'
    

}

process nearbyInPromoter_analysis {
publishDir "${params.project_folder}/${bed_file.simpleName}/${dataset.simpleName}/${model_pickle.simpleName}/plots/${params.mode}/nearbyInPromoter/${cancer_type}/genes/", mode: 'copy', pattern: '*.pdf'
publishDir "${params.project_folder}/${bed_file.simpleName}/${dataset.simpleName}/${model_pickle.simpleName}/plots/${params.mode}/nearbyInPromoter/${cancer_type}/", mode: 'copy', pattern: 'plotData*.tsv'

    input:
    tuple model_pickle, bed_file, dataset, chromosome, cancer_type from models_nearbyInPromoter_in.combine(subtypes_chromosome_selection_nearbyInPromoter)
    each gene_selection from subtypes_gene_selection_nearbyInPromoter

    output:
    path '*.pdf'
    path '*_nearbyInPromoter_out.tsv'
    path '*_nearbyInPromoter.tsv' into nearbyInPromoter_out
    
    when:
    params.stage == 'run_analyses' && params.analyses =~ /nearbyInPromoter/
    
    script:
    template 'nearbyInPromoter_analysis.R'
    

}

process samePromoter_analysis {
publishDir "${params.project_folder}/${bed_file.simpleName}/${dataset.simpleName}/${model_pickle.simpleName}/plots/${params.mode}/samePromoter/${cancer_type}/genes/", mode: 'copy', pattern: '*.pdf'
publishDir "${params.project_folder}/${bed_file.simpleName}/${dataset.simpleName}/${model_pickle.simpleName}/plots/${params.mode}/samePromoter/${cancer_type}/", mode: 'copy', pattern: 'plotData*.tsv'

    input:
    tuple model_pickle, bed_file, dataset, chromosome, cancer_type from models_samePromoter_in.combine(subtypes_chromosome_selection_samePromoter)
    each gene_selection from subtypes_gene_selection_samePromoter

    output:
    path '*.pdf'
    path '*_samePromoter_out.tsv'
    path '*_samePromoter.tsv' into same_promoter_out
    
    when:
    params.stage == 'run_analyses' && params.analyses =~ /samePromoter/
    
    script:
    template 'samePromoter_analysis.R'
    

}

subtypes_gene_selection_otherPromoters
.tap {subtypes_gene_selection_otherPromoters_allGenes}
.splitCsv(header: true)
.map {row -> tuple(row.associated_gene, row.chromosome, row.cancer_type)}
.filter { it[1] in "$params.chromosomes".tokenize( ',' ) }
.set {subtypes_gene_selection_otherPromoters_perGene}

process otherPromoters_analysis {

    input:
    tuple model_pickle, bed_file, dataset, associated_gene, chromosome, cancer_type from models_otherPromoters_in.combine(subtypes_gene_selection_otherPromoters_perGene)
    each gene_selection from subtypes_gene_selection_otherPromoters_allGenes

    output:
    path '*_otherPromoters.tsv' into otherPromoters_out
    
    when:
    params.stage == 'run_analyses' && params.analyses =~ /otherPromoters/
    
    script:
    template 'otherPromoters_analysis.R'
    

}

alternative_alleles_out
.mix(same_promoter_out)
.mix(otherPromoters_out)
.mix(nearbyInPromoter_out)
	.collectFile(keepHeader: true)
	.map{ it -> tuple(it.baseName, it)}
	.set{ alternative_alleles_all }
	
models_wgene
    .map {it -> it.baseName}
    .combine(cancer_types_wgene)
    .combine(analyses_wgene)
    .map {it -> tuple(it.join('_'), it[0], it[1], it[2])}
    .join(alternative_alleles_all)
    .combine(bedDataset_tissue)
    .set{all_grouped_analyses_tissue}
    
process plot_wilcoxon {
    publishDir "${params.project_folder}/${bed_file.simpleName}/${dataset.simpleName}/${model_pickle}/plots/${params.mode}/${analysis_type}/${cancer_type}/", mode: 'copy', saveAs: { filename -> "${params.mode}_${analysis_type}_${cancer_type}_vs_frequency.tsv" }, pattern: '*_table.tsv'
    publishDir "${params.project_folder}/${bed_file.simpleName}/${dataset.simpleName}/${model_pickle}/plots/${params.mode}/${analysis_type}/${cancer_type}/", mode: 'copy', saveAs: { filename -> "$filename" }, pattern: '*_vs_frequency.png'

    input:
    tuple _, model_pickle, cancer_type, analysis_type, gene_level_wilcoxon, bed_file, dataset from all_grouped_analyses_tissue
    file cancer_drivers from file("$params.cancer_drivers")
    
    output:
    path '*_table.tsv' into analyses_to_group_out
    path '*_vs_frequency.png'

    script:
    template 'plot_wilcoxon.R'

}


analyses_to_group_out
.collectFile(keepHeader: true)
.map{ it -> tuple(it.baseName, it)}
.set{analyses_all}

models_wAnalysis
    .map {it -> it.baseName}
    .combine(analyses_wAnalysis)
    .map {it -> tuple(it.join('_') + '_table', it[0], it[1])}
    .join(analyses_all)
    .combine(bedDataset_Analyses)
    .set{all_grouped_analyses}
    
process plot_wilcoxon_all_types {
publishDir "${params.project_folder}/${bed_file.simpleName}/${dataset.simpleName}/${model_pickle}/plots/${params.mode}/${analysis_type}/", mode: 'copy', saveAs: { filename -> "${model_pickle}_${analysis_type}_acrossTypes.tsv" }, pattern: '*_table.tsv'
    publishDir "${params.project_folder}/${bed_file.simpleName}/${dataset.simpleName}/${model_pickle}/plots/${params.mode}/${analysis_type}/", mode: 'copy', saveAs: { filename -> "$filename" }, pattern: '*_acrossTypes.png'
    
    input:
    tuple _, model_pickle, analysis_type, gene_level_wilcoxon, bed_file, dataset from all_grouped_analyses
    
    output:
    path '*_table.tsv' into models_to_group_out
    path '*_acrossTypes.png'
    
    script:
    template 'plot_wilcoxonAnalyses.R'

}

models_to_group_out
  .collectFile(keepHeader: true)
  .map{ it -> tuple(it.baseName, it)}
  .set{model_all}

models_wModel
    .map {it -> tuple("${it.baseName}_table", it.baseName)}
    .join(model_all)
    .combine(bedDataset_models)
    .set{all_grouped_models}
    
process plot_wilcoxon_model_results {
    publishDir "${params.project_folder}/${bed_file.simpleName}/${dataset.simpleName}/${model_pickle}/plots/${params.mode}/", mode: 'copy', saveAs: { filename -> "$filename"}
    input:
    tuple _, model_pickle, gene_level_wilcoxon, bed_file, dataset from all_grouped_models
    each cancer_drivers from Channel.fromPath("$params.cancer_drivers")
    each gene_annotation from Channel.fromPath("$params.gene_annotation")
    
    output:
    path 'summary_of_significant_genes.tsv'
    path 'summary_of_significant_gene_cancertype_hits.tsv'
    
    script:
    template 'plot_wilcoxonModel.R'

}

process plot_ism {
    publishDir "${params.project_folder}/${bed_file.simpleName}/${dataset.simpleName}/${model_pickle}/plots/${params.mode}/genes/", mode: 'copy', saveAs: { filename -> "$filename"}
    input:
    tuple model_pickle, bed_file, dataset from models_plotISM_in
    each gencode_annotation from Channel.fromPath("$params.gencode_annotation")
    path gene_selection from gene_selection_plot
    
    output:
    path 'ism.pdf'
    
    when:
    params.stage == "plot_genes"
    
    script:
    template 'plot_ism.R'

}
