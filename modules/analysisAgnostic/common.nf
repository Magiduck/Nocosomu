include { Borzoi_predict} from './../borzoi'
include { SureResNet_predict} from './../sureResnet'
include { Alphagenome_predict} from './../alphagenome'
include { Cadd_predict} from './../cadd'
include { transform_to_vcflike; prediction_post_process_alphagenome; prediction_post_process_borzoi; prediction_post_process_sure; prediction_post_process_cadd } from './processes.nf'
include { SureResNet_predict as SureResNet_predict_k562 } from './../sureResnet'
include { SureResNet_predict as SureResNet_predict_hepg2 } from './../sureResnet'
include { SureResNet_predict as SureResNet_predict_caco2 } from './../sureResnet'
include { SureResNet_predict as SureResNet_predict_colo320 } from './../sureResnet'
include { SureResNet_predict as SureResNet_predict_ht29 } from './../sureResnet'
include { SureResNet_predict as SureResNet_predict_ito66t } from './../sureResnet'
include { SureResNet_predict as SureResNet_predict_hct116 } from './../sureResnet'
include { SureResNet_predict as SureResNet_predict_hct116_mh04 } from './../sureResnet'
include { prediction_post_process_sure as prediction_post_process_sure_k562 } from './processes.nf'
include { prediction_post_process_sure as prediction_post_process_sure_hepg2 } from './processes.nf'
include { prediction_post_process_sure as prediction_post_process_sure_caco2 } from './processes.nf'
include { prediction_post_process_sure as prediction_post_process_sure_colo320 } from './processes.nf'
include { prediction_post_process_sure as prediction_post_process_sure_ht29 } from './processes.nf'
include { prediction_post_process_sure as prediction_post_process_sure_ito66t } from './processes.nf'
include { prediction_post_process_sure as prediction_post_process_sure_hct116 } from './processes.nf'
include { prediction_post_process_sure as prediction_post_process_sure_hct116_mh04 } from './processes.nf'

include { filter_basepairs as filter_repeats; filter_basepairs as filter_lowComplexity; filter_basepairs as filter_on_coverage; filter_basepairs as filter_on_mappability; filter_basepairs as filter_coding; remove_germline} from './../importDataset/processes.nf'

workflow Models_predict {

   take:
   all_results

   main:
   vcflike = transform_to_vcflike(all_results)
   
   parmhepg2_results = Channel.empty()
   if ("$params.model" =~ /PARM-HEPG2/){
    parmhepg2_results = prediction_post_process_sure_hepg2(SureResNet_predict_hepg2(vcflike, Channel.fromPath("$params.parm_hepg2")))
   }
   parmk562_results = Channel.empty()
   if ("$params.model" =~ /PARM-K562/){
    parmk562_results = prediction_post_process_sure_k562(SureResNet_predict_k562(vcflike, Channel.fromPath("$params.parm_k562")))
   }
   parmcaco2_results = Channel.empty()
   if ("$params.model" =~ /PARM-CACO2/){
    parmcaco2_results = prediction_post_process_sure_caco2(SureResNet_predict_caco2(vcflike, Channel.fromPath("$params.parm_caco2")))
   }
   parmcolo320_results = Channel.empty()
   if ("$params.model" =~ /PARM-COLO320/){
    parmcolo320_results = prediction_post_process_sure_colo320(SureResNet_predict_colo320(vcflike, Channel.fromPath("$params.parm_colo320")))
   }
   parmht29_results = Channel.empty()
   if ("$params.model" =~ /PARM-HT29/){
    parmht29_results = prediction_post_process_sure_ht29(SureResNet_predict_ht29(vcflike, Channel.fromPath("$params.parm_ht29")))
   }
   parmito66t_results = Channel.empty()
   if ("$params.model" =~ /PARM-ITO66T/){
    parmito66t_results = prediction_post_process_sure_ito66t(SureResNet_predict_ito66t(vcflike, Channel.fromPath("$params.parm_ito66t")))
   }
   parmhct116_results = Channel.empty()
   if ("$params.model" =~ /PARM-HCT116/){
    parmhct116_results = prediction_post_process_sure_hct116(SureResNet_predict_hct116(vcflike, Channel.fromPath("$params.parm_hct116")))
   }
   parmhct116_mh04_results = Channel.empty()
   if ("$params.model" =~ /PARM-HCT116-MH04/){
    parmhct116_mh04_results = prediction_post_process_sure_hct116_mh04(SureResNet_predict_hct116_mh04(vcflike, Channel.fromPath("$params.parm_hct116_mh04")))
   }
   borzoi_results = Channel.empty()
   if ("$params.model" =~ /BORZOI/){
    borzoi_results = prediction_post_process_borzoi(Borzoi_predict(vcflike))
   .map { tuple(it[0], 'borzoi_' + it[1], it[2])}
   }
   cadd_results = Channel.empty()
   if ("$params.model" =~ /CADD/){
    cadd_results = prediction_post_process_cadd(Cadd_predict(vcflike, Channel.fromPath("$params.cadd_file")))
   }
   alphagenome_results = Channel.empty()
   if ("$params.model" =~ /ALPHAGENOME/){
    alphagenome_results = prediction_post_process_alphagenome(Alphagenome_predict(vcflike), Channel.fromPath("$params.gene_annotation"))
   .map { tuple(it[0], 'alphagenome_' + it[1], it[2])}
   }
   
   prediction_files = alphagenome_results
   .mix(borzoi_results)
   .mix(parmhepg2_results)
   .mix(parmk562_results)
   .mix(parmcaco2_results)
   .mix(parmcolo320_results)
   .mix(parmht29_results)
   .mix(parmito66t_results)
   .mix(parmhct116_results)
   .mix(parmhct116_mh04_results)
   .mix(cadd_results)
   .map { tuple(it[0].split('_')[0], it[0].split('_')[1], it[1], it[2])} // label chr model file
   
   emit:
   prediction_files

}

workflow SimpleBackground {

   take:
   background_test
   donorMapping
   background_test_id
   realMutation_predictions
   region_selection
   
   main:
   
   //Do QC, then if no elegible mutations remain then drop the file
   background_results = background_test 
   | map { tuple('_dummyBED', '_dummyParquet', it[0], it[1], "${params.below_coverage}/chr${it[0]}_positions_below_30.csv.gz", "$params.skipQC" ==~ /COV30/) }
   | filter_on_coverage
   | map { tuple(it[0], it[1], it[2], it[3], "${params.within_repeat}/chr${it[2]}_harmonizedDietlein.txt.gz", "$params.skipQC" ==~ /REPEAT/) }
   | filter_repeats
   | map { tuple(it[0], it[1], it[2], it[3], "${params.within_repeat}/chr${it[2]}_repeatMaskerLowComplexity.txt.gz", "$params.skipQC" ==~ /COMPLEXITY/ ) }
   | filter_lowComplexity
   | map { tuple(it[0], it[1], it[2], it[3], "${params.low_mappable}/chr${it[2]}_hg38.k36.umap.multiread_blacklist.txt.gz", "$params.skipQC" ==~ /UMAPK36/ ) }
   | filter_on_mappability
   | map { tuple(it[0], it[1], it[2], it[3], "${params.coding}/chr${it[2]}_CDS.txt.gz", "$params.skipQC" ==~ /CDS/ ) }
   | filter_coding
   | map { tuple(it[0], it[1], it[2], it[3], "${params.germline_variants}/gnomad.joint.v4.1.sites.chr${it[2]}_FAF.tsv.gz", params.germlineAlleleFrequency) }
   | remove_germline
   | map { tuple(it[2], it[3]) }
   | filter { it[1].countLines() > 1 }
   
   all_results = background_results.map { tuple(background_test_id + '_' + it[0], it[1]) }
     
   prediction_files = Models_predict(all_results)
   
   analysis_files = prediction_files
    .join(realMutation_predictions, by: [1,2])
    .combine(background_results, by:0)
    .combine(donorMapping, by:0)
    
   //analysis_files = prediction_files
  //  .join(realMutation_predictions, by: [1,2])
  //  .combine(region_selection, by: 0)
  //  .combine(background_results, by:0)
  //  .combine(donorMapping, by:0)
    
    //agreggated_analysis_files = analysis_files.groupTuple(by:[1,2]).view()

   emit:
   analysis_files

}
