include { Borzoi_predict} from './../borzoi'
include { SureResNet_predict} from './../sureResnet'
include { no_background; alternativeAlleles_background; samePromoter_background; otherTumors_background; region_metadata} from './../backgroundModels'
include { transform_to_vcflike; wilcoxon_compare; prediction_post_process_alphagenome; prediction_post_process_borzoi; prediction_post_process_sure; add_trinucleotide; region_subselection; mutation_subselection; agreggate_analysis_files } from './processes.nf'
include { publish_file2 as publish_summary } from './../general'
include { Models_predict } from './common.nf'
include { SimpleBackground as SimpleBackground_alternativeAlleles } from './common.nf'
include { SimpleBackground as SimpleBackground_samePromoter } from './common.nf'
include { SimpleBackground as SimpleBackground_otherTumors } from './common.nf'

workflow {
Analyses_predict()

Analyses_predict.out
| map {tuple("/${params.bedID}/${params.datasetID}/${it[0]}/${it[2]}/${it[1]}", it[3])}
| publish_summary
}

workflow Analyses_predict {
   
   main:
   chromosomes = Channel.fromList("$params.chromosomes".tokenize( ',' ))
   .map {tuple(it, "${params.outdir}/${params.bedID}/${params.datasetID}/$it/elegible_regions.csv")}
   .filter { file(it[1]).exists() } //If no elegible regions are found then drop the chromosome
   .map { it[0] }
   
   //Determine regions to run 
   elegible_regions = chromosomes
      .map {tuple(it, "${params.outdir}/${params.bedID}/${params.datasetID}/$it/elegible_regions.csv")}

   region_selection = region_subselection(elegible_regions, "${params.region_subselection}")
   | filter { it[1].countLines() > 1 } // no eligible mutations are found then drop
   
   //Determine saturation mutations
   mut_file = chromosomes
      .map {tuple(it, "${params.outdir}/${params.bedID}/${params.datasetID}/$it/merged_mutations_in_region.tsv")}

   ism_file = region_metadata(mut_file.join(region_selection))
   ism_file_wtrinucleotide =  add_trinucleotide(ism_file, "$params.reference_genome")
   
   //Limit mutations to subselection
   mutation_files_selection = chromosomes
      .map {tuple(it, "${params.outdir}/${params.bedID}/${params.datasetID}/$it/merged_mutations_in_region.tsv")}
      .join(region_selection)
      | mutation_subselection
      | filter { it[1].countLines() > 1 } //If no elegible mutations are found then drop the file
   mutation_files = mutation_files_selection.join(ism_file_wtrinucleotide)
      
   
   //Make backgrounds
   no_background_sbs_mutations = no_background(mutation_files)
   .filter { it[1].countLines() > 1 } //If no elegible mutations are found then drop the file
   .map { tuple('noBackground' + '_' + it[0], it[1]) }
   realMutation_predictions = Models_predict(no_background_sbs_mutations)
   
   alternativeAlleles_analysis_files = Channel.empty()
   if ("$params.analyses" =~ /alternativeAlleles/){
       alternativeAlleles_background(mutation_files)
       alternativeAlleles_analysis_files = SimpleBackground_alternativeAlleles(alternativeAlleles_background.out.background,
                                                            alternativeAlleles_background.out.donorMapping,
                                                           'alternativeAlleles',
                                                           realMutation_predictions,
                                                           region_selection
                                                           )
   }
   
   samePromoter_analysis_files = Channel.empty()
   if ("$params.analyses" =~ /samePromoter/){
       samePromoter_background(mutation_files)
       samePromoter_analysis_files = SimpleBackground_samePromoter(samePromoter_background.out.background,
                                                      samePromoter_background.out.donorMapping,
                                                           'samePromoter',
                                                           realMutation_predictions,
                                                           region_selection
                                                           )
   }
   
   otherTumors_analysis_files = Channel.empty()
   if ("$params.analyses" =~ /otherTumors/){
       otherTumors_background(mutation_files)
       otherTumors_analysis_files = SimpleBackground_otherTumors(otherTumors_background.out.background,
                                                      otherTumors_background.out.donorMapping,
                                                           'otherTumors',
                                                           realMutation_predictions,
                                                           region_selection
                                                           )
   }

   analysis_files = alternativeAlleles_analysis_files
      .mix(samePromoter_analysis_files)
      .mix(otherTumors_analysis_files)
    
    analysis_files_per_chromosome = analysis_files.combine(region_selection, by: 0)
    analysis_files_aggregation = agreggate_analysis_files(analysis_files.combine(region_selection, by: 0)
                                                          .groupTuple(by:[1,2]), "$params.region_aggregation")
    //Run comparisons
    
    all_analysis_files = analysis_files_per_chromosome
      .mix(analysis_files_aggregation)
    wilcoxon_compare(all_analysis_files, "$params.minimum_donors")
   
   prediction_files = wilcoxon_compare.out.summary
   .mix(wilcoxon_compare.out.summary_abs)
   .mix(wilcoxon_compare.out.tested_mutations)
   
   emit:
   prediction_files

}

