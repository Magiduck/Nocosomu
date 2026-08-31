include { borzoi_format; borzoi_tsv_like; borzoi_extract_sed; borzoi_merge_folds; borzoi_splice_effect} from './processes.nf'
include { name_file2; publish_file; delimitedFileSplitter } from './../general'

workflow {
  if("$params.borzoiMode" == "splice"){
     Borzoi_predict(Channel.fromPath("$params.vcflike").map {tuple(it.baseName, it)}.filter { it[1].countLines() > 1 })
     splice_effect = borzoi_splice_effect(Borzoi_predict.out, "$params.borzoi_conversionTable")
     name_file2(splice_effect)
     publish_file(name_file2.out, 'borzoi_out_splice')
  } else{
     Borzoi_predict(Channel.fromPath("$params.vcflike").map {tuple(it.baseName, it)}.filter { it[1].countLines() > 1 })
     name_file2(Borzoi_predict.out)
     publish_file(name_file2.out, 'borzoi_out')
  }
}

workflow Borzoi_predict {
  take:
    vcflikeFiles
  
  main:
    borzoi_format(vcflikeFiles, params.borzoi_columns)
    
    vcflikeFileShards_ch = delimitedFileSplitter(borzoi_format.out, params.chunkSize)
    | transpose
      
    borzoi_models = Channel.fromList("$params.borzoi_folds".tokenize( ',' ))
      .map {"${params.borzoi_dir}/saved_models/${it}/model0_best.h5"}
//    borzoi_models = Channel.fromList("$params.borzoi_folds".tokenize( ',' ))
//      .map {"${params.borzoi_dir}/mini_models/human_gtex/${it}/model0_best.h5"}

   borzoi_tsv_like(vcflikeFileShards_ch,
                  borzoi_models,
                  "$params.reference_genome", 
                  "${params.borzoi_targets}",
                  "${params.borzoi_annotation}",
                  "${params.borzoi_dir}/params_pred.json",
                  "${params.borzoi_flags}")
    
//                  "${params.borzoi_dir}/mini_models/human_gtex/params.json",

    vcflikeFile_results = borzoi_extract_sed(borzoi_tsv_like.out,
                                            "${params.borzoi_targets}")
    | transpose // ID with list of all borzoi results ->> ID with single borzoi result (flattened)
    | map {tuple(it[0], it[1].simpleName, it[1])} // ID with single borzoi result = ID groupID file (based on file name)
    | groupTuple(by:[0,1]) // ID groupID file (based on file name) -> ID groupID list of all matching pieces

    
    vcflikeFile_merged = borzoi_merge_folds(vcflikeFile_results,
                                            "$params.borzoi_folds")
                                          
    emit:
    vcflikeFile_merged
}

