include { cadd_tsv_like; cadd_merge_folds } from './processes.nf'
include { name_file2; publish_file; delimitedFileSplitter } from './../general'

workflow {
  Cadd_predict(Channel.fromPath("$params.vcflike").map {tuple(it.baseName, it)})
  name_file2(Cadd_predict.out)
  publish_file(name_file2.out, 'borzoi_out')
}

workflow Cadd_predict {
  take:
    vcflikeFiles
    caddModelPath_ch
  
  main:
    vcflikeFileShards_ch = delimitedFileSplitter(vcflikeFiles, 1500)
    | transpose
    
    vcflikeFile_results = cadd_tsv_like(vcflikeFileShards_ch, caddModelPath_ch)
    | groupTuple(by:[0,1]) // ID groupID file (based on file name) -> ID groupID list of all matching pieces
    
    vcflikeFile_merged = cadd_merge_folds(vcflikeFile_results)

    emit:
    vcflikeFile_merged
}

