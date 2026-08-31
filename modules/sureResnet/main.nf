include { sureResNet_tsv_like; parm_merge } from './processes.nf'
include { name_file2; publish_file; delimitedFileSplitter } from './../general'

workflow {
  SureResNet_predict(Channel.fromPath("$params.vcflike").map {tuple(it.simpleName, it)}, Channel.fromPath("$params.sure_model"))
  name_file2(SureResNet_predict.out)
  publish_file(name_file.out, 'sureResnet_out')
}

workflow SureResNet_predict {
  take:
    vcflikeFiles_ch
    sureModelPath_ch
  
  main:
    vcflikeFileShards_ch = delimitedFileSplitter(vcflikeFiles_ch, 15000)
    | transpose

    vcflikeFile_results = sureResNet_tsv_like(vcflikeFileShards_ch,
                                              sureModelPath_ch,
                                              "$params.reference_genome")
    | groupTuple(by:[0,1]) // ID groupID file (based on file name) -> ID groupID list of all matching pieces
    
    vcflikeFile_merged = parm_merge(vcflikeFile_results)
    | transpose // ID groupID with list of all parm strand results ->> ID groupID with single parm strand prediction (flattened)
    | map { tuple(it[0], it[1] + '_' + it[2].simpleName, it[2]) } // ID groupID with single parm strand prediction = ID groupID+strand file (based on file name)

  emit:
    vcflikeFile_merged
}

