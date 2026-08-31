include { pangolin_tsv_like; pangolin_merge; pangolin_format } from './processes.nf'
include { name_file2; publish_file; delimitedFileSplitter } from './../general'

workflow {
  Pangolin_predict(Channel.fromPath("$params.vcflike").map {tuple(it.simpleName, it)}, Channel.fromPath("$params.pangolin_annotation"))
  name_file2(Pangolin_predict.out)
  publish_file(name_file2.out, 'Pangolin_out')
}

workflow Pangolin_predict {
  take:
    vcflikeFiles_ch
    pangolinAnnotation_ch
  
  main:
    if("$params.fmt" == "vcf"){
      formatted_file = pangolin_format(vcflikeFiles_ch)
    } else{
      formatted_file = vcflikeFiles_ch
    }
  
    vcflikeFileShards_ch = delimitedFileSplitter(formatted_file, 15000)
    | transpose
    | filter { it[1].countLines() > 1 }

    vcflikeFile_results = pangolin_tsv_like(vcflikeFileShards_ch,
                                              pangolinAnnotation_ch,
                                              "$params.reference_genome",
                                              "$params.pangolin_params")
    | groupTuple(by:[0,1]) // ID groupID file (based on file name) -> ID groupID list of all matching pieces
    
    vcflikeFile_merged = pangolin_merge(vcflikeFile_results)
    | transpose // ID groupID with list of all parm strand results ->> ID groupID with single parm strand prediction (flattened)
    | map { tuple(it[0], it[1] + '_' + it[2].simpleName, it[2]) } // ID groupID with single parm strand prediction = ID groupID+strand file (based on file name)

  emit:
    vcflikeFile_merged
}

