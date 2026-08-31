include { Borzoi_predict} from './../borzoi'
include { SureResNet_predict} from './../sureResnet'
include { name_file2; publish_file2 } from './../general'

include { createISM } from './processes.nf'


workflow {
  Region_predict(Channel.fromPath("$params.bed").map {tuple(it.baseName, it)})
  //name_file2(Region_predict.out)
  //publish_file2(name_file2.out, 'borzoi_out')
}

workflow Region_predict {

   take:
   bed_file

   main:
   
   createISM(bed_file)
   
   prediction_files = Borzoi_predict(createISM.out.map {tuple(it.baseName, it)})
   
   emit:
   prediction_files

}
