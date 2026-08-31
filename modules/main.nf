include { Borzoi_predict} from './../borzoi'
include { SureResNet_predict} from './../sureResNet'
include { name_file; publish_file } from './../general'


workflow {
  SureResNet_predict("$params.vcflike")
  publish_file(sureResNet_predict.out, 'sureResnet_out')
}