//include { SureResNet_predict} from './../sureResnet'
include { Borzoi_predict} from './../borzoi'

workflow Models_predict {
   take:
   vcflikeFiles
   
   main:
   variantFiles_w_models = vcflikeFiles
   
   Borzoi_results = Borzoi_predict(variantFiles_w_models)
    .map { tuple(it[0] + '_' + 'borzoi', it[1]) }
   
   //SureResNetK562_results = SureResNet_predict(variantFiles_w_models, "locationofmodel")
   //.map { tuple('SureResNetK562', it) }
   //SureResNetHEPG2_results = SureResNet_predict(variantFiles_w_models, "locationofmodel")
     
   all_results = Borzoi_results
   //.mix(SureResNetK562_results)
   //.mix(SureResNetHEPG2_results)
   
   emit:
   all_results

}

