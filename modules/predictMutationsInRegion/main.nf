include { Models_predict} from './../modelAgnostic'

//borzoi_models = Channel.fromList("$params.borzoi_folds".tokenize( ',' ))

workflow PredictRegion {
   take:
   bedFile
   
   main:
   variantFiles = variants_in_region(bedFile)
   Models_predict(variantFiles, Channel.fromList("params.models".tokenize(',')))
   
   
   emit:
   Models_predict.out

}

