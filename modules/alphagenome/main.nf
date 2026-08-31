include { ALPHAGENOME_PREDICT; alphagenome_format; MERGE_ALPHAGENOME_CHUNKS } from './processes.nf'
include { name_file2; publish_file; delimitedFileSplitter } from './../general'

workflow Alphagenome_predict {
    take:
        vcflikeFiles

    main:
        // Format the generic vcflike files specifically for AlphaGenome
        alphagenome_format(vcflikeFiles)

        // Split the formatted files into smaller chunks
        vcflikeFileShards_ch = delimitedFileSplitter(alphagenome_format.out, params.chunkSize)
        | transpose

        // Pass the delimited mutation files (shards) to the prediction process
        ALPHAGENOME_PREDICT(vcflikeFileShards_ch)

        // Collect the resulting CSV outputs
        // Group all predicted chunks by the original ID (0) so they can be merged.
        // Map a dummy 'alphagenome' groupID because AlphaGenome does not use model folds.
        predictions_grouped = ALPHAGENOME_PREDICT.out.predictions
        | groupTuple(by: 0)
        | map { tuple(it[0], "alphagenome", it[1]) } 

        // Merge the collected chunks into a single consolidated file per group
        predictions = MERGE_ALPHAGENOME_CHUNKS(predictions_grouped)

    emit:
        predictions
}
