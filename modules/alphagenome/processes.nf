process ALPHAGENOME_PREDICT {
  label 'gpu'
  tag "${tsv_file.baseName}"

  input:
  tuple val(ID), path(tsv_file)

  output:
  tuple val(ID), path("${tsv_file.baseName}_alphagenome_scores.csv"), emit: predictions

  script:
  """
  # Call the python script with the explicit module directory path
  ${moduleDir}/bin/alphagenome_from_vcf.py \\
      ${tsv_file} \\
      ${tsv_file.baseName}_alphagenome_scores.csv
  """
}

process alphagenome_format {
  label 'rscript'
  tag "${vcflikeFile.baseName}"

  input:
  tuple val(ID), path(vcflikeFile)

  output:
  tuple val(ID), path("${vcflikeFile.baseName}_alphagenome_format.tsv")

  script:
  template 'format_alphagenome.R'
}

process MERGE_ALPHAGENOME_CHUNKS {
  label 'rscript'
  tag "${groupID}"
  
  input:
  tuple val(ID), val(groupID), path(chunks)

  output:
  tuple val(ID), val(groupID), path("${groupID}_merged.csv")

  script:
  """
  #!/usr/bin/env Rscript
  library(data.table)
  
  all_files = list.files()
  csv_files = all_files[endsWith(all_files, ".csv")]
  
  dts = lapply(csv_files, function(x) {
      if(file.info(x)[["size"]] > 0) {
          return(fread(x))
      }
      return(NULL)
  })
  
  valid_dts = Filter(Negate(is.null), dts)
  
  if(length(valid_dts) > 0) {
      dt = rbindlist(valid_dts, fill=TRUE)
      fwrite(dt, "${groupID}_merged.csv", sep=",", quote=FALSE)
  } else {
      dt_empty = data.table(variant_name=character(), raw_score=numeric(), gene_name=character())
      fwrite(dt_empty, "${groupID}_merged.csv", sep=",", quote=FALSE)
  }
  """
}
