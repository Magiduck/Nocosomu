process name_file {
  
  input:
  tuple val(fileID), path(some_file) 
  
  output:
  path "${fileID}_${some_file}"
  
  script:
  """
  cp ${some_file} ${fileID}_${some_file}
  """
  
}

process name_file2 {
  
  input:
  tuple val(fileID), val(groupID), path(some_file) 
  
  output:
  path "${fileID}_${groupID}_${some_file}"
  
  script:
  """
  cp ${some_file} ${fileID}_${groupID}_${some_file}
  """
  
}

process publish_file {
  
  publishDir "${params.outdir}/$sectionID", mode: 'copy', saveAs: { filename -> "$filename" }
  
  input:
    path tsv_file
    val sectionID
  
  output:
    path tsv_file
  
  script:
    """
  """
  
}

process publish_file2 {
  
  publishDir "${params.outdir}/${sectionID}", mode: 'copy', saveAs: { filename -> "$filename" }
  
  input:
    tuple val(sectionID), path(tsv_file)
  
  output:
    path tsv_file
  
  script:
  """
  echo working on $tsv_file
  """
  
}

process collectFileList {
  
  input:
  tuple val(ID), path(dummy) 
  
  output:
  tuple val(ID), path('merged_*', arity: '1')
  
  script:
  """
  cat * > merged_${dummy[0]}
  """
  
}

process delimitedFileSplitter {
  
  input:
  tuple val(ID), path(big_file) 
  val nlines
  
  output:
  tuple val(ID), path('piece*')
  
  script:
  template 'split_delimited_file.R'
  
}
