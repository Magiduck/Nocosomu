include { name_file; name_file2; publish_file; collectFileList; publish_file2; delimitedFileSplitter} from './processes.nf'

workflow CollectFileTuple {
  
  take:
  Files_with_IDs
  
  main:
  groupedFiles = Files_with_IDs
  | groupTuple
  collectedFiles = collectFileList(groupedFiles)
  
  emit:
  collectedFiles
  
}
