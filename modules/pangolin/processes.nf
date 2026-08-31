process pangolin_tsv_like {
  label 'gpu'
  
  input:
  tuple val(ID), path(tsv_file)
  each path(pangolin_annotation)
  each path(reference_genome)
  each pangolin_params

  output:
  tuple val(ID), val('pangolin'), path("prediction.csv")
  
  script:
    """
        cat $tsv_file | tr '\t' ',' > input.csv
        pangolin $pangolin_params\
        input.csv \
        $reference_genome \
        $pangolin_annotation \
        prediction
    """
  
}

process pangolin_merge {
  label 'rscript'
  
  input:
  tuple val(ID), val(groupID), path('?_predictions.tsv')
  
  output:
  tuple val(ID), val(groupID), path("*Splice.tsv", arity: '3')
  
  script:
    template 'merge_pangolin_data.R'
  
}

process pangolin_format {
  
  input:
  tuple val(ID), path(vcf_file)
  
  output:
  tuple val(ID), path('snvs.tsv')
  
  shell:
  '''
  ml BCFtools
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t\n' !{vcf_file} | awk '{ if ($3~/^.$/ && $4~/^.$/) { print } }' | cut -f1-4 | awk '{ printf $1"_"$2"_"$3"_"$4"\\t"$0"\\n"}' > snvs.tsv
  '''
  
}
