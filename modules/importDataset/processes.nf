process import_ICGC {
  
  input:
  tuple path(parquet_file), path(bed_file), val(chromosome)
  path(subgroup_file)

  output:
  tuple val("${bed_file.simpleName}"), val("${parquet_file.simpleName}"), val(chromosome), path("merged_mutations_in_region.tsv")

  script:
    template 'add_region_icgc.R'
  
}

process import_GE {
  input:
  tuple path(tsv_file), path(bed_file), val(chromosome)
  path(subgroup_file)

  output:
  tuple val("${bed_file.simpleName}"), val("${tsv_file.simpleName}"), val(chromosome), path("merged_mutations_in_region.tsv")

  script:
    """
    prepare_tsv_for_bed.py -i $tsv_file -c $chromosome -o gel_all_pseudodonor.bed 
    bedtools intersect -a $bed_file -b gel_all_pseudodonor.bed -wo | cut -f2,3,4,10 | sort | uniq > mutations_ids_in_region.tsv
    add_region_ge.py -i $tsv_file -r mutations_ids_in_region.tsv -c $chromosome -m $bed_file -o merged_mutations_in_region.tsv
    """ 

}

process remove_germline {
  label 'pyscript'

  input:
  tuple val(bed_id), val(parquet_id), val(chromosome), path(merged_mutations_in_regions), path(germline_variants), val(threshold)
  
  output:
  tuple val(bed_id), val(parquet_id), val(chromosome), path(merged_mutations_in_regions)

  script:
    """
    mv $merged_mutations_in_regions merged_mutations_in_region_prefiltered.tsv
    remove_germline_icgc.py -i merged_mutations_in_region_prefiltered.tsv -g $germline_variants -o merged_mutations_in_region.tsv -a $threshold
    """
}

process filter_basepairs {
  label 'pyscript'

  input:
  tuple val(bed_id), val(parquet_id), val(chromosome), path(merged_mutations_in_regions, stageAs:'merged_mutations_in_regions_prefiltered_before_coverage.tsv'), path(below_coverage), val(skip_qc)

  output:
  tuple val(bed_id), val(parquet_id), val(chromosome), path("merged_mutations_in_region.tsv")

  script:
    if( skip_qc == false )
    """
    remove_below_coverage.py -i $merged_mutations_in_regions -c $below_coverage -o merged_mutations_in_region.tsv
    """
    
    else
    """
    cp $merged_mutations_in_regions merged_mutations_in_region.tsv
    """
}

process keep_basepairs {
  label 'pyscript'

  input:
  tuple val(bed_id), val(parquet_id), val(chromosome), path(merged_mutations_in_regions, stageAs:'merged_mutations_in_regions_prefiltered_before_coverage.tsv'), path(below_coverage), val(skip_qc)

  output:
  tuple val(bed_id), val(parquet_id), val(chromosome), path("merged_mutations_in_region.tsv")

  script:
    if( skip_qc == false )
    """
    whitelist_basepairs.py -i $merged_mutations_in_regions -c $below_coverage -o merged_mutations_in_region.tsv
    """
    
    else
    """
    cp $merged_mutations_in_regions merged_mutations_in_region.tsv
    """
}

process filter_on_blacklist {
  label 'pyscript'

  input:
  path(blacklist)
  path(bed_file)

  output:
  path("${bed_file.baseName}")

  script:
    """
    cp $bed_file regions_before_blacklist.bed
    bedtools intersect -a regions_before_blacklist.bed -b $blacklist -v -header > ${bed_file.baseName}
    """
}

process elegible {
  
  input:
  tuple path(parquet_file), path(bed_file), val(chromosome), path(merged_mutations)
  path(subgroup_file)
  each minimum_donors

  output:
  tuple val(bed_file.simpleName), val(parquet_file.simpleName), val(chromosome), path('elegible_regions.csv')
  
  script:
    template 'elegible_ICGC.R'
}
