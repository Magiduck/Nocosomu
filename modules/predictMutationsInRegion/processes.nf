
process variants_in_region {
	input:
	path bed_file

	output:
	tuple val("${bed_file.simpleName}"), path('variants_in_region.tsv')

	script:
	template 'bed2tsv.R'
	
}

