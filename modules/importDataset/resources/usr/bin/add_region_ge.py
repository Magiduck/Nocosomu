#!/usr/bin/env python3
"""
<A single line describing this program goes here.>

MIT License

Copyright (c) 2022 Tijs van Lieshout

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Uses:
#TODO: remove
singularity \
exec \
-B \
/groups/umcg-fg/tmp02/projects/non-coding-somatic/ \
/groups/umcg-fg/tmp02/projects/non-coding-somatic/shared/singularity/containers_cgut/sequence_based_models.sif \
python3 /groups/umcg-fg/tmp02/projects/non-coding-somatic/umcg-tvanlieshout/repos_tvanlieshout/nocosomu_iDriver/modules/importDataset/bin/add_region_ge.py \
-i /groups/umcg-fg/tmp02/projects/non-coding-somatic/cancer/genomics_england/hmf_regions_frequency_matrices/gel_all_pseudodonor.tsv.gz \
-r mutations_ids_in_region.tsv \
-c 5 \
-m /groups/umcg-fg/tmp01/projects/non-coding-somatic/shared/reftss/artifacts/2024-05-03_refTSS_v4.1_human_coordinate_srtdb_overlap_for_cutoff_22_900_upstream_300_downstream_1bpTSS_sorted_nonoverlapping.bed \
-o TMP

singularity \
exec \
-B \
/groups/umcg-fg/tmp01/projects/non-coding-somatic/ \
/groups/umcg-fg/tmp01/projects/non-coding-somatic/shared/singularity/containers_cgut/sequence_based_models.sif \
python3 /groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-tvanlieshout/repos_tvanlieshout/nocosomu_iDriver/modules/importDataset/bin/add_region_ge.py \
-i /groups/umcg-fg/tmp01/projects/non-coding-somatic/cancer/genomics_england/hmf_regions_frequency_matrices/gel_all_pseudodonor.tsv.gz \
-r mutations_ids_in_region.tsv \
-c 5 \
-m /groups/umcg-fg/tmp01/projects/non-coding-somatic/shared/reftss/artifacts/2024-05-03_refTSS_v4.1_human_coordinate_srtdb_overlap_for_cutoff_22_900_upstream_300_downstream_1bpTSS_sorted_nonoverlapping.bed \
-o TMP
"""

# Metadata
__title__ = "Template for a CLI python script" 
__author__ = "Tijs van Lieshout"
__created__ = "2022-05-04"
__updated__ = "2022-05-27"
__maintainer__ = "Tijs van Lieshout"
__email__ = "t.van.lieshout@umcg.nl"
__version__ = 0.2
__license__ = "GPLv3"
__description__ = f"""{__title__} is a python script created on {__created__} by {__author__}.
                      Last update (version {__version__}) was on {__updated__} by {__maintainer__}.
                      Under license {__license__} please contact {__email__} for any questions."""

# Imports
import argparse
import os

import pandas as pd

def main(args):
  #TODO: output with same columns
  #TODO: all_subtypes (duplicate)
  df = pd.read_csv(args.inputPath, sep="\t")
  df['mutation_id'] = df['#CHROM'] + "_" + df["POS"].astype(str) + "_" + df["REF"] + "_" + df["ALT"]
  allowed_mutations = pd.read_csv(args.regionsPath, sep="\t", header=None, names=['min_pos_region', 'max_pos_region', 'associated_region', 'mutation_id'])
  df['chromosome'] = df['#CHROM'].str.replace("chr", "")
  df = df[df['chromosome'] == args.chromosome]
  df = pd.merge(df, allowed_mutations, how='inner', on='mutation_id')
  df['chromosome_start'] = df["POS"] - 1
  df['chromosome_end'] = df["POS"]
  df['strand'] = "." #TODO: discuss if we need strand information from genomics england
  df['reference_genome_allele'] = "" #TODO: TMP, probably not used in pipeline...
  df['mutated_from_allele'] = df["REF"]
  df['mutated_to_allele'] = df["ALT"]
  df = annotate_mutation_type(df)
  df = df[df['mutation_type'] == "single base substitution"] #TODO: TMP in order to recreate the way the ICGC script works, but should be handled in different part of the pipeline
  df['subgroup'] = df['gel_cancer'].str.replace("_frequency_table.tsv", "")
  df['region_metadata'] = args.metadata
  df = df[['mutation_id', 'donor_id', 'chromosome', 'chromosome_start', 'chromosome_end', 'strand', 'mutation_type', 'reference_genome_allele', 'mutated_from_allele', 'mutated_to_allele', 'subgroup', 'associated_region', 'region_metadata', 'min_pos_region', 'max_pos_region']]
  df2 = df.copy()
  df2['subgroup'] = 'all_subtypes'
  df = pd.concat([df, df2])
  df.to_csv(args.outputPath, sep="\t", index=False)
  return


def annotate_mutation_type(df):
  sbs_mask = (df['mutated_from_allele'].str.len() == df['mutated_to_allele'].str.len()) & (df['mutated_from_allele'].str.len() == 1)
  insertion_mask = (df['mutated_from_allele'].str.len() < df['mutated_to_allele'].str.len()) & (df['mutated_to_allele'].str.len() <= 200)
  deletion_mask = (df['mutated_from_allele'].str.len() > df['mutated_to_allele'].str.len()) & (df['mutated_from_allele'].str.len() <= 200)
  mbs_mask = (df['mutated_from_allele'].str.len() == df['mutated_to_allele'].str.len()) & (df['mutated_from_allele'].str.len() >= 2) & (df['mutated_from_allele'].str.len() <= 200)
  df.loc[sbs_mask, "mutation_type"] = "single base substitution"
  df.loc[insertion_mask, "mutation_type"] = "insertion of <=200bp"
  df.loc[deletion_mask, "mutation_type"] = "deletion of <=200bp"
  df.loc[mbs_mask, "mutation_type"] = "multiple base substitution (>=2bp and <=200bp)"
  return df


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--inputPath", type=str, required=True, help="Help goes here") 
  parser.add_argument("-r", "--regionsPath", type=str, required=True, help="Help goes here") 
  parser.add_argument("-c", "--chromosome", type=str, required=True, help="Help goes here") 
  parser.add_argument("-m", "--metadata", type=str, required=True, help="Help goes here") 
  parser.add_argument("-o", "--outputPath", type=str, required=True, help="Help goes here")
  args = parser.parse_args()
  
  main(args)
