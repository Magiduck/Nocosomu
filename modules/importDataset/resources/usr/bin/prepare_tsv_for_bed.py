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
singularity \
exec \
-B \
/groups/umcg-fg/tmp01/projects/non-coding-somatic/ \
/groups/umcg-fg/tmp01/projects/non-coding-somatic/shared/singularity/containers_cgut/sequence_based_models.sif \
python3 /groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-tvanlieshout/repos_tvanlieshout/nocosomu_iDriver/modules/importDataset/bin/prepare_tsv_for_bed.py \
-i /groups/umcg-fg/tmp01/projects/non-coding-somatic/cancer/genomics_england/hmf_regions_frequency_matrices/gel_all_pseudodonor.tsv.gz \
-c 5 \
-o gel_all_pseudodonor.bed
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
  df = pd.read_csv(args.inputPath, sep="\t")
  df['#CHROM'] = df['#CHROM'].str.replace("chr", "")
  df = df[df['#CHROM'] == args.chromosome]
  df['#CHROM'] = "chr" + df['#CHROM']
  df['START'] = df['POS'] - 1
  df['STOP'] = df['POS']
  df['NAME'] = df['#CHROM'] + "_" + df["POS"].astype(str) + "_" + df["REF"] + "_" + df["ALT"]
  df['SCORE'] = "."
  df['STRAND'] = "."
  df = df[['#CHROM', 'START', 'STOP', 'NAME', 'SCORE', 'STRAND']]
  print(df)
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
  parser.add_argument("-c", "--chromosome", type=str, required=True, help="Help goes here") 
  parser.add_argument("-o", "--outputPath", type=str, required=True, help="Help goes here")
  args = parser.parse_args()
  
  main(args)
