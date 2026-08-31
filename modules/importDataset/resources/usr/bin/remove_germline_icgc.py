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
--nv \
-B \
/groups/umcg-fg/tmp02/projects/non-coding-somatic/ \
/groups/umcg-fg/tmp02/projects/non-coding-somatic/shared/singularity/containers_cgut/sequence_based_models.sif \
python3 /groups/umcg-fg/tmp02/projects/non-coding-somatic/umcg-tvanlieshout/repos_tvanlieshout/nocosomu_iDriver/modules/importDataset/bin/remove_germline_icgc.py \
-i /groups/umcg-fg/tmp02/projects/non-coding-somatic/umcg-tvanlieshout/scratch_tvanlieshout/nocosomu_iDriver/work/1e/ca1f3fe1a3375d389d74f9d6f6f329/merged_mutations_in_region.tsv \
-g /groups/umcg-fg/tmp02/projects/non-coding-somatic/shared/gnomad/artifacts/hg19/chrY_TrueSNV_TruePASS.csv \
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
  df = pd.read_csv(args.inputPath, sep="\t")
  df["var_id"] = df["chromosome"].astype(str) + "_" + df["chromosome_end"].astype(str) + "_" + df['mutated_from_allele'] + "_" + df['mutated_to_allele']
  if os.path.exists(args.germlinePath):
    germline = pd.read_csv(args.germlinePath, sep="\t", header=None, names=["CHROM", "POS", "REF", "ALT", "FILTER", "FAF"])
    germline = germline[germline["FAF"] >= args.alleleFrequency]
    forbidden = germline["CHROM"].str[3:] + "_" + germline["POS"].astype(str) + "_" + germline['REF'] + "_" + germline['ALT']
    df = df[~df["var_id"].isin(forbidden)]
  df = df.drop(columns=['var_id'])
  df.to_csv(args.outputPath, sep="\t", index=False)
  return


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--inputPath", type=str, required=True, help="Help goes here") 
  parser.add_argument("-g", "--germlinePath", type=str, required=True, help="Help goes here") 
  parser.add_argument("-o", "--outputPath", type=str, required=True, help="Help goes here")
  parser.add_argument("-a", "--alleleFrequency", type=float, required=True, help="Help goes here")
  args = parser.parse_args()
  
  main(args)
