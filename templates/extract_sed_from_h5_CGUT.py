#!/usr/bin/env python

# Imports
import argparse

import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm

/groups/umcg-fg/tmp02/projects/non-coding-somatic/models/borzoi_stable//extract_sed_from_h5_CGUT.py         -i ./sed.h5         -t targets_gtex.txt         -x RNA:blood         -o ENSG00000164362_ism_f0_predictions.tsv

def main(args):
  hf = h5py.File(args.inputPath, 'r')

  hf_si = pd.Series(hf['si'])
  hf_genes = pd.Series(hf['gene']).str.decode('utf-8').str.split(".").str[0]
  
  targets_df = pd.read_csv(args.targetsPath, index_col=0, sep='\t')
  targets_df['local_index'] = np.arange(len(targets_df))

  if args.target == 'RNA:all':
    index_of_interest = targets_df.loc[:]['local_index'].tolist()
  else:
    index_of_interest = targets_df.loc[targets_df['description'] == args.target]['local_index'].tolist()

  hf_logSED = np.array(hf['logSED'])
  hf_logSED = pd.Series(np.mean(hf_logSED[:, index_of_interest], axis=1))
  
  hf_REF = np.array(hf['REF'])
  hf_REF = pd.Series(np.mean(hf_REF[:, index_of_interest], axis=1))
  
  hf_ALT = np.array(hf['ALT'])
  hf_ALT = pd.Series(np.mean(hf_ALT[:, index_of_interest], axis=1))
  
  hf_snps = pd.Series(hf['snp']).str.decode('utf-8')
  
  snps_df = pd.DataFrame(data={'si': np.arange(len(hf_snps)), 'ID': hf_snps})
  sed_df = pd.DataFrame(data={'si': hf_si, 'gene': hf_genes, 'hf_logSED': hf_logSED, 'hf_REF': hf_REF, 'hf_ALT': hf_ALT})
  
  all_sed_snps = pd.merge(snps_df, sed_df, on='si', how="inner")
  
  #annotations = pd.read_csv('qa_sanity_test_unsorted_GRCh37_annotation.tsv', sep="\t")
  #queried_sed_snps = pd.merge(all_sed_snps, annotations, on=['ID', 'gene'], how="inner")
  
  all_sed_snps.to_csv(args.outputPath, sep="\t", index=False) 

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--inputPath", type=str, required=True, help="Help goes here") 
  #parser.add_argument("-a", "--annotationsPath", type=str, required=True, help="Help goes here") 
  #parser.add_argument("-f", "--fileToPredictPath", type=str, required=True, help="Help goes here") 
  parser.add_argument("-t", "--targetsPath", type=str, required=True, help="Help goes here") 
  parser.add_argument("-x", "--target", type=str, required=True, help="Help goes here") 
  parser.add_argument("-o", "--outputPath", type=str, required=True, help="Help goes here")
  args = parser.parse_args()
  
  main(args)
