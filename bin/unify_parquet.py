#!/usr/bin/env python3
import argparse
import glob

import duckdb
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import duckdb
import numpy as np
import pysam


def main(args):
  all_parquet_files = glob.glob('./*.parquet')
  all_parquet_files.sort()

  parquet_schema = pq.read_schema(all_parquet_files[1], memory_map=True)
  parquet_schema = parquet_schema.append(pa.field('tri_context_ref', pa.string())).append(pa.field('tri_context_alt', pa.string()))
  parquet_writer2 = pq.ParquetWriter("unified.parquet", parquet_schema, compression='snappy')
  run_all(all_parquet_files, parquet_writer2, parquet_schema, args.genomePath)
  parquet_writer2.close()


def run_all(all_parquet_files, parquet_writer2, parquet_schema, genome_file):
  for pq_file in all_parquet_files:
    print(f"{pq_file}")
    print(" Reading ")
    donor_mutation_file = duckdb.query(f"SELECT * from '{pq_file}'").df()
    print("adding trinucleotide context")
    donor_mutation_file = add_tri_nucleotide_context(genome_file, donor_mutation_file)
    print(" Writing ")
    write_parquet(donor_mutation_file, parquet_writer2, parquet_schema)
  return(0)


def add_tri_nucleotide_context(genome_file, mutations):
  ref = pysam.FastaFile(genome_file)
  mutations['tri_context_ref'] = ""
  mutations['tri_context_alt'] = ""

  if len(mutations) == 0:
    return mutations
 
# commented out code is for a test that we do indeed grab the correct position
#  correct = 0
  for it, row in mutations.iterrows():
    tri_ref = ref.fetch(row['chromosome'], row['chromosome_start'] - 1, row['chromosome_start'] + 2)
    if len(tri_ref) > 0:
#      if tri_ref[1] == row['mutated_from_allele']:
#        correct += 1
      mutations.loc[it, 'tri_context_ref'] = tri_ref[0] + row['mutated_from_allele'] + tri_ref[2]
      mutations.loc[it, 'tri_context_alt'] = tri_ref[0] + row['mutated_to_allele'] + tri_ref[2]
  ref.close()
#  print(f"Correct reference in {100 / len(mutations) * correct}% of the commands")
  return mutations


def write_parquet(df, parquet_writer2, parquet_schema):
  table = pa.Table.from_pandas(df, schema=parquet_schema)
  del(df)
  parquet_writer2.write_table(table)
  del(table)


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-g", "--genomePath", type=str, required=True, help="Help goes here") 
  args = parser.parse_args()
  main(args)


