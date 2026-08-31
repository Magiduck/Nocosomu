#!/usr/bin/env python3

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import duckdb
import glob

test_csv_file = "/groups/umcg-fg/tmp01/projects/non-coding-somatic/shared/cadd/whole_genome_SNVs_peek.txt"
parquet_file = "/groups/umcg-fg/tmp01/projects/non-coding-somatic/shared/cadd/whole_genome_SNVs.parquet"
chunksize = 100_000


columns_to_use = [
  "#Chrom",
  "Pos",
  "Ref",
  "Alt",
  "RawScore",
  "PHRED"
]

test_pd = pd.read_csv(test_csv_file, usecols= columns_to_use, skiprows=1, sep="\t")
# Guess the schema of the CSV file from the schema file
parquet_schema = pa.Table.from_pandas(df=test_pd).schema
# Adapt the fields that are not strings to their expected types
parquet_schema = parquet_schema.set(parquet_schema.get_field_index("#Chrom"), pa.field("#Chrom","string"))
parquet_schema = parquet_schema.set(parquet_schema.get_field_index("Pos"), pa.field("Pos","uint32"))

# Open a Parquet file for writing
parquet_writer = pq.ParquetWriter(parquet_file, parquet_schema, compression='zstd')

all_ssm_files = glob.glob('./whole_genome_SNVs.tsv.gz')
all_ssm_files.sort()

for ssm_file in all_ssm_files:
  print(ssm_file)
  csv_stream = pd.read_csv(ssm_file,
    compression='gzip', 
    sep='\t', 
    chunksize=chunksize, 
    low_memory=False, 
    usecols=columns_to_use,
    skiprows=1,
    dtype={
      '#Chrom': object,
      'Pos': int,
      'Ref': object,
      'Alt': object,
      'RawScore': float,
      'PHRED': float
    })
  for i, chunk in enumerate(csv_stream):
    print("Chunk", i)
    # Write CSV chunk to the parquet file
    table = pa.Table.from_pandas(chunk, schema=parquet_schema)
    parquet_writer.write_table(table)


parquet_writer.close()
