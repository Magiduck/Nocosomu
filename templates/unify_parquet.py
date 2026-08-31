#!/usr/bin/env python3

import duckdb
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import duckdb
import numpy as np
import glob

all_parquet_files = glob.glob('./*.parquet')
all_parquet_files.sort()

parquet_schema = pq.read_schema(all_parquet_files[1], memory_map=True)

def write_parquet(df):
  table = pa.Table.from_pandas(df, schema=parquet_schema)
  del(df)
  parquet_writer2.write_table(table)
  del(table)

def run_all():
  for pq_file in all_parquet_files:
    print(f"{pq_file}")
    print(" Reading ")
    donor_mutation_file = duckdb.query(f"SELECT * from '{pq_file}'").df()
    print(" Writing ")
    write_parquet(donor_mutation_file)
  return(0)

parquet_writer2 = pq.ParquetWriter("unified.parquet", parquet_schema, compression='snappy')
run_all()
parquet_writer2.close()
