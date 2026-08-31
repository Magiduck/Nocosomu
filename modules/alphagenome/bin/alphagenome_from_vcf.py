#!/usr/bin/env python3
import os
import sys
import json
import pandas as pd
from tqdm import tqdm

print("NEW SCRIPT")
if len(sys.argv) != 3:
    print("Usage: alphagenome_predict.py <input_file.tsv> <output_file.csv>")
    sys.exit(1)

input_file = sys.argv[1]
output_filename = sys.argv[2]
trace_log = "alphagenome_variant_trace.log"

os.environ['KAGGLE_CONFIG_DIR'] = '/home/umcg-tvanlieshout/.kaggle'
cache_dir = '/groups/umcg-fg/tmp04/projects/non-coding-somatic/kaggle_cache'
os.environ['KAGGLEHUB_CACHE'] = cache_dir
os.makedirs(cache_dir, exist_ok=True)

kaggle_json_path = '/home/umcg-tvanlieshout/.kaggle/kaggle.json'
try:
    with open(kaggle_json_path, "r") as f:
        kaggle_creds = json.load(f)
        os.environ['KAGGLE_USERNAME'] = kaggle_creds['username']
        os.environ['KAGGLE_KEY'] = kaggle_creds['key']
        print(f"Successfully loaded Kaggle credentials for user: {kaggle_creds.get('username')}")
except Exception as e:
    print(f"WARNING: Could not read {kaggle_json_path}. Error details: {e}")

os.environ['XLA_FLAGS'] = '--xla_gpu_deterministic_ops'
os.environ['XLA_PYTHON_CLIENT_PREALLOCATE'] = 'false'
os.environ['XLA_PYTHON_CLIENT_ALLOCATOR'] = 'platform'
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

print("Importing modules...")
import kagglehub
kagglehub.login = lambda *args, **kwargs: print("Bypassed interactive kagglehub.login() prompt.")

original_download = kagglehub.model_download
def offline_model_download(handle, *args, **kwargs):
    import glob
    handle_path = os.path.join(cache_dir, "models", *handle.split('/'))
    versions = [v for v in glob.glob(os.path.join(handle_path, "*")) if os.path.isdir(v)]
    if versions:
        latest = sorted(versions)[-1]
        print(f"Network bypass active. Loading {handle} directly from {latest}")
        return latest
    
    print(f"FATAL ERROR: Offline model weights missing for {handle}.")
    print("You must execute the following command on the head node before resuming:")
    print(f"KAGGLEHUB_CACHE={cache_dir} python3 -c \"import kagglehub; kagglehub.model_download('{handle}')\"")
    sys.exit(1)

kagglehub.model_download = offline_model_download

import jax
from alphagenome_research.model import dna_model
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
from alphagenome_research.io import fasta

local_fasta_path = os.path.join(cache_dir, 'GRCh38.p13.genome.fa')
if not os.path.exists(local_fasta_path):
    print(f"FATAL ERROR: Offline FASTA file missing at {local_fasta_path}.")
    print("You must execute the following commands on the head node before resuming:")
    print(f"wget https://storage.googleapis.com/alphagenome/reference/gencode/hg38/GRCh38.p13.genome.fa -O {local_fasta_path}")
    print(f"wget https://storage.googleapis.com/alphagenome/reference/gencode/hg38/GRCh38.p13.genome.fa.fai -O {local_fasta_path}.fai")
    sys.exit(1)

original_extractor_init = fasta.FastaExtractor.__init__
def offline_extractor_init(self, fasta_path, *args, **kwargs):
    if 'storage.googleapis.com' in str(fasta_path):
        print(f"Network bypass active. Loading FASTA directly from {local_fasta_path}")
        fasta_path = local_fasta_path
    original_extractor_init(self, fasta_path, *args, **kwargs)

fasta.FastaExtractor.__init__ = offline_extractor_init

original_read_feather = pd.read_feather
def offline_read_feather(path, *args, **kwargs):
    if isinstance(path, str) and 'storage.googleapis.com' in path:
        filename = path.split('/')[-1]
        local_file_path = os.path.join(cache_dir, filename)
        if not os.path.exists(local_file_path):
            print(f"FATAL ERROR: Offline dependency missing at {local_file_path}.")
            print("You must execute the following command on the head node before resuming:")
            print(f"wget {path} -O {local_file_path}")
            sys.exit(1)
        print(f"Network bypass active. Loading dependency directly from {local_file_path}")
        path = local_file_path
    return original_read_feather(path, *args, **kwargs)

pd.read_feather = offline_read_feather

print(f"Loading AlphaGenome model weights into {cache_dir}...")
model = dna_model.create_from_kaggle('all_folds')

vcf = pd.read_csv(input_file, sep='\t')
required_columns = ['variant_id', 'CHROM', 'POS', 'REF', 'ALT']
for column in required_columns:
    if column not in vcf.columns:
        raise ValueError(f'Input file is missing required column: {column}.')

organism = dna_client.Organism.HOMO_SAPIENS
sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_1MB']

score_rna_seq = True
scorer_selections = {'rna_seq': score_rna_seq}

all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
selected_scorers = [
    all_scorers[key] for key in all_scorers if scorer_selections.get(key.lower(), False)
]

unsupported_scorers = [
    scorer for scorer in selected_scorers
    if organism.value not in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
]
for unsupported_scorer in unsupported_scorers:
    selected_scorers.remove(unsupported_scorer)

subset_to_target_tissues = True
target_biosamples = [
    'colonic mucosa',
    'large intestine',
    'transverse colon',
    'left colon',
    'sigmoid colon',
    'mucosa of descending colon',
    'Caco-2',
    'HCT116',
    'HT-29',
]

chunk_size = 10
total_rows = len(vcf)

with open(trace_log, "w") as f:
    f.write("variant_id\tstatus\n")

for i in tqdm(range(0, total_rows, chunk_size), total=total_rows // chunk_size + (total_rows % chunk_size > 0)):
    chunk_vcf = vcf.iloc[i:i + chunk_size]
    results = []
    status_dict = {}
    
    for _, vcf_row in chunk_vcf.iterrows():
        var_id = vcf_row.variant_id
        chrom_str = str(vcf_row.CHROM)
        if not chrom_str.startswith("chr"):
            chrom_str = f"chr{chrom_str}"
        
        ag_id = f"{chrom_str}:{int(vcf_row.POS)}:{vcf_row.REF}>{vcf_row.ALT}"
        
        try:
            variant = genome.Variant(
                chromosome=chrom_str,
                position=int(vcf_row.POS),
                reference_bases=vcf_row.REF,
                alternate_bases=vcf_row.ALT,
                name=var_id,
            )
            interval = variant.reference_interval.resize(sequence_length)
            
            variant_scores = model.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=selected_scorers,
                organism=organism,
            )
            results.append(variant_scores)
            status_dict[var_id] = {"status": "SCORED", "ag_id": ag_id}
        except Exception as e:
            status_dict[var_id] = {"status": f"ERROR: {str(e).replace(chr(10), ' ')}", "ag_id": ag_id}

    df_scores = pd.DataFrame()
    pre_filter_variants = set()
    seen_biosamples = []
    id_col = 'variant_id'

    if results:
        df_scores = variant_scorers.tidy_scores(results)
        id_col = 'variant_id' if 'variant_id' in df_scores.columns else 'variant_name'
        
        if not df_scores.empty and id_col in df_scores.columns:
            df_scores[id_col] = df_scores[id_col].apply(lambda x: x.name if hasattr(x, 'name') else str(x))
            pre_filter_variants = set(df_scores[id_col].unique())
            if 'biosample_name' in df_scores.columns:
                seen_biosamples = df_scores['biosample_name'].dropna().unique().tolist()
                
        if subset_to_target_tissues and not df_scores.empty:
            df_scores = df_scores[df_scores['biosample_name'].isin(target_biosamples)]
    
    surviving_variants = set()
    if not df_scores.empty and id_col in df_scores.columns:
        surviving_variants = set(df_scores[id_col].unique())

    with open(trace_log, "a") as f:
        for _, vcf_row in chunk_vcf.iterrows():
            var_id = vcf_row.variant_id
            info = status_dict.get(var_id, {"status": "UNKNOWN", "ag_id": "UNKNOWN"})
            status = info["status"]
            ag_id = info["ag_id"]
            
            if status == "SCORED":
                if var_id in surviving_variants:
                    f.write(f"{var_id}\tSUCCESS\n")
                elif var_id in pre_filter_variants:
                    biosample_str = ", ".join(seen_biosamples[:5])
                    f.write(f"{var_id}\tDROPPED_BY_TISSUE_FILTER\tSeen tissues: {biosample_str}...\n")
                else:
                    f.write(f"{var_id}\tDROPPED_NO_GENES_NEARBY\n")
            else:
                f.write(f"{var_id}\t{status}\n")

    if not df_scores.empty and 'biosample_name' in df_scores.columns:
        group_cols = [c for c in [id_col, 'gene_id', 'gene_name', 'gene_type', 'gene_strand'] if c in df_scores.columns]
        agg_dict = {c: 'mean' for c in ['raw_score', 'score'] if c in df_scores.columns}
        other_cols = {c: 'first' for c in df_scores.columns if c not in group_cols and c not in agg_dict}
        
        df_scores = df_scores.groupby(group_cols, as_index=False, observed=True).agg({**agg_dict, **other_cols})
        df_scores['biosample_name'] = 'colorectal_average'

    if not df_scores.empty:
        df_scores.to_csv(output_filename, mode='a', header=not os.path.exists(output_filename), index=False)

if not os.path.exists(output_filename):
    pd.DataFrame(columns=['variant_id', 'raw_score', 'gene_name']).to_csv(output_filename, index=False)
