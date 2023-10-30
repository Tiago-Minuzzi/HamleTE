#!/usr/bin/env python3

import time
import tomli
import shutil
import pandas as pd
import utils.hamlete_help as helper
from pathlib import Path
from utils.find_repeats import red_repeat_finder
from utils.get_repeats import repeats_to_fasta
from utils.clustering import cluster_sequences
from utils.te_predict import Predictor, get_seq_from_pred, te_count
from utils.te_predict import get_TE_table, concat_pred_tables, concat_fastas
from utils.get_fasta import get_selected_sequences

# Get date and time
time_label = time.strftime('%y%m%d%H%M%S')
hamlete_start = time.perf_counter()

# Get HamleTE directory
hamlete_dir = Path(__file__).parent

# Model info TOML
# load TOML file
models_toml = open(hamlete_dir/'models/models_info.toml', 'rb')
models_info = tomli.load(models_toml)

# models info
model_01 = models_info['te_model']
model_02 = models_info['class_model']
model_03 = models_info['retro_model']
model_04 = models_info['dna_model']
model_05 = models_info['ltr_model']
model_06 = models_info['nonltr_model']

# k-mer length
kmer_length = helper.args.len_kmer

# minimum repeat sequence length
repeat_length = helper.args.min_len

# cutoff values
te_cutoff = helper.args.cutoff
sfam_cutoff = helper.args.label_cutoff

# batch size value
MAX_PRED_BATCH = 250
batch_value = helper.args.batch_value
batch_value = batch_value if batch_value <= MAX_PRED_BATCH else MAX_PRED_BATCH

# Other command-line arguments
progress_bar = helper.args.nobar  # disable progress bar
clustering = helper.args.noclust  # disable clustering
mode = helper.args.mode

# Genome/library in fasta format
input_fasta = Path(helper.args.fasta)
base_name = input_fasta.stem
renamed_fasta = Path(f'{base_name}.fa')

# Out directories
temp_dir = Path('tmp')
redout_dir = temp_dir/'redout'
output_directory = Path(helper.args.output_dir)

# Output files
# Output prefix
hamlete_prefix = f'HamleTE_{base_name}_{time_label}'

# Red masked genome
masked_fasta = f'{base_name}.msk'
masked_fasta_location = temp_dir / masked_fasta

# Red repeat coordinates
repeats = f'{base_name}.rpt'
repeats_location = temp_dir / repeats

# Repeats to fasta file
repeats_fasta = 'tmp_repeats.fasta'
repeats_fasta_location = temp_dir / repeats_fasta

# Clustered fasta
clustered_fasta = 'tmp_clustered.fasta'
clustered_fasta_location = temp_dir / clustered_fasta

# TE prediction dataframe
step01_te_pred_df = temp_dir / 'tmp_TE.tsv'
step01_te_fasta = temp_dir / 'tmp_TE.fasta'
tmp_te_table = temp_dir / 'tmp_TE_TABLE.tsv'

# Class prediction
step02_te_pred_df = temp_dir / 'tmp_CLASS.tsv'
step02_te_fasta = temp_dir / 'tmp_CLASS.fasta'
pred_retro_fasta = temp_dir / 'tmp_RETRO.fasta'
pred_dna_fasta = temp_dir / 'tmp_DNA.fasta'

# Retro LTR/non-LTR prediction
step03_te_pred_df = temp_dir / 'tmp_RETRO_SORD.tsv'
pred_retro_fasta = temp_dir / 'tmp_RETRO_SORD.fasta'
pred_ltr_fasta = temp_dir / 'tmp_RETRO_LTR.fasta'
pred_nonltr_fasta = temp_dir / 'tmp_RETRO_nonLTR.fasta'

# DNA TE prediction
step04_te_pred_df = temp_dir / 'tmp_DNA_FINAL.tsv'
dna_final_fasta = temp_dir / 'tmp_DNA_FINAL.fasta'

# Retro non-LTR prediction
step05_te_pred_df = temp_dir / 'tmp_LTR_FINAL.tsv'
ltr_final_fasta = temp_dir / 'tmp_LTR_FINAL.fasta'

# Retro non-LTR prediction
step06_te_pred_df = temp_dir / 'tmp_nonLTR_FINAL.tsv'
nonltr_final_fasta = temp_dir / 'tmp_nonLTR_FINAL.fasta'

# Final files
final_prediction_table = output_directory / f'{hamlete_prefix}_PRD.tsv'
final_prediction_fasta = output_directory / f'{hamlete_prefix}_SQS.fasta'
final_prediction_counts = output_directory / f'{hamlete_prefix}_CNT.tsv'

# -----------//-----------

# Check if input file exists
if not input_fasta.exists():
    print('>>> ERROR: File not found')
    exit(1)

if mode == 'a':
    # Find repeats using Red
    print("\n### Starting repeat detector ###")
    red_repeat_finder(input_fasta,
                      temp_dir,
                      redout_dir,
                      klen=kmer_length)

    # Get masked repeats
    repeats_to_fasta(masked_fasta_location,
                     repeats_location,
                     repeats_fasta_location,
                     repeat_length)

    # Cluster sequences using cd-hit-est
    if clustering:
        cluster_sequences(repeats_fasta_location,
                          clustered_fasta_location)
    else:
        clustered_fasta_location = repeats_fasta_location
elif mode == 'c':
    clustered_fasta_location = input_fasta
    temp_dir.mkdir(exist_ok=True)
# Predict TEs from clustered repeats
print(f'\n### Running model {model_01["name"]} ###')
pred_01 = Predictor(hamlete_dir / model_01['location'], model_01['labels'])
pred_01.label_prediction(clustered_fasta_location,
                         step01_te_pred_df,
                         batch_size_value=batch_value,
                         no_bar=progress_bar)

get_seq_from_pred(step01_te_pred_df,
                  'TE',
                  clustered_fasta_location,
                  step01_te_fasta,
                  cut_value=te_cutoff)

if step01_te_fasta.exists():
    # Predict TE class
    print(f'\n### Running model {model_02["name"]} ###')
    pred_02 = Predictor(hamlete_dir / model_02['location'], model_02['labels'])
    pred_02.label_prediction(step01_te_fasta,
                             step02_te_pred_df,
                             batch_size_value=batch_value,
                             no_bar=progress_bar)

    get_seq_from_pred(step02_te_pred_df,
                      'Retro',
                      step01_te_fasta,
                      pred_retro_fasta,
                      cut_value=te_cutoff)

    get_seq_from_pred(step02_te_pred_df,
                      'DNA',
                      step01_te_fasta,
                      pred_dna_fasta,
                      cut_value=te_cutoff)

    # Predict LTR/non-LTR
    print(f'\n### Running model {model_03["name"]} ###')
    pred_03 = Predictor(hamlete_dir / model_03['location'], model_03['labels'])
    pred_03.label_prediction(pred_retro_fasta,
                             step03_te_pred_df,
                             batch_size_value=batch_value,
                             no_bar=progress_bar)

    get_seq_from_pred(step03_te_pred_df,
                      'LTR',
                      pred_retro_fasta,
                      pred_ltr_fasta,
                      cut_value=te_cutoff)

    get_seq_from_pred(step03_te_pred_df,
                      'nonLTR',
                      pred_retro_fasta,
                      pred_nonltr_fasta,
                      cut_value=te_cutoff)

    # Predict DNA TE label
    print(f'\n### Running model {model_04["name"]} ###')
    pred_04 = Predictor(hamlete_dir / model_04['location'], model_04['labels'])
    pred_04.label_prediction(pred_dna_fasta,
                             step04_te_pred_df,
                             batch_size_value=batch_value,
                             no_bar=progress_bar)
    pred_04.filter(step04_te_pred_df, cut_value=sfam_cutoff)
    get_selected_sequences(pred_dna_fasta, step04_te_pred_df, dna_final_fasta)

    # Predict LTR label
    print(f'\n### Running model {model_05["name"]} ###')
    pred_05 = Predictor(hamlete_dir / model_05['location'], model_05['labels'])
    pred_05.label_prediction(pred_ltr_fasta,
                             step05_te_pred_df,
                             batch_size_value=batch_value,
                             no_bar=progress_bar)
    pred_05.filter(step05_te_pred_df, cut_value=sfam_cutoff)
    get_selected_sequences(pred_ltr_fasta, step05_te_pred_df, ltr_final_fasta)

    # Predict nonLTR label
    print(f'\n### Running model {model_06["name"]} ###')
    pred_06 = Predictor(hamlete_dir / model_06['location'], model_06['labels'])
    pred_06.label_prediction(pred_nonltr_fasta,
                             step06_te_pred_df,
                             batch_size_value=batch_value,
                             no_bar=progress_bar)
    pred_06.filter(step06_te_pred_df, cut_value=sfam_cutoff)
    get_selected_sequences(pred_nonltr_fasta,
                           step06_te_pred_df,
                           nonltr_final_fasta)

    if not output_directory.exists():
        print(f'\n### Creating directory {output_directory}')
        output_directory.mkdir()

    # Get TE and class predictions
    get_TE_table(step01_te_pred_df,
                 step02_te_pred_df,
                 step03_te_pred_df,
                 tmp_te_table,
                 cut_value=te_cutoff)

    # Concatenate final predictions
    concat_pred_tables(final_prediction_table,
                       tmp_te_table,
                       step05_te_pred_df,
                       step06_te_pred_df,
                       step04_te_pred_df,
                       mode)

    # Concatenate final fastas
    concat_fastas(final_prediction_fasta,
                  ltr_final_fasta,
                  nonltr_final_fasta,
                  dna_final_fasta)

    if final_prediction_table.exists() and final_prediction_fasta.exists():
        # Remove temporary directory
        shutil.rmtree(temp_dir)
        # Create prediction count table
        counts = te_count(final_prediction_table, mode)
        counts.to_csv(final_prediction_counts, index=False, sep='\t')

    hamlete_end = time.perf_counter()
    hamlete_total = hamlete_end - hamlete_start  # compute total run time
    print(f'\n>>> HamleTE {"classifier" if mode=="c" else "annotation"} mode finished in {hamlete_total:.2f} seconds.')
else:
    print('>>> No TEs found.')
models_toml.close()
