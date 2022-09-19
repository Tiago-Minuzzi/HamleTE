from pathlib import Path
import utils.flow_te_help as helper
from utils.find_repeats import red_repeat_finder
from utils.get_repeats import repeats_to_fasta
from utils.clustering import cluster_sequences
from utils.te_predict import label_pred_dataframe, label_prediction

# Genome/library in fasta format
input_fasta = Path(helper.args.fasta)
base_name = input_fasta.stem
renamed_fasta = Path(f'{base_name}.fa')

# Out directories
temp_dir = Path('tmp')
redout_dir = temp_dir/'redout'

# Output files
## Red masked genome
masked_fasta = f'{base_name}.msk'
masked_fasta_location = temp_dir / masked_fasta

## Red repeat coordinates
repeats = f'{base_name}.rpt'
repeats_location = temp_dir / repeats

## Repeats to fasta file
repeats_fasta = f'{base_name}_repeats.fasta'
repeats_fasta_location = temp_dir / repeats_fasta

## Clustered fasta
clustered_fasta = f'{base_name}_clustered.fasta'
clustered_fasta_location = temp_dir / clustered_fasta

## TE prediction dataframe
step01_te_pred_df = temp_dir / f'{base_name}_te_prediction.tsv'

if helper.args.mode == 'g':
# Find repeats using Red
    red_repeat_finder(input_fasta, temp_dir, redout_dir)

# Get masked repeats
    repeats_to_fasta(masked_fasta_location, repeats_location, repeats_fasta_location)

# Cluster sequences using cd-hit-est
    cluster_sequences(repeats_fasta_location, clustered_fasta_location)

if not clustered_fasta_location.exists():
    clustered_fasta_location = input_fasta
    temp_dir.mkdir(exist_ok=True)
# Predict TEs from clustered repeats
    print('### STEP 01: TE prediction ###')
    label_prediction(clustered_fasta_location, step01_te_pred_df)
    print('Done!')