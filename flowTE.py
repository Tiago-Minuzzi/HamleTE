import tomli
from pathlib import Path
import utils.flow_te_help as helper
from utils.find_repeats import red_repeat_finder
from utils.get_repeats import repeats_to_fasta
from utils.clustering import cluster_sequences
from utils.te_predict import Predictor

# Model info TOML
## load TOML file
models_toml = open('models/models_info.toml','rb')
models_info = tomli.load(models_toml)

## models info
model_01 = models_info['te_model']
model_02 = models_info['class_model']
model_03 = models_info['retro_model']
model_04 = models_info['dna_model']
model_05 = models_info['ltr_model']
model_06 = models_info['nonltr_model']

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

elif helper.args.mode == 'c':
    clustered_fasta_location = input_fasta
    temp_dir.mkdir(exist_ok=True)
    
# Predict TEs from clustered repeats
print(f'### Running model {model_01["name"]} ###')

pred_01 = Predictor(model_01['location'], model_01['labels'])
pred_01.label_prediction(clustered_fasta_location, step01_te_pred_df)

print('Done!')

models_toml.close()