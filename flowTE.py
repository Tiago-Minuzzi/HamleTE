import tomli
from pathlib import Path
import multiprocessing
import utils.flow_te_help as helper
from utils.find_repeats import red_repeat_finder
from utils.get_repeats import repeats_to_fasta
from utils.clustering import cluster_sequences
from utils.te_predict import Predictor, get_seq_from_pred

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
step01_te_pred_df = temp_dir / f'{base_name}_TE_prediction.tsv'
step01_te_fasta = temp_dir / f'{base_name}_TE.fasta'

## Class prediction
step02_te_pred_df = temp_dir / f'{base_name}_CLASS_prediction.tsv'
step02_te_fasta = temp_dir / f'{base_name}_CLASS.fasta'
pred_retro_fasta = temp_dir / f'{base_name}_RETRO.fasta'
pred_dna_fasta = temp_dir / f'{base_name}_DNA.fasta'


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
## Get sequences
get_seq_from_pred(step01_te_pred_df, 'TE', clustered_fasta_location, step01_te_fasta)

# Predict TE class
print(f'### Running model {model_02["name"]} ###')
pred_02 = Predictor(model_02['location'], model_02['labels'])
pred_02.label_prediction(step01_te_fasta, step02_te_pred_df)

retro_mp = multiprocessing.Process(target=get_seq_from_pred,args=[step02_te_pred_df, 'Retro', step01_te_fasta, pred_retro_fasta])
dna_mp = multiprocessing.Process(target=get_seq_from_pred,args=[step02_te_pred_df, 'DNA', step01_te_fasta, pred_dna_fasta])

retro_mp.start()
dna_mp.start()

print('### DONE! ###')

models_toml.close()