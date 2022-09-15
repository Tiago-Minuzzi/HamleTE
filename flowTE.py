from pathlib import Path
import utils.flow_te_help as helper
from utils.find_repeats import red_repeat_finder
from utils.get_repeats import repeats_to_fasta

# Genome/library in fasta format
input_fasta = Path(helper.args.fasta)
temp_dir = Path('tmp')
redout_dir = temp_dir/'redout'
renamed_fasta = Path(f'{input_fasta.stem}.fa')

# Output files
## Red output
masked_fasta = f'{renamed_fasta.stem}.msk'
masked_fasta_location = temp_dir / masked_fasta

repeats = f'{renamed_fasta.stem}.rpt'
repeats_location = temp_dir / repeats

## fasta from Red repeat positions
repeats_fasta = f'{renamed_fasta.stem}_repeats.fasta'
repeats_fasta_location = temp_dir / repeats_fasta

# Find repeats using Red
red_repeat_finder(input_fasta, temp_dir, redout_dir)

# Get masked repeats
repeats_to_fasta(masked_fasta_location, repeats_location, repeats_fasta_location)