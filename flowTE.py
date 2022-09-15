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
masked_fasta = f'{renamed_fasta.stem}.msk'
repeats = f'{renamed_fasta.stem}.rpt'
repeats_fasta = f'{renamed_fasta.stem}_repeats.fasta'


red_repeat_finder(input_fasta, temp_dir, redout_dir)

# repeats_to_fasta(masked_fasta, repeats, repeats_fasta)