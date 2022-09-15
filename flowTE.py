from pathlib import Path
import utils.flow_te_help as helper
from utils.find_repeats import red_repeat_finder


# Genome/library in fasta format
input_fasta = Path(helper.args.fasta)
temp_dir = Path('tmp')
redout_dir = temp_dir/'redout'
renamed_fasta = Path(f'{input_fasta.stem}.fa')
# Output files
masked_fasta = f'{renamed_fasta.stem}.msk'
repeats = f'{renamed_fasta.stem}.rpt'

red_repeat_finder(input_fasta, temp_dir, redout_dir)