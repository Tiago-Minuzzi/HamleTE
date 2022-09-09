import sys
import shutil
import subprocess
from pathlib import Path

# Input file
input_fasta = Path(sys.argv[1])
renamed_fasta = Path(f'{input_fasta.stem}.fa')
# Temporary directory and Red output directory
temp_dir = Path('tmp')
redout_dir = temp_dir/'redout'
# Output files
masked_fasta = f'{renamed_fasta.stem}.msk'
repeats = f'{renamed_fasta.stem}.rpt'

def location_handler():
    # Create temporary directory, if not exists.
    if not temp_dir.exists():
        temp_dir.mkdir()

    # Create Red directory, if not exists.
    if not redout_dir.exists():
        redout_dir.mkdir()

    # Move input fasta to temp directory
    if input_fasta.exists():
        shutil.move(input_fasta, temp_dir/renamed_fasta)
    
def red_repeat_finder():
    """Run Red and find repeats"""
    location_handler()
    # Run Red
    red_sftw = subprocess.run(['Red',
                             '-gnm',temp_dir,
                             '-msk',redout_dir,
                             '-rpt',redout_dir])
                             
    # If Red succeeds, move output files out of Red directory 
    if red_sftw.returncode == 0:
        # move input fasta to original location
        shutil.move(temp_dir/renamed_fasta,input_fasta)
        
        # move masked fasta to tmp folder
        masked_fasta_location = redout_dir / masked_fasta
        masked_new_location = temp_dir / masked_fasta
        shutil.move(masked_fasta_location,masked_new_location)

        # move repeat locations file to tmp folder
        repeats_location = redout_dir / repeats
        repeats_new_location = temp_dir / repeats
        shutil.move(repeats_location,repeats_new_location)
        
        # remove redout directory
        redout_dir.rmdir()


if __name__ == "__main__":
    red_repeat_finder()