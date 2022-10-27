import sys
import shutil
import subprocess
from pathlib import Path

    
def red_repeat_finder(input_fasta: str, temp_dir: str, redout_dir: str) -> None:
    """Run Red and find repeats"""
    # Change extensions to 'fa'
    renamed_fasta = Path(f'{input_fasta.stem}.fa')

    # Create temporary directory, if not exists.    
    if not temp_dir.exists():
        temp_dir.mkdir()

    # Create Red directory, if not exists.
    if not redout_dir.exists():
        redout_dir.mkdir()

    # Move input fasta to temp directory
    if input_fasta.exists():
        fasta_link = temp_dir/renamed_fasta
        fasta_link.symlink_to(input_fasta.absolute())

    # Run Red
    red_sftw = subprocess.run(['Red',
                             '-gnm',temp_dir,
                             '-msk',redout_dir,
                             '-rpt',redout_dir])

    # If Red succeeds, move output files out of Red directory 
    if red_sftw.returncode == 0:
        # Output files
        masked_fasta = f'{renamed_fasta.stem}.msk'
        repeats = f'{renamed_fasta.stem}.rpt'
        
        # move masked fasta to tmp folder
        masked_fasta_location = redout_dir / masked_fasta
        masked_new_location = temp_dir / masked_fasta
        shutil.move(masked_fasta_location,masked_new_location)

        # move repeat locations file to tmp folder
        repeats_location = redout_dir / repeats
        repeats_new_location = temp_dir / repeats
        shutil.move(repeats_location,repeats_new_location)
        
    elif red_sftw.returncode != 0:
        print(f'Red execution error: {red_sftw.returncode}')

    fasta_link.unlink()
    redout_dir.rmdir()