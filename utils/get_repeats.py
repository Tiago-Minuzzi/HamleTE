from pyfaidx import Fasta
from typing import TextIO


def repeats_to_fasta(masked_fasta: str, repeats_file: str, repeats_fasta: str) -> TextIO:
    """ Read repeat coordinates from Red output and save the repeats in a fasta file."""
    masked_fasta = Fasta(masked_fasta, read_long_names=True)
    with open(repeats_file) as rpts, open(repeats_fasta,'w') as sd:
        for linha in rpts:
            linha=linha.strip()[1:] # remove '>' from line
            # get sequence name and repeat coordinates
            fid, coords = linha.rsplit(':',1)
            sstart, send = coords.split('-')
            sstart, send = int(sstart), int(send)
            # get sequence
            if 50 <= (send - sstart) <= 30_000: # filter sequences by length
                fsq = masked_fasta[fid][sstart:send].seq
                record = f'>{linha}\n{fsq}\n'
                # write sequences to file
                sd.write(record)
                
