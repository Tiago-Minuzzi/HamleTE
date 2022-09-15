import sys
from typing import TextIO
from Bio import SeqIO


def get_selected_sequences(in_fasta: str, in_ids_file: str, out_fasta: str) -> TextIO:
    with open(in_fasta) as fa, open(in_ids_file) as fids, open(out_fasta,'w') as sd:
        fids = [ i.strip() for i in fids.readlines() ]
        for record in SeqIO.parse(fa, 'fasta'):
            if record.description in fids:
                SeqIO.write(record,sd,'fasta')
