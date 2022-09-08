import sys
from typing import TextIO
from Bio import SeqIO

in_fasta = sys.argv[1]
in_ids_file = sys.argv[2]
out_fasta = sys.argv[3]

def get_selected_sequences(in_fasta: str, in_ids_file: str, out_fasta: str) -> TextIO:
    with open(in_fasta) as fa, open(in_ids_file) as fids, open(out_fasta,'w') as sd:
        fids = [ i.strip() for i in fids.readlines() ]
        for record in SeqIO.parse(fa, 'fasta'):
            if record.description in fids:
                SeqIO.write(record,sd,'fasta')


if __name__ == "__main__":
    get_selected_sequences(in_fasta, in_ids_file, out_fasta)