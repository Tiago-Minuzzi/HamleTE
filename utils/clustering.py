import subprocess
from typing import TextIO


def cluster_sequences(repeats_fasta: str, clustered_fasta: str) -> TextIO:
    """Run cd-hit-est to cluster most similar predicted TE sequences."""
    #cdhit = subprocess.run(['cd-hit-est',
    #                        '-i',   repeats_fasta,
    #                        '-o',   clustered_fasta,
    #                        '-G',   '0',    # use local sequence identity.
    #                        '-g',   '1',    # chooses best representative sequence.
    #                        '-aS',  '0.8',  # alignment coverage for the shorter sequence.
    #                        '-c',   '0.8',  # global sequence identity to comply with the 80-80-80 rule.
    #                        '-b',   '500',  # alignment’s bandwidth.
    #                        '-M',   '0',    # Memory limit. 0 means no limit.
    #                        '-T',   '0'])   # Number of threads to use. 0 means all threads.
    vsearch = subprocess.run(['vsearch',
                              '--cluster_fast',
                              repeats_fasta,
                              '--id',           '0.8',
                              '--centroids',    clustered_fasta,
                              ])
