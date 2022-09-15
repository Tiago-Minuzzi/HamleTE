#!/usr/bin/env python3
import os
import sys
import pathlib
import numpy as np
import pandas as pd
from Bio import SeqIO
# Hide warnings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
stderr = sys.stderr
sys.stderr = open(os.devnull, 'w')
from keras.preprocessing.text import Tokenizer
sys.stderr = stderr


class InputFile:
    def __init__(self,in_file):
        # file path
        self.in_file = pathlib.Path(in_file)
        # file stem
        self.base_name = self.in_file.stem
        # file parent directory
        self.in_file_dir = self.in_file.absolute().parent
        # crete directory to store output files
        self.results_dir = self.in_file_dir / f'MITE_bdeep_results_{self.base_name}'
        # resulting file from predictions
        self.pred_results_tsv = self.base_name + '_candidates.tsv'


def fasta_frame(fasta_file):
    # Initialize fasta ids and fasta sequences lists
    fids = []
    fseq = []
    with open(fasta_file) as fasta:
        # Parse fasta file
        for record in SeqIO.parse(fasta, 'fasta'):
            fids.append(record.description) # append ids to list
            fseq.append(str(record.seq).lower()) # append sequences to list
    # lists to pandas series
    s1 = pd.Series(fids, name = 'id')
    s2 = pd.Series(fseq, name = 'sequence')
    # create dictionary
    data = {'id': s1, 'sequence': s2}
    # create dataframe
    df = pd.concat(data, axis=1)
    return df


def tokenize_sequences(sequencias):
    # Initialize tokenizer
    tkz_seq = Tokenizer(num_words = None, split = ' ', char_level = True, lower = True)
    # fit fasta sequences to text
    tkz_seq.fit_on_texts(sequencias)
    # tokenize sequences
    x_seq_arrays = tkz_seq.texts_to_sequences(sequencias)
    return x_seq_arrays
