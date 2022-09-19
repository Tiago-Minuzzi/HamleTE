#!/usr/bin/env python3
import os
import sys
import pathlib
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
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


def fasta_reader(fasta_file: str):
    fids = []
    fsqs = []
    with open(fasta_file) as fa:
        for fid, fsq in SimpleFastaParser(fa):
            fids.append(fid)
            fsqs.append(fsq)
    return fids, fsqs


def tokenize_sequences(sequencias):
    # Initialize tokenizer
    tkz_seq = Tokenizer(num_words = None, split = ' ', char_level = True, lower = True)
    # fit fasta sequences to text
    tkz_seq.fit_on_texts(sequencias)
    # tokenize sequences
    x_seq_arrays = tkz_seq.texts_to_sequences(sequencias)
    return x_seq_arrays


