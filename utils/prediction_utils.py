#!/usr/bin/env python3
import os
import sys
import pathlib
import numpy as np
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
# Hide warnings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
stderr = sys.stderr
sys.stderr = open(os.devnull, 'w')
from keras.preprocessing.text import Tokenizer
sys.stderr = stderr


def fasta_reader(fasta_file: str):
    '''Read fasta file.'''
    fids = [] # fasta identifiers list
    fsqs = [] # fasta sequences list
    with open(fasta_file) as fa:
        for fid, fsq in SimpleFastaParser(fa):
            fids.append(fid)
            fsqs.append(fsq)
    return fids, fsqs


def batch_iterator(iterator, size):
    '''Split large files in batches.'''
    batch = []
    for entry in iterator:
        batch.append(entry)
        if len(batch) == size:
            yield batch
            batch = []
    if batch:
        yield batch


def replace_nnt(sequence: str) -> str:
    '''Replace invalid nucleotides on sequence.'''
    not_nt = 'bdefhijklmopqrsuvwxyz'
    sequence = sequence.lower()
    for nnt in not_nt:
        if nnt in sequence:
            sequence = sequence.replace(nnt,'n')
    return sequence
    

def tokenize_sequences(sequencias):
    '''Transform nucleotides in tokens in a vector.'''
    # Initialize tokenizer
    tkz_seq = Tokenizer(num_words = None, split = ' ', char_level = True, lower = True)
    # fit fasta sequences to text
    tkz_seq.fit_on_texts(sequencias)
    # tokenize sequences
    x_seq_arrays = tkz_seq.texts_to_sequences(sequencias)
    return x_seq_arrays


def label_pred_dataframe(fasta_ids, prediction_results, colunas) -> pd.DataFrame:
    '''Return prediction values as dataframe.'''
    # Labels
    prediction_results = prediction_results[:,1:]
    # Create dataframe
    label_pred_df = pd.DataFrame(prediction_results, columns = colunas)
    label_pred_df = label_pred_df.round(3)
    label_pred_df = pd.concat([fasta_ids,label_pred_df],axis=1)
    # Create label column
    label_pred_df['prediction'] = label_pred_df[colunas].idxmax(axis=1)
    return label_pred_df

