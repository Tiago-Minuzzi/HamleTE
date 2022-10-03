#!/usr/bin/env python3

import os
import sys
import pathlib
import numpy as np
import pandas as pd
from utils.prediction_utils import batch_iterator, label_pred_dataframe
from Bio.SeqIO.FastaIO import SimpleFastaParser
from numpy import array
from numpy import argmax
from pathlib import Path
from dataclasses import dataclass
# Hide warning messages
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
stderr = sys.stderr
sys.stderr = open(os.devnull, 'w')
import tensorflow as tf
from tensorflow.keras.models import load_model
from keras.preprocessing.sequence import pad_sequences
sys.stderr = stderr

@dataclass
class Predictor:
    location: str
    labels: list

    def label_prediction(self, in_fasta: str, out_table: str, batch_size_value: int = 500):
        '''Run model to predict classes and return the predictions as a tsv file'''
        modelo = self.location
        colunas = self.labels
        PADVALUE = 30_000
        MAX_PRED_BATCH = 250
        # Load model
        modelo = load_model(modelo)
        predictions = []
        nt_to_token = { 'a':1, 't':2, 'g':3, 'c':4 }
        if in_fasta.exists():
            with open(in_fasta) as fa:
                for record in batch_iterator(SimpleFastaParser(fa), batch_size_value):
                    identifiers = []
                    sequences = []
                    for fid, fsq in record:
                        identifiers.append(fid)
                        # Tokenize sequences
                        sequences.append([ nt_to_token[nt] if nt in nt_to_token.keys() else 5 for nt in fsq.lower() ])
                    # Pad sequences
                    padded_seqs = pad_sequences(sequences, padding='post', maxlen = PADVALUE)
                    pred_values = modelo.predict(padded_seqs,
                                                batch_size = batch_size_value if batch_size_value <= MAX_PRED_BATCH else MAX_PRED_BATCH,
                                                verbose = 1)
                    # Predict labels
                    identifiers = pd.Series(identifiers, name='id')
                    results_df = label_pred_dataframe(identifiers, pred_values, colunas)
                    predictions.append(results_df)
            predictions = pd.concat(predictions)
            predictions.to_csv(out_table, index=False, sep='\t')


def get_seq_from_pred(pred_table: str, label: str, reference_fasta: str, out_fasta: str) -> None:
    df = pd.read_table(pred_table)
    label_ids = df.loc[df['prediction']==label]['id'].to_list()
    reference_fasta = Path(reference_fasta)
    if label_ids:
        print(f'    Retrieving {label} sequences...')
        with open(reference_fasta) as fa, open(out_fasta,'w') as sd:
            for fid, fsq in SimpleFastaParser(fa):
                    if fid in label_ids:
                        record = f'>{fid}|{label}\n{fsq}\n'
                        sd.write(record)