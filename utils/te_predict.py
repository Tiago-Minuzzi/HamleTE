#!/usr/bin/env python3

import os
import sys
import pathlib
import numpy as np
import pandas as pd
from utils.prediction_utils import tokenize_sequences, batch_iterator, label_pred_dataframe
from Bio.SeqIO.FastaIO import SimpleFastaParser
from numpy import array
from numpy import argmax
from dataclasses import dataclass
# Hide warning messages
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
stderr = sys.stderr
sys.stderr = open(os.devnull, 'w')
import tensorflow as tf
from keras.preprocessing.text import Tokenizer
from tensorflow.keras.models import load_model
from keras.preprocessing.sequence import pad_sequences
sys.stderr = stderr

@dataclass
class Predictor:
    location: str
    labels: list

    def label_prediction(self, in_fasta: str, out_table: str, batch_size_value: int = 4):
        '''Run model to predict classes and return the predictions as a tsv file'''
        modelo = self.location
        colunas = self.labels
        PADVALUE = 30_000
        # Load model
        modelo = load_model(modelo)
        predictions = []
        with open(in_fasta) as fa:
            for record in batch_iterator(SimpleFastaParser(fa), 50):
                identifiers = []
                sequences = []
                for fid, fsq in record:
                    identifiers.append(fid)
                    sequences.append(fsq)
                # Tokenize sequences
                tokenized_seqs = tokenize_sequences(sequences)
                # Pad sequences
                padded_seqs = pad_sequences(tokenized_seqs, padding='post', maxlen = PADVALUE)
                pred_values = modelo.predict(padded_seqs,
                                            batch_size = batch_size_value,
                                            verbose = 1)
                # Predict labels
                identifiers = pd.Series(identifiers)
                results_df = label_pred_dataframe(identifiers, pred_values, colunas)
                predictions.append(results_df)
        predictions = pd.concat(predictions)
        predictions.to_csv(out_table, index=False, sep='\t')

