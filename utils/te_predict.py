#!/usr/bin/env python3

import os
import sys
import pathlib
import numpy as np
import pandas as pd
from prediction_utils import fasta_frame, tokenize_sequences
from Bio import SeqIO
from numpy import array
from numpy import argmax
# Hide warning messages
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
stderr = sys.stderr
sys.stderr = open(os.devnull, 'w')
import tensorflow as tf
from keras.preprocessing.text import Tokenizer
from tensorflow.keras.models import load_model
from keras.preprocessing.sequence import pad_sequences
sys.stderr = stderr


def label_pred_dataframe(fasta_ids, prediction_results):
    # Labels
    prediction_results = prediction_results[:,1:]
    # Create dataframe
    colunas = ['coding', 'noncod','te']
    label_pred_df = pd.DataFrame(prediction_results, columns = colunas)*100
    label_pred_df = label_pred_df.round(3)
    label_pred_df = pd.concat([fasta_ids,label_pred_df],axis=1)
    # Create label column
    label_pred_df['prediction'] = label_pred_df[colunas].idxmax(axis=1)
    return label_pred_df

def label_prediction(in_fasta, batch_size_value=4):
    modelo = '../models/te_identifier.hdf5'
    label_model = modelo
    PADVALUE = 30_000
    # Read fasta as dataframe
    fas_df = fasta_frame(in_fasta)
    identifiers = fas_df['id']
    sequences = fas_df['sequence']
    # Tokenize sequences
    tokenized_seqs = tokenize_sequences(sequences)
    # Pad sequences
    padded_seqs = pad_sequences(tokenized_seqs, padding='post', maxlen = PADVALUE)
    # Load model
    modelo = load_model(label_model)
    pred_values = modelo.predict(padded_seqs, batch_size = batch_size_value, verbose = 1)
    # Predict labels
    results_df = label_pred_dataframe(identifiers, pred_values)
    return results_df


if __name__ == '__main__':
    entrada = sys.argv[1]
    print('Running...')
    df = label_prediction(entrada)
    print(df)
    saida = sys.argv[2]
    df.to_csv(saida, index=False, sep='\t')
    print('Done!')

