#!/usr/bin/env python3

import os
import sys
import pathlib
import numpy as np
import pandas as pd
from utils.prediction_utils import batch_iterator, label_pred_dataframe
from Bio.SeqIO.FastaIO import SimpleFastaParser
from math import ceil
from tqdm import tqdm
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

    def label_prediction(self, in_fasta: str, out_table: str, batch_size_value: int = 32, no_bar: bool = False) -> None:
        '''Run model to predict classes and return the predictions as a tsv file.'''
        modelo = self.location
        colunas = self.labels
        PADVALUE = 30_000
        modelo = load_model(modelo)
        predictions = []
        nt_to_token = { 'a':1, 't':2, 'g':3, 'c':4 }
        if in_fasta.exists():
            n_seqs = len([ i for i in open(in_fasta) if i.startswith('>') ])
            total_batches = ceil(n_seqs / batch_size_value)
            print(f'    - Predicting on {n_seqs} sequences divided in {total_batches} batch(es).')
            with open(in_fasta) as fa:
                for record in tqdm(batch_iterator(SimpleFastaParser(fa), batch_size_value), desc="    - Status", total=total_batches, unit='', ascii=' =',disable=no_bar):
                    identifiers = []
                    sequences = []
                    for fid, fsq in record:
                        identifiers.append(fid)
                        # Tokenize sequences
                        sequences.append([ nt_to_token[nt] if nt in nt_to_token.keys() else 5 for nt in fsq.lower() ])
                    # Pad sequences
                    padded_seqs = pad_sequences(sequences, padding='post',
                                                maxlen = PADVALUE, 
                                                truncating='post', 
                                                dtype='uint8')
                    # Predict labels
                    pred_values = modelo.predict(padded_seqs,
                                                 batch_size = batch_size_value,
                                                 verbose = 0)
                    identifiers = pd.Series(identifiers, name='id')
                    results_df = label_pred_dataframe(identifiers, pred_values, colunas)
                    predictions.append(results_df)
            predictions = pd.concat(predictions)
            predictions.to_csv(out_table, index=False, sep='\t')


    def filter(self, pred_table: str, cut_value: int = None) -> None:
        '''If cutoff value is set, predictions below it are classified as unknown.'''
        filter_table = pd.read_table(pred_table)
        if cut_value:
            filter_table['prediction'] = np.where(filter_table.select_dtypes('float').max(axis=1) < cut_value,
                                                 'Unknown',
                                                  filter_table.select_dtypes('float').idxmax(axis=1))
        filter_table.to_csv(pred_table,index=False,sep='\t')


def get_seq_from_pred(pred_table: str, label: str, reference_fasta: str, out_fasta: str, cut_value=None) -> None:
    '''Get sequences exact predicted label.'''
    if pred_table.exists():
        df = pd.read_table(pred_table)
        if cut_value:
            df = df.loc[df[label] >= cut_value]
        label_ids = df.loc[df['prediction']==label]['id'].to_list()
        reference_fasta = Path(reference_fasta)
        if label_ids:
            print(f'    - Retrieving {label} sequences...')
            with open(reference_fasta) as fa, open(out_fasta,'w') as sd:
                for fid, fsq in SimpleFastaParser(fa):
                        if fid in label_ids:
                            record = f'>{fid}|{label}\n{fsq}\n'
                            sd.write(record)


def prediction_processing(dataframe: pd.DataFrame) -> pd.DataFrame:
    '''Format final prediction dataframes.'''
    df = pd.read_table(dataframe)
    df['accuracy'] = df.select_dtypes('float').max(axis=1)
    df = df[['id','prediction','accuracy']]
    prev_pred = df['id'].str.split('|').to_list() # get previous prediction label
    df['id'] = df['id'].str.split('|',expand=True)[0]
    prev_pred = pd.Series([ i[-1] for i in prev_pred ])
    df['prediction'] = prev_pred + '|' + df['prediction'] # concatenate previous and last prediction
    return df


def te_count(dataframe: pd.DataFrame) -> pd.DataFrame:
    '''Return dataframe with counts for each label.'''
    df = pd.read_table(dataframe)
    counts = df['prediction'].value_counts().rename_axis('prediction').reset_index(name='count')
    return counts