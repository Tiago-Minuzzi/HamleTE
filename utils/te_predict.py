#!/usr/bin/env python3

import os
import sys
import shutil
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
pd.options.mode.chained_assignment = None  # default='warn'


@dataclass
class Predictor:
    location: str
    labels: list

    def label_prediction(self, in_fasta: str, out_table: str, batch_size_value: int = 32, no_bar: bool = False) -> None:
        '''Run model to predict classes and return the predictions as a tsv file.'''
        modelo      = self.location
        colunas     = self.labels
        PADVALUE    = 25_000
        try:
            modelo      = load_model(modelo)
            predictions = []
            nt_to_token = { 'a' : 1, 'c' : 2, 'g' : 3, 't' : 4 }
            if in_fasta.exists():
                n_seqs          = len([ i for i in open(in_fasta) if i.startswith('>') ])
                total_batches   = ceil(n_seqs / batch_size_value)
                print(f'    - Predicting on {n_seqs} sequences divided in {total_batches} batch(es).')
                with open(in_fasta) as fa:
                    for record in tqdm(batch_iterator(SimpleFastaParser(fa), batch_size_value), desc="    - Status", total=total_batches, unit='', ascii=' =',disable=no_bar):
                        identifiers = []
                        sequences   = []
                        for fid, fsq in record:
                            identifiers.append(fid)
                            # Tokenize sequences
                            sequences.append([ nt_to_token.get(nt,5) for nt in fsq.lower() ])
                        # Pad sequences
                        padded_seqs = pad_sequences(sequences,
                                                    maxlen      = PADVALUE,
                                                    padding     = 'post',
                                                    truncating  = 'post',
                                                    dtype       = 'uint8')
                        # Predict labels
                        pred_values = modelo.predict(padded_seqs,
                                                    batch_size  = batch_size_value,
                                                    verbose     = 0)
                        identifiers = pd.Series(identifiers, name='id')
                        results_df  = label_pred_dataframe(identifiers, pred_values, colunas)
                        predictions.append(results_df)
                predictions = pd.concat(predictions)
                predictions.to_csv(out_table, index=False, sep='\t')
        except OSError as eos:
            print(eos)
            print("  - Please, make sure the you've extracted the models in the 'models' directory.")
            print("  - Read the 'README.md' file for more information.",end='\n\n')
            exit(1)

    def filter(self, pred_table: str, cut_value: int = None) -> None:
        '''If cutoff value is set, predictions below it are classified as unknown.'''
        filter_table = pd.read_table(pred_table)
        if cut_value:
            filter_table['prediction'] = np.where(filter_table.select_dtypes('float').max(axis=1) < cut_value,
                                                  'Unknown',
                                                  filter_table.select_dtypes('float').idxmax(axis=1))
        filter_table.to_csv(pred_table, index=False, sep='\t')


def table_filter(df: pd.DataFrame) -> pd.DataFrame:
    """Filter columns of full prediction table resulting on id,
    prediction and accuracy columns only."""

    df['accuracy']  = df.select_dtypes('float').max(axis=1)
    df              = df[['id', 'prediction', 'accuracy']]
    df['id']        = df['id'].str.split('|').str[0]
    return df


def get_TE_table(te_pred_table, class_pred_table, ret_table, out_table, cut_value=None) -> None:
    """Get TE predicted sequences and class predictions
    and merge to further create final table."""

    df          = table_filter(pd.read_table(te_pred_table))
    df          = df.loc[df['prediction'] == 'TE']
    te_class    = table_filter(pd.read_table(class_pred_table))
    ret         = table_filter(pd.read_table(ret_table))
    merged      = df.merge(te_class, on='id', suffixes=['_1', '_2']).merge(ret, on='id', how='outer').fillna(te_class)

    if cut_value:
        merged  = merged.loc[merged['accuracy'] >= cut_value]

    merged.to_csv(out_table, sep='\t', index=False)


def get_seq_from_pred(pred_table: str, label: str, reference_fasta: str, out_fasta: str, cut_value=None) -> None:
    '''Get sequences exact predicted label.'''
    if pred_table.exists():
        df = pd.read_table(pred_table)
        if cut_value:
            df = df.loc[df[label] >= cut_value]
        label_ids       = df.loc[df['prediction'] == label]['id'].to_list()
        reference_fasta = Path(reference_fasta)
        fd = pd.DataFrame({ fid:fsq for fid,fsq in SimpleFastaParser(open(reference_fasta)) },index=[0]).T.rename_axis('id').reset_index()
        fd = fd.loc[fd['id'].isin(label_ids)].to_records(index=False)
        if label_ids:
            print(f'    - Retrieving {label} sequences...')
            with open(out_fasta, 'w') as sd:
                for fid, fsq in fd:
                    record = f'>{fid}|{label}\n{fsq}\n'
                    sd.write(record)


def te_count(dataframe: pd.DataFrame, mod: str) -> pd.DataFrame:
    '''Return dataframe with counts for each label.'''
    df      = pd.read_table(dataframe)
    counts  = df[['prediction_3', 'prediction_final']].value_counts().reset_index(name='count')
    if mod == "a":
        bases           = df.groupby('prediction_final')['length'].sum().reset_index(name='base_count')
        counts          = pd.merge(counts, bases, on='prediction_final')
        counts['id']    = counts['prediction_3'] + '|' + counts['prediction_final']
        counts          = counts.loc[:, ['id', 'count', 'base_count']]
    else:
        counts['id']    = counts['prediction_3'] + '|' + counts['prediction_final']
        counts          = counts.loc[:, ['id', 'count']]

    return counts


def concat_pred_tables(dfs, merged_0102, stp05, stp06, stp04, mod) -> None:
    """Concatenate final prediction tables in one"""
    final_dfs   = []
    merg        = pd.read_table(merged_0102)
    for ft in [stp05, stp06, stp04]:
        if ft.exists():
            df = table_filter(pd.read_table(ft))
            final_dfs.append(df)
    final_dfs = pd.concat(final_dfs)
    final_dfs = merg.merge(final_dfs, on='id', suffixes=['_3', '_final'])

    if mod == "a":
        fids            = final_dfs['id'].str.rsplit(':', n=1).str[0]
        start_end_col   = final_dfs['id'].str.rsplit(':', n=1).str[1]
        final_dfs['id'] = fids
        final_dfs.insert(1, 'start-end', start_end_col)
        final_dfs.insert(2, 'length', -(final_dfs['start-end'].map(eval)-1))

    final_dfs.to_csv(dfs, index=False, sep='\t')


def concat_fastas(final_fasta, ltr, nonltr, dna) -> None:
    """Concatenate final fastas files in one representing the TE library"""
    with open(final_fasta, 'wb') as wfd:
        for f in [ltr, nonltr, dna]:
            if f.exists():
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
