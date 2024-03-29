#!/usr/bin/env python3
import os
import sys
import pathlib
import numpy as np
import pandas as pd
from orffinder import orffinder
from seqshannon import shannon_entropy
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
# Hide warnings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
stderr      = sys.stderr
sys.stderr  = open(os.devnull, 'w')
from keras.preprocessing.text import Tokenizer
sys.stderr  = stderr


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


def label_pred_dataframe(fasta_ids, prediction_results, colunas) -> pd.DataFrame:
    '''Return prediction values as dataframe.'''
    # Create dataframe
    label_pred_df   = pd.DataFrame(prediction_results,
                                   columns = colunas)
    label_pred_df   = label_pred_df.round(3)
    label_pred_df   = pd.concat([fasta_ids, label_pred_df], axis=1)
    # Create label column
    label_pred_df['prediction'] = label_pred_df[colunas].idxmax(axis=1)
    return label_pred_df


def filter_by_len(fasta: str, temp_dir: str, filtered: str) -> None:
    '''Filter sequences by length. Keeps only sequences between 200 and 25_000 bp.'''
    with open(fasta) as fa, open(filtered, "w") as sd:
        for fid, fsq in SimpleFastaParser(fa):
            if 200 <= len(fsq) <= 25_000:
                sd.write(f">{fid}\n{fsq}\n")


def calc_entropy(fasta: str, outfasta: str) -> None:
    '''Calculates sequence entropy to filter out low complexity sequences.'''
    with open(fasta) as fa, open(outfasta, "w") as sd:
        for fid, fsq in SimpleFastaParser(fa):
            entropia: float = shannon_entropy(fsq)
            if entropia >= 1.9:
                sd.write(f">{fid}\n{fsq}\n")


def orf_checker(fasta: str, outfasta: str) -> None:
    fasta       = pathlib.Path(fasta)
    basename    = fasta.stem
    with open(fasta) as fa, open(outfasta, "w") as sd:
        for record in SeqIO.parse(fa, 'fasta'):
            res: list = orffinder.getORFs(record)
            if res:
                SeqIO.write(record, sd, 'fasta')
