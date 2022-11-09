import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser


def get_selected_sequences(in_fasta: str, prediction_table: str, out_fasta: str) -> None:
    if prediction_table.exists():
        df = pd.read_table(prediction_table)
        sf_labels = ', '.join(df['prediction'].unique().tolist())
        fids = df[['id','prediction']].values.tolist()
        print(f'    Retrieving sequences for {sf_labels}')
        with open(in_fasta) as fa, open(out_fasta,'w') as sd:
            for fid, fsq in SimpleFastaParser(fa):
                for i in fids:
                    if fid == i[0]:
                        record = f'>{fid}|{i[1]}\n{fsq}\n'
                        sd.write(record)
