import argparse


def restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{x} not a floating-point literal.")

    if not 0.0 < x < 1.0:
        raise argparse.ArgumentTypeError(f"{x} not in range [0.0, 1.0].")
    return x


parser=argparse.ArgumentParser(prog="hamleTE.py", description="Find repeats in eukaryotic genomes and classify them using deep learning.")

parser.add_argument('-f', '--fasta',
                    help='Genome or repeats/TEs fasta file.',
                    required=True)

parser.add_argument('-m', '--mode',
                    type=str,
                    default='a',
                    help="Type (without quotation marks) 'a' for annotation mode or 'c' for classifier mode. Default = a.")

parser.add_argument('-c', '--cutoff',
                    type=restricted_float,
                    default=0.5,
                    help="Cutoff value for TE identification. Value must be between 0 and 1. Default = 0.5.")

parser.add_argument('-k', '--label_cutoff',
                    type=restricted_float,
                    default=None,
                    help="Cutoff value for TE superfamily classification. Value must be between 0 and 1.")

parser.add_argument('-b', '--batch_value',
                    type=int,
                    default=32,
                    help="Set batch size. Default = 32, max = 250.")

parser.add_argument('-o', '--output_dir',
                    type=str,
                    default='.',
                    help="Set output directory to save results.")

parser.add_argument('-l', '--len_kmer',
                    type=int,
                    default=13,
                    help="Length of k-mer to find repeats in genomes. Default = 13.")

parser.add_argument('--min_len',
                    type=int,
                    default=50,
                    help="Minimum repeat sequence length. Default=50.")

parser.add_argument('--noclust',
                    action='store_false',
                    help="Do not cluster repeats. Results on more repeats to be classified as TEs.")

parser.add_argument('--nobar',
                    action='store_true',
                    help="Disable progress bar.")

args=parser.parse_args()
