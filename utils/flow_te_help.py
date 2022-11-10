import argparse


def restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{x} not a floating-point literal.")

    if not 0.0 < x < 1.0:
        raise argparse.ArgumentTypeError(f"{x} not in range [0.0, 1.0].")
    return x


parser = argparse.ArgumentParser(prog = "flowTE.py", description = "Find repeats in eukaryotic genomes and classify them using deep learning.")

parser.add_argument('-f', '--fasta',
                    help = 'Genome or repeats/TEs fasta file.',
                    required = True)

parser.add_argument('-m','--mode',
                    type = str,
                    default = 'g',
                    help = "Type (without quotation marks) 'g' for genome mode or 'c' for classifier mode (default = g).")

parser.add_argument('-c', '--cutoff',
                    type = restricted_float,
                    default = None,
                    help = "Cutoff value for TE identification. Value must be between 0 and 1.")

parser.add_argument('-k', '--label_cutoff',
                    type = restricted_float,
                    default = None,
                    help = "Cutoff value for TE superfamily classification. Value must be between 0 and 1.")

# parser.add_argument('-b','--batch_value',
#                     type = int,
#                     default = 32,
#                     help = "Set batch size (Default = 32).")

parser.add_argument('-o','--output_dir',
                    type = str,
                    default = '.',
                    help = "Set output directory to save results.")

args = parser.parse_args()
