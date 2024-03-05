import sys
import subprocess
from pathlib import Path


def repeat_scout_runner(fasta: str, tmp_dir: str, outfasta: str) -> None:
    # Manipulate input info
    fasta = Path(fasta)
    base_name       = fasta.stem
    frequence_file  = tmp_dir / f"{base_name}.freq"
    output_fasta    = outfasta
    # Run build_lmer_table
    print("### Building l-mer table")
    build_table     = subprocess.run(["build_lmer_table",
                                      "-sequence",  fasta,
                                      "-freq",      frequence_file])
    # Run RepeatScout
    print("### Starting RepeatScout")
    rscout          = subprocess.run(["RepeatScout",
                                      "-sequence",  fasta,
                                      "-freq",      frequence_file,
                                      "-output",    output_fasta,
                                      "-v"])
