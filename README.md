# HamleTE: a deep learning powered tool to find and classify transposable elements

## Table of contents

- [Introduction](#intro)

- [Install](#installation)
  
  - [Conda environment](#conda)
  - [Docker](#docker)
  - [Manually](#depends)

- [Usage](#usage)
  
  - [Basic usage](#basic)
  - [Docker container](#udocker)
    
- [Output example](#output)
    - [Annotaion mode](#annot)
    - [Classification mode](#class)

    <img title="HamleTE" src="https://github.com/Tiago-Minuzzi/HamleTE/blob/dev/img/hamleTE.png" width="200" height="200" align="right">

## Introduction<a name="intro"></a>

HamleTE is a deep learning-based tool with a workflow for finding and classifying transposable elements (TEs) in eukaryotic genomes.
It uses [Red](https://github.com/BioinformaticsToolsmith/Red) to find genomic repeats and, by using the power of convolutional neural networks feature extraction, 6 models to classify sequences as either being a TE or not, and then, the ones classified as TEs to the level of superfamily.

<img title="" src="https://github.com/Tiago-Minuzzi/HamleTE/blob/dev/img/workflow.png" width="600">

## Install <a name="installation"></a>

HamleTE can be installed either by creating a conda environment or manually. The first step is to [download](https://github.com/Tiago-Minuzzi/HamleTE/archive/refs/heads/main.zip) or clone this repo. To clone it run:

`git clone https://github.com/Tiago-Minuzzi/HamleTE`

Decompress the `hamlete_models.tar.xz`  in the `models` directory using your favorite application or via command-line using:

`tar xJvf hamlete_models.tar.xz`

After that, you can install the dependecies through a conda environment or manually.

Here is a video tutorial to install HamleTE: [Installation video](https://youtu.be/atULrCs3Pow).

### Conda environment <a name="conda"></a>

If you don't have conda installed, you can check how to install on [miniconda's webpage](https://docs.conda.io/en/latest/miniconda.html) or you can watch the installation tutorial for Linux here: [Video tutorial link](https://youtu.be/KI3yrW6VJIc). With conda installed on your system you can easily create a conda environment containing all dependecies by running:

`conda env create -f hamlete_env.yml`

Then, you can enter the conda environment with the command:

`conda activate hamlete`

### Docker <a name="docker"></a>

To run HamleTE using a docker container, first build the image using the Dockerfile. Inside HamleTE's directory run:

`docker build -t hamlete .`

### Manually <a name="depends"></a>

If want a conda-free installation, it can be done manually by installing the depencies below:

- Python=3.7
- biopython=1.79
- h5py=2.10.0
- Keras=2.3.1
- numpy=1.19.5
- pandas=1.3.4
- scikit-learn=1.0.1
- scipy=1.4.1
- tensorflow=2.1.0
- cd-hit-est=4.8.1
- tomli=2.0.1
- tqdm=4.64.1
- protobuf=3.20.0
- Red=2.0

For the manual installation, it's suggest to use a Python version management tool such as [Pyenv](https://github.com/pyenv/pyenv) and use it through a virtual environment to avoid depency conflicts. You can run `pip install -r requirements.txt` to install the Python packages needed.

To install cd-hit, you need to [download it](https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz) and add the cd-hit directory [to your path](https://linuxize.com/post/how-to-add-directory-to-path-in-linux/).

To install Red, clone [Red's github repository](https://github.com/BioinformaticsToolsmith/Red), change the name of your `C++` compiler inside Red's makefile and [compile the program](https://github.com/BioinformaticsToolsmith/Red/blob/master/src_2.0/HowToCompile.txt).

## Usage <a name="usage"></a>

The `annotation` mode is the default, which is used to find TE's in genomes or transcriptomes. If you have a set of sequences/TE library that you just want to classify, you can use the `classifier` mode by changing the mode flag. Below are the available options.

```
usage: hamleTE.py [-h] -f FASTA [-m MODE] [-c CUTOFF] [-k LABEL_CUTOFF]
                 [-b BATCH_VALUE] [-o OUTPUT_DIR] [-l LEN_KMER] [--noclust]
                 [--nobar]

Find repeats in eukaryotic genomes and classify them using deep learning.

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Genome or repeats/TEs fasta file.
  -m MODE, --mode MODE  Type (without quotation marks) 'a' for annotation mode or
                        'c' for classifier mode. Default = a.
  -c CUTOFF, --cutoff CUTOFF
                        Cutoff value for TE identification. Value must be
                        between 0 and 1. Default = 0.5.
  -k LABEL_CUTOFF, --label_cutoff LABEL_CUTOFF
                        Cutoff value for TE superfamily classification. Value
                        must be between 0 and 1.
  -b BATCH_VALUE, --batch_value BATCH_VALUE
                        Set batch size. Default = 32, max = 250.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Set output directory to save results.
  -l LEN_KMER, --len_kmer LEN_KMER
                        Length of k-mer to find repeats in genomes. Default =
                        13.
  --min_len MIN_LEN     Minimum repeat sequence length. Default = 50.
  --noclust             Do not cluster repeats. Results on more repeats to be
                        classified as TEs.
  --nobar               Disable progress bar.
```

### Basic usage <a name="basic"></a>

After activating HamleTE's conda environment, for genomes or transcriptomes, you can simply run:

`python3 hamleTE.py -f genome.fasta`

Video tutorial: [Running annotation mode](https://youtu.be/Ubl5oaI-HBw).

To use the classifier mode change the mode flag as follows:

`python3 hamleTE.py -m c -f my_TE_set.fasta`

Video tutorial: [Running classification mode](https://youtu.be/sEqVi2mauu8).

### Docker container <a name="udocker"></a>

To run the docker container version, mount the directory containing your fasta files inside the container using the `-v` flag.

`docker run -v /path/to/my/directory:/mnt -it hamlete hamleTE.py -f /mnt/genome.fasta -o /mnt/out_flow`

## Example output <a name="output"></a>

### Annotation mode <a name="annot"></a>

| id            | start-end | prediction   | accuracy |
| ------------- | --------- | ------------ | -------- |
| chrom_0100472 | 531-1258  | DNA\|Mutator | 0.983    |
| chrom_0079084 | 3410-4583 | nonLTR\|LINE | 0.828    |
| chrom_0330414 | 4081-4836 | LTR\|Gypsy   | 0.897    |
| chrom_0089582 | 1-175     | DNA\|hAT     | 0.895    |
| chrom_0071622 | 4394-4493 | LTR\|Copia   | 0.907    |

### Classification mode <a name="class"></a>

| id    | prediction    | accuracy |
| ----- | ------------- | -------- |
| seq_1 | nonLTR\|DIRS  | 0.808    |
| seq_2 | LTR\|Bel-Pao  | 0.943    |
| seq_3 | DNA\|Helitron | 0.835    |
| seq_4 | nonLTR\|SINE  | 0.843    |
| seq_5 | LTR\|Copia    | 0.705    |

---

## To-do

- [x] Return non-TE table

- [x] Return the probability value for the Class/Order

- [ ] Add GPU support
