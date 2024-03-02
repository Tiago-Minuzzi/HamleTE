# HamleTE: a deep learning powered tool to annotate and classify transposable elements

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
    - [Count table](#counts)

- [Questions, issues and requests](#issues)

- [Updating](#updates)

    <img title="HamleTE" src="https://github.com/Tiago-Minuzzi/HamleTE/blob/dev/img/hamleTE.png" width="200" height="200" align="right">

## Latest updates

```
- 2024-03-02: New non-LTR model - classification between LINE and SINE"

- 2024-02-13: New ClassI/ClassII model.

- 2024-02-07: New TE/non-TE model.

- 2024-02-02: Clustering software changed to vsearch from cd-hit-est;no cluster by default.

- 2024-02-01: Python version to 3.10.12; Tensorflow to 2.13; Models update; k-mer length to 14.

- 2023-10-24: Total base counts by TE group added to CNT table.

- 2023-10-24: Function to replace non "ACTGN" bases optimized.

- 2023-09-08: Video tutorials for installation and usage.

- 2023-08-04: Default cutoff value set to 0.5.

- 2023-08-04: Prediction table containing the accuracy for all classification levels.
```

## Introduction<a name="intro"></a>

HamleTE is a deep learning-based tool with a workflow for finding and classifying transposable elements (TEs) in eukaryotic genomes.
It uses [Red](https://github.com/BioinformaticsToolsmith/Red) to find genomic repeats and, by using the power of convolutional neural networks feature extraction, 6 models to classify sequences as either being a TE or not, and then, the ones classified as TEs to the level of superfamily.

<img title="" src="https://github.com/Tiago-Minuzzi/HamleTE/blob/dev/img/workflow.png" width="600">

## Install <a name="installation"></a>

HamleTE can be installed either by creating a conda environment or manually. The first step is to [download](https://github.com/Tiago-Minuzzi/HamleTE/archive/refs/heads/main.zip) or clone this repo. To clone it run:

`git clone --depth 1 https://github.com/Tiago-Minuzzi/HamleTE`

Decompress the `hamlete_models.tar.xz`  in the `models` directory using your favorite application or via command-line using:

`tar xJvf hamlete_models.tar.xz`

After that, you can install the dependecies through a conda environment or manually.

Here is a video tutorial to install HamleTE: [Installation video](https://youtu.be/atULrCs3Pow).

### Conda environment <a name="conda"></a>

If you don't have conda installed, you can check how to install on [miniconda's webpage](https://docs.conda.io/en/latest/miniconda.html) or you can watch the installation tutorial for Linux here: [Video tutorial link](https://youtu.be/KI3yrW6VJIc). With conda installed on your system you can easily create a conda environment containing all dependecies by running:

`conda env create -f hamlete_env.yml`

Then, you can enter the conda environment with the command:

`conda activate hamleTE`

Download conda for linux [clicking here](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh).

### Docker <a name="docker"></a>

To run HamleTE using a docker container, first build the image using the Dockerfile. Inside HamleTE's directory run:

`docker build -t hamlete .`

### Manually <a name="depends"></a>

If you prefer a conda-free installation, it can be done manually by installing the depencies below:

- Python=3.10.12
- biopython=1.81
- h5py=3.9.0
- Keras=2.13.1
- numpy=1.24.3
- pandas=1.3.4
- scikit-learn=1.2.2
- scipy=1.10.1
- tensorflow=2.13.0
- vsearch=2.27.0
- tomli=2.0.1
- tqdm=4.64.1
- protobuf=4.24.0
- Red=2.0

For the manual installation, it's suggest to use a Python version management tool such as [Pyenv](https://github.com/pyenv/pyenv) and use it through a virtual environment to avoid depency conflicts. You can run `pip install -r requirements.txt` to install the Python packages needed.

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
                        between 0 and 1. Default = 0.9.
  -k LABEL_CUTOFF, --label_cutoff LABEL_CUTOFF
                        Cutoff value for TE superfamily classification. Value
                        must be between 0 and 1. Default = 0.5.
  -b BATCH_VALUE, --batch_value BATCH_VALUE
                        Set batch size. Default = 32, max = 250.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Set output directory to save results.
  -l LEN_KMER, --len_kmer LEN_KMER
                        Length of k-mer to find repeats in genomes. Default =
                        14.
  --min_len MIN_LEN     Minimum repeat sequence length. Default = 200.
  --clust               Cluster repeats. Slows down analysis cosiderably,
                        but reduces redundancy.
  --nobar               Disable progress bar.
```

### Basic usage <a name="basic"></a>

After activating HamleTE's conda environment (using `conda activate hamlete`), for genomes or transcriptomes, you can simply run:

`python3 hamleTE.py -f genome.fasta`

Video tutorial: [Running annotation mode](https://youtu.be/Ubl5oaI-HBw).

Clustering of repeats is disabled by default. If you would like to cluster the repeats to reduce redundancy, please, use the flag `--clust`. Example:

`python3 hamleTE.py -f genome.fasta --clust`

To use the classifier mode change the mode flag as follows:

`python3 hamleTE.py -m c -f my_TE_set.fasta`

Video tutorial: [Running classification mode](https://youtu.be/sEqVi2mauu8).

### Docker container <a name="udocker"></a>

To run the docker container version, mount the directory containing your fasta files inside the container using the `-v` flag.

`docker run -v /path/to/my/directory:/mnt -it hamlete hamleTE.py -f /mnt/genome.fasta -o /mnt/out_flow`

## Output example <a name="output"></a>

### Annotation mode <a name="annot"></a>

| id | start-end | length | prediction_1 | accuracy_1 | prediction_2 | accuracy_2 | prediction_3 | accuracy_3 | prediction_final | accuracy_final |
| :-------------: | :---------: | :------: | :------------: | :--------: | :------------: | :--------: | :------------: | :--------: | :------------: | :--------: | 
| chrom1 | 4852-4968 | 117 | TE | 0.999 | Retro | 0.998 | LTR | 1.0 | Gypsy | 0.809 |
| chrom2 | 88-1423 | 1336 | TE | 0.907 | Retro | 0.956 | nonLTR | 1.0 | LINE | 0.841 |
| chrom3 | 1-1906 | 1906 | TE | 0.983 | DNA | 0.994 | DNA | 0.994 | Tc1-Mariner | 0.952 |
| chrom4 | 1-1579 | 1579 | TE | 0.941 | DNA | 0.966 | DNA | 0.966 | Helitron | 0.979 |

### Classification mode <a name="class"></a>

| id | prediction_1 | accuracy_1 | prediction_2 | accuracy_2 | prediction_3 | accuracy_3 | prediction_final | accuracy_final |
| :-------------: | :------------: | :--------: | :------------: | :--------: | :------------: | :--------: | :------------: | :--------: |
| Seq_430 | TE | 1.0 | Retro | 0.582 | LTR | 0.516 | Gypsy | 0.999 |
| Seq_835 | TE | 0.792 | DNA | 0.89 | DNA | 0.89 | Tc1-Mariner | 1.0 |
| Seq_328 | TE | 0.966 | Retro | 1.0 | LTR | 1.0 | Copia | 0.705 |
| Seq_102 | TE | 0.99 | Retro | 0.9 | nonLTR | 1.0 | LINE | 0.966|

### Count table <a name="counts"></a>

The file ending with `CNT.tsv` contains the total count of TE groups and the total count of bases for each group in the annotation mode.

| id | count | base_count |
| :---: | :---: | :---: |
| LTR\|Gypsy | 166 | 414936 |
| nonLTR\|LINE | 106 | 318573 |
| DNA\|Helitron | 97 | 280230 |
| nonLTR\|Penelope | 51 | 75123 |

## Questions, issues and requests <a name="issues"></a>

If you have any questions about the usage, issues found during usage or feature requests, please, feel free to open an issue on the [issues section](https://github.com/Tiago-Minuzzi/HamleTE/issues) of HameleTE's github page.

## Updating <a name="updates"></a>

If you don't have the latest features, run the following command from the command-line inside HamleTE's folder on your machine:

`git pull`

---

## To-do

- [ ] Option to use `RepeatScount` instead of `Red`.

- [ ] Add Docker install and usage tutorial.

- [ ] Model for classifying LINEs in superfamilies.

- [ ] Add model for Class II subclass classification.

- [ ] Return log file.

- [ ] Add GPU support.
