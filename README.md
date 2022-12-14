# FlowTE

FlowTE is a tool with a workflow to find and classify transposable elements (TEs) in eukaryotic genomes using deep learning.
It utilizes [Red](https://github.com/BioinformaticsToolsmith/Red) to find genomic repeats and, using the power of convolutional neural networks feature extraction, 6 models to classify sequences as either being a TE or not, and then, the ones classified as TEs to the level of super family.

- [Install](#installation)
  - [Conda environment](#conda)
  - [Docker](#docker)
  - [Manually](#depends)

## Install <a name="installation"></a>

FlowTE can be installed either by creating a conda environment or manually. The first step is to [download](https://github.com/Tiago-Minuzzi/FlowTE/archive/refs/heads/main.zip) or clone this repo. To clone it run:

`git clone https://github.com/Tiago-Minuzzi/FlowTE`

Decompress the `flowte_models.tar.xz`  in the `models` directory using your favorite application or via command-line using:

`tar xJvf flowte_models.tar.xz`

After that, you can install the dependecies through a conda environment or manually.

### Conda environment <a name="conda"></a>

If you don't have conda installed, you can check how to install on [miniconda's webpage](https://docs.conda.io/en/latest/miniconda.html). With conda installed on your system you can easily create a conda environment containing all dependecies by running:

`conda env create -f flowte_env.yml`

Then, you can enter the conda environment with the command:

`conda activate flowte`

### Docker <a name="docker"><\a>

To run FlowTE using a docker container, first build the image using the Dockerfile. Inside FlowTE's directory run:

`docker build -t flowte .`

### Manually <a name="depends"></a>

If want a conda-free installation, it can be done manually by installing the depencies below:

- Python=3.7.12
- biopython=1.79
- h5py=2.10.0
- Keras=2.3.1
- numpy=1.19.5
- pandas=1.3.4
- scikit-learn=1.0.1
- scipy=1.4.1
- tensorflow=2.1.0
- cd-hit-est=4.8.1
- Red

For the manual installation, it's suggest to use a Python version management tool such as [Pyenv](https://github.com/pyenv/pyenv) and use it through a virtual environment to avoid depency conflicts. You can run `pip install -r requirements.txt` to install the Python packages needed.

To install Red, you will need to clone [Red's github repository](https://github.com/BioinformaticsToolsmith/Red), change the name of your `C++` compiler inside Red's makefile and [compile the program](https://github.com/BioinformaticsToolsmith/Red/blob/master/src_2.0/HowToCompile.txt).
