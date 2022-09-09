# FlowTE

FlowTE is a tool with a workflow to find and classify transposable elements (TEs) in eukaryotic genomes using deep learning.
It utilizes [Red](https://github.com/BioinformaticsToolsmith/Red) to find genomic repeats and, using the power of convolutional neural networks feature extraction, 6 models to classify sequences as either being a TE or not, and then, the ones classified as TEs to the level of super family.

- [Installation](#installation)
	- [Conda](#conda)
	- [Manually](#depends)
## Installation <a name="installation"></a>

FlowTE can be installed either by creating a conda environment or manually.

### Conda <a name="conda"></a>

`conda env create -f flowte_env.yml`

Then, you can enter the conda environment running the command:

`conda activate flowte`

### Manually<a name="depends"></a>

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
- Red (see on Red's github page)





