# FlowTE

FlowTE is a tool with a workflow to find and classify transposable elements (TEs) in eukaryotic genomes using deep learning.
It utilizes [Red](https://github.com/BioinformaticsToolsmith/Red) to find genomic repeats and, using the power of convolutional neural networks feature extraction, 6 models to classify sequences as either being a TE or not, and then, the ones classified as TEs to the level of super family.
