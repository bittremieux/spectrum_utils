# Install

spectrum_utils, including all its required dependencies, can be easily
[installed using conda](https://anaconda.org/bioconda/spectrum_utils) from the
Bioconda channel:

    conda install -c conda-forge -c bioconda spectrum_utils

## Detailed installation information

### Supported Python versions

spectrum_utils supports Python version 3.6 and above.

### Dependencies

spectrum_utils has the following dependencies:

- [Altair](https://altair-viz.github.io/)
- [Matplotlib](https://matplotlib.org/)
- [Numba](http://numba.pydata.org/)
- [NumPy](https://www.numpy.org/)
- [Pandas](https://pandas.pydata.org/)
- [Pyteomics](https://pyteomics.readthedocs.io/)
- [RDKit](https://www.rdkit.org/)

Missing dependencies will be automatically installed when you install
spectrum_utils using conda.

## Alternative installation options

The recommended way to install spectrum_utils is using conda. Alternatively,
spectrum_utils can also be installed using pip:

    pip install spectrum_utils

To install the basic spectrum_utils version. Or:

    pip install spectrum_utils[iplot]

To include the interactive plotting functionality (requires Pandas and Altair).

When installing using pip it is recommended to explicitly install any
dependencies (listed above or in the
[environment file](https://github.com/bittremieux/spectrum_utils/blob/master/environment.yml))
in advance. Any missing dependencies will be automatically installed from PyPI
when you install spectrum_utils, _except_ RDKit. Please refer to the
[RDKit installation notes](https://www.rdkit.org/docs/Install.html) for
information on how to install RDKit.

In contrast, when installing spectrum_utils using conda all dependencies will
be automatically installed, including RDKit.
