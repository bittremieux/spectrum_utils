# Install

spectrum_utils requires Python version 3.6+ and can be installed with pip or
conda.

Using conda:

    conda install -c bioconda spectrum_utils

Using pip (default; only static plotting using Matplotlib):

    pip install spectrum_utils

To include interactive plotting functionality (requires Pandas and Altair) when
installing using pip:

    pip install spectrum_utils[iplot]

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

Missing dependencies will be automatically installed when you install
spectrum_utils using conda or pip.
