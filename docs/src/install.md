# Install

spectrum_utils requires Python version 3.8+ and can be installed with pip or conda.

Using pip:

    pip install spectrum_utils[iplot]

Using conda:

    conda install -c bioconda spectrum_utils

## Supported Python versions

spectrum_utils supports Python version 3.8 and above.

## Dependencies

spectrum_utils has the following third-party dependencies:

- [fastobo](https://fastobo.readthedocs.io/)
- [Lark](https://lark-parser.readthedocs.io/)
- [Matplotlib](https://matplotlib.org/)
- [Numba](http://numba.pydata.org/)
- [NumPy](https://www.numpy.org/)
- [Pandas](https://pandas.pydata.org/)
- [platformdirs](https://github.com/platformdirs/platformdirs)
- [Pyteomics](https://pyteomics.readthedocs.io/)
- [Vega-Altair](https://altair-viz.github.io/)

Missing dependencies will be automatically installed when you install spectrum_utils using pip or conda.

Additionally, we recommend manually installing [pyteomics.cythonize](https://pypi.org/project/pyteomics.cythonize/) as a plug-in replacement for faster fragment ion mass calculations.

## Advanced installation instructions

spectrum_utils provides modular installation capabilities to minimize the number of third-party dependencies that will be installed when only a subset of the spectrum_utils functionality is required.
The previous pip and conda commands will install all optional spectrum_utils extensions (excluding developer and documentation dependencies).
Power users can customize their spectrum_utils installation by specifying one or more of the following sets of dependencies:

- `dev`: Developer dependencies for automatic linting and testing.
- `docs`: Dependencies to generate these documentation pages.
- `iplot`: Interactive spectrum plotting using [Vega-Altair](https://altair-viz.github.io/).
