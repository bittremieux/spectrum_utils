# spectrum_utils

[![Build status](https://travis-ci.org/bittremieux/spectrum_utils.svg?master)](https://travis-ci.org/bittremieux/spectrum_utils)
![Python 3.6](https://img.shields.io/badge/python-3.6-brightgreen.svg)
![Python 3.7](https://img.shields.io/badge/python-3.7-brightgreen.svg)
![Python 3.8](https://img.shields.io/badge/python-3.8-brightgreen.svg)

spectrum_utils is a Python package for efficient MS/MS spectrum processing and
visualization.

spectrum_utils contains the following features:

- Spectrum processing
    - Precursor & noise peak removal
    - Intensity filtering
    - Intensity scaling
    - Peak annotations
        - Modification-aware (static & variable) peptide fragments
        - SMILES-based molecules
        - Custom strings
- Spectrum plotting
    - Fully customizable individual spectrum plots
    - Mirror plot of matching spectra
    - Interactive spectrum plots

![spectrum_utils](spectrum_utils.png)

## Installation

spectrum_utils, including all its required dependencies, can be easily
[installed using conda](https://anaconda.org/bioconda/spectrum_utils) from the
Bioconda channel:

    conda install -c conda-forge -c bioconda -c defaults spectrum_utils

## Documentation

Please see the [documentation](https://spectrum-utils.readthedocs.io/) for
detailed installation instructions, usage examples, the API reference, and more
information.

## Citation
 
spectrum_utils is freely available as open source under the
[Apache 2.0 license](http://opensource.org/licenses/Apache-2.0).

When using spectrum_utils please cite the following manuscript:
 
Wout Bittremieux. "spectrum_utils: A Python package for mass spectrometry data
processing and visualization." _Analytical Chemistry_ 92 (1) 659-661 (2020)
doi:[10.1021/acs.analchem.9b04884](https://doi.org/10.1021/acs.analchem.9b04884).
