# spectrum_utils

[![conda](https://img.shields.io/conda/vn/bioconda/spectrum_utils?color=green)](http://bioconda.github.io/recipes/spectrum_utils/README.html)
[![PyPI](https://img.shields.io/pypi/v/spectrum_utils?color=green)](https://pypi.org/project/spectrum_utils/)
[![Build status](https://github.com/bittremieux/spectrum_utils/workflows/tests/badge.svg)](https://github.com/bittremieux/spectrum_utils/actions?query=workflow:tests)
[![docs](https://readthedocs.org/projects/spectrum_utils/badge/?version=latest)](https://spectrum_utils.readthedocs.io/en/latest/?badge=latest)

## About spectrum_utils

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
 
 See the documentation for more information and detailed examples on how to use
 this functionality.
 
 ## Citation
 
spectrum_utils is freely available as open source under the
[Apache 2.0 license](http://opensource.org/licenses/Apache-2.0).

When using spectrum_utils please cite the following manuscript:
 
Wout Bittremieux. "spectrum_utils: A Python package for mass spectrometry data
processing and visualization." _Analytical Chemistry_ 92 (1) 659-661 (2020)
doi:[10.1021/acs.analchem.9b04884](https://doi.org/10.1021/acs.analchem.9b04884).

## Contents

- [Install](install.md)
- [Quickstart](quickstart.md)
- [Spectrum processing](processing.md)
- [Spectrum visualization](plotting.md)
- [Computational efficiency](runtime.md)
- [API reference](api.md)
- [Contact](contact.md)
