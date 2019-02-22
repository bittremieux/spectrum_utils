# spectrum_utils

[![Build status](https://travis-ci.org/bittremieux/spectrum_utils.svg?master)](https://travis-ci.org/bittremieux/spectrum_utils)
![Python 3.6](https://img.shields.io/badge/python-3.6-brightgreen.svg)
![Python 3.7](https://img.shields.io/badge/python-3.7-brightgreen.svg)

Simple MS/MS spectrum preprocessing and visualization in Python.

## Example

```
import matplotlib.pyplot as plt

from spectrum_utils import plot
from spectrum_utils import spectrum


# Initialize spectrum information first...

spec = spectrum.MsmsSpectrum(identifier, precursor_mz, precursor_charge,
                             mz, intensity, retention_time=retention_time,
                             peptide=peptide)

# Preprocess the MS/MS spectrum.
fragment_tol_mass = 10
fragment_tol_mode = 'ppm'
spec = (spec.set_mz_range(min_mz=100, max_mz=1400)
            .remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
            .filter_intensity(min_intensity=0.05, max_num_peaks=150)
            .scale_intensity(scaling='root')
            .annotate_peaks(fragment_tol_mass, fragment_tol_mode,
                            ion_types='aby'))

# Plot the MS/MS spectrum.
plot.spectrum(spec)

plt.show()
plt.close()
```
(Condensed example. See [here](https://github.com/bittremieux/spectrum_utils/blob/master/notebooks/preprocess_and_plot.ipynb) for the full code.)

![spectrum_utils](spectrum_utils.png)
