==================
``spectrum_utils``
==================

.. image:: https://travis-ci.org/bittremieux/spectrum_utils.svg?master
    :target: https://travis-ci.org/bittremieux/spectrum_utils
.. image:: https://img.shields.io/badge/python-3.6-brightgreen.svg
.. image:: https://img.shields.io/badge/python-3.7-brightgreen.svg

Simple MS/MS spectrum preprocessing and visualization in Python.

Features
--------

- Spectrum (pre)processing
    - Precursor & noise peak removal
    - Intensity filtering
    - Intensity scaling
    - Modification-aware fragment ion annotating (powered by `Pyteomics <https://pyteomics.readthedocs.io/>`_)
- Spectrum plotting
    - Single spectrum with annotated ions
    - Mirror plot of matching spectra
    - Interactive spectrum plots

Example
-------

::

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


(Condensed example. See `here <https://github.com/bittremieux/spectrum_utils/blob/master/notebooks/quickstart.ipynb>`_ for the full code to generate the figure below.)

.. image:: spectrum_utils.png

Installation
------------

``spectrum_utils`` can be installed easily via pip::

    pip install spectrum_utils

Or via conda::

    conda install -c bioconda spectrum_utils

Dependencies
~~~~~~~~~~~~

``spectrum_utils`` has the following dependencies:

- `Altair <https://altair-viz.github.io/>`_
- `Matplotlib <https://matplotlib.org/>`_
- `Numba <http://numba.pydata.org/>`_
- `NumPy <https://www.numpy.org/>`_
- `Pandas <https://pandas.pydata.org/>`_
- `Pyteomics <https://pyteomics.readthedocs.io/>`_

Missing dependencies will be automatically installed when you install ``spectrum_utils``.

API documentation
-----------------

.. toctree::

    spectrum_utils

Contact
-------

For more information you can visit the `official code website <https://github.com/bittremieux/spectrum_utils/>`_ or send an `email <wout.bittremieux@uantwerpen.be>`_.
