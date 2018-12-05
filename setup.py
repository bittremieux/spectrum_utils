import setuptools

import spectrum_utils


DISTNAME = 'spectrum_utils'
VERSION = spectrum_utils.__version__
DESCRIPTION = 'Mass spectrometry utility functions'
AUTHOR = 'Wout Bittremieux'
AUTHOR_EMAIL = 'wout.bittremieux@uantwerpen.be'
URL = 'https://github.com/bittremieux/spectrum_utils'
LICENSE = 'Apache 2.0'

setuptools.setup(
        name=DISTNAME,
        version=VERSION,
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        description=DESCRIPTION,
        url=URL,
        license=LICENSE,
        platforms=['any'],
        classifiers=[
            'Development Status :: 2 - Pre-Alpha',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache Software License',
            'Operating System :: MacOS',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: Unix',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3 :: Only',
            'Topic :: Scientific/Engineering :: Bio-Informatics'],
        packages=['spectrum_utils'],
        install_requires=[
            'matplotlib',
            'numpy',
            'pyteomics'],
)
