# Spectrum processing

See the [quickstart](quickstart.md) for a brief introduction to how to start
using spectrum_utils. Here we will describe the spectrum processing
functionality provided by spectrum_utils in more detail.

## Peak annotations

Fragment ions can be annotated as follows:

- Using `MsmsSpectrum.annotate_peptide_fragments(...)` to annotate a, b, c, x,
y, or z ions for peptide spectra.
- Using `MsmsSpectrum.annotate_molecule_fragment(...)` using SMILES strings to
annotate peaks with molecule (sub)structures.
- Using `MsmsSpectrum.annotate_mz_fragment(...)` to annotate peaks with their
_m_/_z_ value.

Peak annotations can be visualized using the spectrum_utils plotting
functionality.

The example in the [quickstart](quickstart.md) shows how spectrum peaks can be
annotated with peptide fragments.

The following example shows how spectrum peaks can be annotated with their
_m_/_z_ values:

```python
import matplotlib.pyplot as plt
import pandas as pd
import spectrum_utils.spectrum as sus
import spectrum_utils.plot as sup
import urllib.parse

usi = 'mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00000840351'
peaks = pd.read_csv(
    f'https://metabolomics-usi.ucsd.edu/csv/?usi={urllib.parse.quote(usi)}')
precursor_mz = 633.2680
precursor_charge = 1
spectrum = sus.MsmsSpectrum(usi, precursor_mz, precursor_charge,
                            peaks['mz'].values, peaks['intensity'].values)
spectrum.filter_intensity(0.05)

tol_mass, tol_mode = 0.5, 'Da'
annotate_fragment_mz = [133.102, 147.080, 195.117, 237.164, 267.174, 295.170,
                        313.181, 355.192, 377.172, 391.187, 451.209, 511.231,
                        573.245, 633.269]
for fragment_mz in annotate_fragment_mz:
    spectrum.annotate_mz_fragment(fragment_mz, tol_mass, tol_mode)

fig, ax = plt.subplots(figsize=(12, 6))
sup.spectrum(spectrum, ax=ax)
plt.show()
plt.close()
```

Resulting in the following spectrum plot:

![](annotations.png)

## Peptide modifications

When annotating peptide spectra the masses of the fragment ions will be
automatically calculated to annotate the corresponding peaks. This
functionality is modification-aware: you can specify static modifications to
globally change the mass of specific amino acid residues and variable
modifications to modify amino acids in specific positions for a given peptide.

### Static modifications

Static modifications can be set as follows:

```python
sus.static_modification('C', 57.02146)
```

To, for example, set a static carbamidomethylation of cysteine.

Modification mass differences can be either positive or negative.

All static modifications can be reset:

```python
sus.reset_modifications()
```

### Variable modifications

Variable modifications can be set for an individual spectrum and peptide by
specifying the amino acid index and corresponding mass difference for each
modification.

Modification positions can be the following:

- The index of the amino acid where the modification is present (0-based).
- 'N-term' for N-terminal modifications.
- 'C-term' for C-terminal modifications.

For example, for a spectrum corresponding to peptide 'DLTDYLMK' with an
oxidated methionine:

```python
peptide = 'DLTDYLMK'
modifications = {6: 15.994915}

spectrum = sus.MsmsSpectrum(
    identifier, precursor_mz, precursor_charge, mz, intensity,
    peptide=peptide, modifications=modifications)
```


## Neutral losses

Besides the canonical a, b, c, x, y, and z ions, each of these ions can also be
automatically annotated with a neutral loss (or gain). Neutral losses need to
be specified by their identifier (molecular formula) and their mass difference:

```python
spectrum.annotate_peptide_fragments(0.05, 'Da', ion_types='aby',
                                    neutral_losses={'NH3': -17.026549,
                                                    'H2O': -18.010565})
```

The above example will consider all peaks with an optional ammonia (NH3) or
water (H2O) neutral loss.

Common neutral losses to consider are:

| Neutral loss/gain | Molecular formula | Mass difference |
| --- | --- | --- |
| Hydrogen | H | 1.007825 |
| Ammonia | NH3 | 17.026549 |
| Water | H2O | 18.010565 |
| Carbon monoxide | CO | 27.994915 |
| Carbon dioxide | CO2 | 43.989829 |
| Formamide | HCONH2 | 45.021464 |
| Formic acid | HCOOH | 46.005479 |
| Methanesulfenic acid | CH4OS | 63.998301 |
| Sulfur trioxide | SO3 | 79.956818 |
| Metaphosphoric acid | HPO3 | 79.966331 |
| Mercaptoacetamide | C2H5NOS | 91.009195 |
| Mercaptoacetic acid | C2H4O2S | 91.993211 |
| Phosphoric acid | H3PO4 | 97.976896 |

Note that typically the neutral _loss_ mass difference should be negative.
