# Spectrum annotating

See the [quickstart](quickstart.md) for a brief introduction on how to start using spectrum_utils.
Here we will describe the spectrum annotation functionality provided by spectrum_utils in more detail.

## Fragment ion annotation

As demonstrated in the [quickstart](quickstart.md), fragment ions can be annotated based on the [ProForma 2.0](https://www.psidev.info/proforma) specification.

**TODO:** Brief introduction and a few examples of different ways to encode peptidoforms using ProForma.

## Ion types

During fragment ion annotation, by default peptide b and y ions will be annotated.
Additionally, spectrum_utils supports several other types of fragment ions:

- **TODO:** List of supported ion types.

## Neutral losses

Each of the above ions can also be automatically considered with a neutral loss (or gain).
Neutral losses need to be specified by their identifier (molecular formula) and their mass difference:

```python
spectrum.annotate_proforma(
    "WNQLQAFWGTGK",
    fragment_tol_mass=0.05,
    fragment_tol_mode="Da",
    ion_types="aby",
    neutral_losses={"NH3": -17.026549, "H2O": -18.010565},
)
```

The above example will consider all peaks with an optional ammonia (NH3) or water (H2O) neutral loss.

Common neutral losses that can be used are:

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
