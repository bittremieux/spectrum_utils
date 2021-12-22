import operator
import re
from typing import Any, Dict, List, Optional, Tuple

import numba as nb
import numpy as np

try:
    import pyteomics.cmass as pmass
except ImportError:
    import pyteomics.mass as pmass

from spectrum_utils import proforma, utils


# Amino acid and special amino acid masses.
_aa_mass = {
    **pmass.std_aa_mass,
    # Aspartic acid / asparagine (ambiguous mass).
    # "B": 0,
    # Glutamic acid / glutamine (ambiguous mass).
    # "Z": 0,
    # Leucine / isoleucine.
    "J": 113.08406,
    # Selenocysteine (in Pyteomics).
    # "U": 150.95363,
    # Pyrrolysine (in Pyteomics).
    # "O": 237.14772,
    # Any amino acid, gaps (zero mass).
    "X": 0,
}

# Common neutral losses.
_neutral_loss = {
    # No neutral loss.
    None: 0,
    # Hydrogen.
    "H": -1.007825,
    # Ammonia.
    "NH3": -17.026549,
    # Water.
    "H2O": -18.010565,
    # Carbon monoxide.
    "CO": -27.994915,
    # Carbon dioxide.
    "CO2": -43.989829,
    # Formamide.
    "HCONH2": -45.021464,
    # Formic acid.
    "HCOOH": -46.005479,
    # Methanesulfenic acid.
    "CH4OS": -63.998301,
    # Sulfur trioxide.
    "SO3": -79.956818,
    # Metaphosphoric acid.
    "HPO3": -79.966331,
    # Mercaptoacetamide.
    "C2H5NOS": -91.009195,
    # Mercaptoacetic acid.
    "C2H4O2S": -91.993211,
    # Phosphoric acid.
    "H3PO4": -97.976896,
}


class PeakInterpretation:
    def __init__(self):
        """
        Fragment annotation(s) to interpret a specific peak.
        """
        self.annotations = []

    def __str__(self):
        # If no fragment annotations have been specified, interpret as an
        # unknown ion.
        return ",".join(self.annotations) if len(self.annotations) > 0 else "?"


class FragmentAnnotation:
    def __init__(
        self,
        ion_type: str,
        neutral_loss: Optional[str] = None,
        isotope: int = 0,
        charge: Optional[int] = None,
        adduct: Optional[str] = None,
        analyte_number: Optional[int] = None,
        mz_delta: Optional[Tuple[float, str]] = None,
    ) -> None:
        """
        Individual fragment ion annotation.

        This fragment annotation format is derived from the PSI peak
        interpretation specification:
        https://docs.google.com/document/d/1yEUNG4Ump6vnbMDs4iV4s3XISflmOkRAyqUuutcCG2w/edit?usp=sharing

        Fragment notations have the following format:

        (analyte_number)[ion_type](neutral_loss)(isotope)(charge)(adduct)(mz_delta)

        Examples:

        - "y4-H2O+2i^2[M+H+Na]" : Fragment annotation for a y4 ion, with a
          water neutral loss, the second isotopic peak, charge 2, adduct
          [M+H+Na].

        Parameters
        ----------
        ion_type : str
            Specifies the basic type of ion being described.
            Possible prefixes are:

            - "?": unknown ion
            - "a", "b", "c", "x", "y", "z": corresponding peptide fragments
            - "I": immonium ion
            - "m": internal fragment ion
            - "_": named compound
            - "p": precursor ion
            - "r": reporter ion (isobaric label)
            - "f": chemical formula
        neutral_loss : Optional[str]
            A string of neutral loss(es), described by their molecular formula.
            The default is no neutral loss. Note that the neutral loss string
            must include the sign (typically "-" for a neutral loss).
        isotope : int
            The isotope number above or below the monoisotope. The default is
            the monoisotopic peak (0).
        charge : Optional[int]
            The charge of the fragment. The default is an unknown charge (only
            valid for unknown ions).
        adduct : Optional[str]
            The adduct that ionized the fragment. The default is a hydrogen
            adduct matching the charge ([M+xH]).
        mz_delta : Optional[Tuple[float, str]]
            The m/z delta representing the observed m/z minus the theoretical
            m/z and its unit ("Da" or "ppm").
        """
        if ion_type[0] in "GLXS":
            raise NotImplementedError(
                "Advanced ion types are not yet supported"
            )
        elif ion_type[0] not in "?abcxyzIm_prf":
            raise ValueError("Unknown ion type")
        if ion_type == "?" and (
            neutral_loss is not None
            or isotope != 0
            or charge is not None
            or adduct is not None
            or analyte_number is not None
            or mz_delta is not None
        ):
            raise ValueError(
                "Unknown ions should not contain additional information"
            )
        self.ion_type = ion_type
        self.neutral_loss = neutral_loss
        self.isotope = isotope
        self.charge = charge
        self.adduct = f"[M+{self.charge}H]" if adduct is None else adduct
        self.analyte_number = analyte_number
        self.mz_delta = mz_delta

    @property
    def mz_delta(self) -> Optional[Tuple[float, str]]:
        return self._mz_delta

    @mz_delta.setter
    def mz_delta(self, mz_delta: Optional[Tuple[float, str]]):
        if mz_delta is not None and mz_delta[1] not in ("Da", "ppm"):
            raise ValueError(
                "The m/z delta must be specified in Dalton or ppm units"
            )
        self._mz_delta = mz_delta

    @property
    def charge(self) -> Optional[int]:
        return self._charge

    @charge.setter
    def charge(self, charge: Optional[int]):
        if self.ion_type == "?" and charge is not None:
            raise ValueError("Invalid charge for unknown ions")
        elif self.ion_type != "?" and (charge is None or charge <= 0):
            raise ValueError(
                "The charge must be specified and strictly positive for known "
                "ion types"
            )
        self._charge = charge

    def __str__(self) -> str:
        if self.ion_type == "?":
            return "?"
        else:
            annot_str = []
            if self.analyte_number is not None:
                annot_str.append(f"{self.analyte_number}@")
            annot_str.append(self.ion_type)
            if self.neutral_loss is not None:
                annot_str.append(self.neutral_loss)
            if abs(self.isotope) == 1:
                annot_str.append("+i" if self.isotope > 0 else "-i")
            elif self.isotope != 0:
                annot_str.append(f"{self.isotope:+}i")
            if self.charge is not None and self.charge > 1:
                annot_str.append(f"^{self.charge}")
            if re.match(r"\[M\+\d+H\]", self.adduct) is not None:
                annot_str.append(self.adduct)
            if self.mz_delta is not None:
                annot_str.append(
                    f"/{self.mz_delta[0]}"
                    f"{'ppm' if self.mz_delta[1] == 'ppm' else ''}"
                )
            return "".join(annot_str)

    def __eq__(self, other: Any) -> bool:
        return isinstance(other, FragmentAnnotation) and str(self) == str(
            other
        )
