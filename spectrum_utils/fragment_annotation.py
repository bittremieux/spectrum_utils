import operator
from typing import Any, Dict, List, Optional, Tuple

import numba as nb
import numpy as np

try:
    import pyteomics.cmass as pmass
except ImportError:
    import pyteomics.mass as pmass

from spectrum_utils import proforma, utils


# Amino acid and special amino acid masses.
aa_mass = {
    **pmass.std_aa_mass,
    # "B": 0,             # aspartic acid / asparagine (ambiguous mass)
    # "Z": 0,             # glutamic acid / glutamine (ambiguous mass)
    "J": 113.08406,       # leucine / isoleucine
    # "U": 150.95363,     # selenocysteine (in Pyteomics)
    # "O": 237.14772,     # pyrrolysine (in Pyteomics)
    "X": 0,               # any amino acid, gaps (zero mass)
}


class FragmentAnnotation:
    def __init__(
        self,
        ion_type: str,
        neutral_loss: Optional[str] = None,
        isotope: int = 0,
        charge: int = 0,
        adduct: Optional[str] = None,
        calc_mz: float = None,
    ) -> None:
        """
        Interpretation of a single fragment ion.

        This fragment annotation format is derived from the PSI peak
        interpretation specification:
        https://docs.google.com/document/d/1yEUNG4Ump6vnbMDs4iV4s3XISflmOkRAyqUuutcCG2w/edit?usp=sharing
        Currently a simplified subset of this specification is supported.

        Ion notations have the following format:

        [ion type](neutral loss)(isotope)(charge)(adduct type)

        Examples:

        - "y4-H2O+2i^2[M+H+Na]" : Fragment annotation for a y4 ion, with a
          water neutral loss, the second isotopic peak, charge 2, adduct
          M+H+Na.

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
        neutral_loss : str, optional
            A string of neutral loss(es), described by their molecular formula.
            The default is no neutral loss.
        isotope : int, optional
            The isotope number above or below the monoisotope. The default is
            the monoisotopic peak (0).
        charge : int, optional
            The charge of the fragment. The default is charge 0 (for unknown
            ions).
        adduct : str, optional
            The adduct that ionized the fragment.
        calc_mz : float
            The theoretical m/z value of the fragment.
        """
        if ion_type[0] not in "?abcxyzIm_prf":
            raise ValueError("Unsupported ion type")
        if ion_type == "?" and (
            neutral_loss is not None
            or isotope != 0
            or charge != 0
            or adduct is not None
        ):
            raise ValueError(
                "No information should be specified for unknown ions"
            )
        self.ion_type = ion_type
        self.neutral_loss = neutral_loss
        self.isotope = isotope
        if charge == 0 and ion_type != "?":
            raise ValueError("Invalid charge 0 for annotated fragment")
        elif charge < 0:
            raise ValueError("Invalid negative charge")
        self.charge = charge
        self.adduct = adduct
        self.calc_mz = calc_mz

    def __repr__(self) -> str:
        if self.ion_type == "?":
            fragment_repr = "?"
        else:
            fragment_repr = self.ion_type
            if self.neutral_loss is not None:
                fragment_repr += self.neutral_loss
            if self.isotope != 0:
                fragment_repr += f"{self.isotope:+}i"
            if self.charge > 1:
                fragment_repr += f"^{self.charge}"
            if self.adduct is not None:
                fragment_repr += self.adduct
        return f"FragmentAnnotation({fragment_repr}, mz={self.calc_mz})"

    def __str__(self) -> str:
        if self.ion_type == "?":
            return str(self.calc_mz)
        else:
            annotation = self.ion_type
            if self.neutral_loss is not None:
                annotation += f"{self.neutral_loss}"
            if self.isotope != 0:
                annotation += f"{self.isotope:+}i"
            annotation += "+" * self.charge
            if self.adduct is not None:
                annotation += self.adduct
            return annotation

    def __eq__(self, other: Any) -> bool:
        return isinstance(other, FragmentAnnotation) and repr(self) == repr(
            other
        )


def get_theoretical_fragments(
    proteoform: proforma.Proteoform,
    fragment_types: str = "by",
    max_charge: int = 1,
    neutral_losses: Optional[Dict[Optional[str], float]] = None,
) -> List[FragmentAnnotation]:
    """
    Get theoretical fragment annotations for the given sequence.

    Parameters
    ----------
    proteoform : proforma.Proteoform
        The proteoform for which the fragment annotations will be generated.
    fragment_types : str, optional
        The peptide fragment type. Can be any combination of 'a', 'b', 'c',
        'x', 'y', and 'z' (the default is 'by', which means that b-ions and
        y-ions will be generated).
    max_charge : int, optional
        All fragments up to and including the given charge will be generated
        (the default is 1 to only generate singly-charged fragments).
    neutral_losses : Optional[Dict[Optional[str], float]]
        A dictionary with neutral loss names and (negative) mass differences to
        be considered.

    Returns
    -------
    List[FragmentAnnotation]
        A list of all theoretical fragments in ascending m/z order.
    """
    if "B" in proteoform.sequence:
        raise ValueError(
            "Explicitly specify aspartic acid (D) or asparagine (N) instead of"
            " the ambiguous B to compute the fragment annotations"
        )
    if "Z" in proteoform.sequence:
        raise ValueError(
            "Explicitly specify glutamic acid (E) or glutamine (Q) instead of "
            "the ambiguous Z to compute the fragment annotations"
        )

    neutral_losses = {None: 0} if neutral_losses is None else neutral_losses
    fragments = []
    # FIXME
    # # Single peak annotation.
    # if sequence == 'X':
    #     return [FragmentAnnotation('?', calc_mz=modifications[0])]
    # Get all possible peptide fragments.
    for i in range(1, len(proteoform.sequence)):
        for fragment_type in fragment_types:
            # N-terminal fragment.
            if fragment_type in "abc":
                fragment_index = i
                fragment_sequence = proteoform.sequence[:i]
                if proteoform.modifications is not None:
                    mod_mass = sum(
                        [
                            mod.mass
                            for mod in proteoform.modifications
                            if (
                                isinstance(mod.position, int)
                                and mod.position < i
                            )
                            or mod.position == "N-term"
                        ]
                    )
                else:
                    mod_mass = 0
            # C-terminal fragment.
            elif fragment_type in "xyz":
                fragment_index = len(proteoform.sequence) - i
                fragment_sequence = proteoform.sequence[i:]
                if proteoform.modifications is not None:
                    mod_mass = sum(
                        [
                            mod.mass
                            for mod in proteoform.modifications
                            if (
                                isinstance(mod.position, int)
                                and mod.position >= i
                            )
                            or mod.position == "C-term"
                        ]
                    )
                else:
                    mod_mass = 0
            else:
                raise ValueError(
                    f"Unknown/unsupported ion type: {fragment_type}"
                )
            for charge in range(1, max_charge + 1):
                for nl_name, nl_mass in neutral_losses.items():
                    nl_name = (
                        None
                        if nl_name is None
                        else f'{"-" if nl_mass < 0 else "+"}{nl_name}'
                    )
                    fragments.append(
                        FragmentAnnotation(
                            ion_type=f"{fragment_type}{fragment_index}",
                            neutral_loss=nl_name,
                            isotope=0,
                            charge=charge,
                            calc_mz=pmass.fast_mass(
                                sequence=fragment_sequence,
                                ion_type=fragment_type,
                                charge=charge,
                                aa_mass=aa_mass,
                            )
                            + (mod_mass + nl_mass) / charge,
                        )
                    )
    return sorted(fragments, key=operator.attrgetter("calc_mz"))


@nb.njit(cache=True)
def get_peak_annotation_indexes(
    spectrum_mz: np.ndarray,
    spectrum_intensity: np.ndarray,
    annotation_mz: np.ndarray,
    fragment_tol_mass: float,
    fragment_tol_mode: str,
    peak_assignment: str = "most_intense",
) -> List[Tuple[int, int]]:
    """
    Find matching indexes between observed m/z values and theoretical fragment
    labels.

    JIT helper function for `MsmsSpectrum.annotate_peaks`.

    Parameters
    ----------
    spectrum_mz : np.ndarray
        The m/z values of the spectrum fragment peaks.
    spectrum_intensity : np.ndarray
        The intensities of the spectrum fragment peaks.
    annotation_mz : np.ndarray
        The m/z values of the theoretical fragment labels.
    fragment_tol_mass : float
        Fragment mass tolerance to match spectrum peaks against theoretical
        peaks.
    fragment_tol_mode : {'Da', 'ppm'}
        Fragment mass tolerance unit. Either 'Da' or 'ppm'.
    peak_assignment : {'most_intense', 'nearest_mz'}, optional
        In case multiple peaks occur within the given mass window around a
        theoretical peak, only a single peak will be annotated with the
        fragment type:
        - 'most_intense': The most intense peak will be annotated (default).
        - 'nearest_mz':   The peak whose m/z is closest to the theoretical m/z
                          will be annotated.

    Returns
    -------
    List[Tuple[int, int]]
        A list of matching (peak index, annotation index) tuples.
    """
    annotation_i_map, peak_i_start = [], 0
    for fragment_i, fragment_mz in enumerate(annotation_mz):
        while (
            peak_i_start < len(spectrum_mz)
            and utils.mass_diff(
                spectrum_mz[peak_i_start],
                fragment_mz,
                fragment_tol_mode == "Da",
            )
            < -fragment_tol_mass
        ):
            peak_i_start += 1
        peak_i_stop = peak_i_start
        while (
            peak_i_stop < len(spectrum_mz)
            and utils.mass_diff(
                spectrum_mz[peak_i_stop],
                fragment_mz,
                fragment_tol_mode == "Da",
            )
            <= fragment_tol_mass
        ):
            peak_i_stop += 1
        if peak_i_start != peak_i_stop:
            peak_annotation_i = 0
            if peak_assignment == "nearest_mz":
                peak_annotation_i = np.argmin(
                    np.abs(spectrum_mz[peak_i_start:peak_i_stop] - fragment_mz)
                )
            elif peak_assignment == "most_intense":
                peak_annotation_i = np.argmax(
                    spectrum_intensity[peak_i_start:peak_i_stop]
                )
            annotation_i_map.append(
                (peak_i_start + peak_annotation_i, fragment_i)
            )
    return annotation_i_map
