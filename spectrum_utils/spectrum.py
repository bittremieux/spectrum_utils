import operator
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

import numba as nb
import numpy as np
import numpy.typing as npt
import pyteomics.usi

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
        charge: int = 1,
        adduct: Optional[str] = None,
        calc_mz: float = None,
    ) -> None:
        """
        Interpretation of a single fragment ion.

        This fragment annotation format is derived from the PSI peak
        interpretation specification:
        https://docs.google.com/document/d/1yEUNG4Ump6vnbMDs4iV4s3XISflmOkRAyqUuutcCG2w/edit?usp=sharing

        Parameters
        ----------
        ion_type : str
            Specifies the basic type of ion being described. Examples are
            b ions, y ions, immonium ions, unfragmented precursor ions,
            internal fragmentation ions, isobaric tag ions, etc.
        neutral_loss : str, optional
            A string of 0 to n loss components, described by their molecular
            formula. The default is no neutral loss.
        isotope : int, optional
            The isotope number above or below the monoisotope. The default is
            the monoisotopic peak (0).
        charge : int, optional
            The charge of the fragment. The default is charge 1.
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
            raise ValueError("Information specified for an unknown ion")
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
            ion = "?"
        else:
            ion = self.ion_type
            if self.neutral_loss is not None:
                ion += f"{self.neutral_loss}"
            if self.isotope != 0:
                ion += f"{self.isotope:+}i"
            if self.charge > 1:
                ion += f"^{self.charge}"
            if self.adduct is not None:
                ion += self.adduct
            return ion
        return f"FragmentAnnotation({ion}, mz={self.calc_mz})"

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
        if not isinstance(other, FragmentAnnotation):
            return False
        else:
            return str(self) == str(other)


def _get_theoretical_fragments(
    sequence: str,
    modifications: Optional[List[proforma.Modification]] = None,
    fragment_types: str = "by",
    max_charge: int = 1,
    neutral_losses: Optional[Dict[Optional[str], float]] = None,
) -> List[FragmentAnnotation]:
    """
    Get theoretical fragment annotations for the given sequence.

    Parameters
    ----------
    sequence : str
        The sequence for which the fragment annotations will be generated.
    modifications : Optional[List[proforma.Modification]], optional
        # FIXME
        Mapping of modification positions and mass differences. Valid positions
        are any amino acid index in the peptide (0-based), 'N-term', and
        'C-term'.
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
    if "B" in sequence:
        raise ValueError(
            "Explicitly specify aspartic acid (D) or asparagine (N) instead of"
            " the ambiguous B to compute the fragment annotations"
        )
    if "Z" in sequence:
        raise ValueError(
            "Explicitly specify glutamic acid (E) or glutamine (Q) instead of "
            "the ambiguous Z to compute the fragment annotations"
        )

    # FIXME
    # if modifications is not None:
    #     mods = modifications.copy()
    #     if 'N-term' in modifications:
    #         mods[-1] = mods['N-term']
    #         del mods['N-term']
    #     if 'C-term' in modifications:
    #         mods[len(sequence) + 1] = mods['C-term']
    #         del mods['C-term']
    # else:
    #     mods = {}

    neutral_losses = {None: 0} if neutral_losses is None else neutral_losses
    fragments = []
    # FIXME
    # # Single peak annotation.
    # if sequence == 'X':
    #     return [FragmentAnnotation('?', calc_mz=modifications[0])]
    # Get all possible peptide fragments.
    for i in range(1, len(sequence)):
        for fragment_type in fragment_types:
            # N-terminal fragment.
            if fragment_type in "abc":
                fragment_index = i
                fragment_sequence = sequence[:i]
                mod_mass = sum(
                    [mod.mass for mod in modifications if mod.position < i]
                )
            # C-terminal fragment.
            elif fragment_type in "xyz":
                fragment_index = len(sequence) - i
                fragment_sequence = sequence[i:]
                mod_mass = sum(
                    [mod.mass for mod in modifications if mod.position >= i]
                )
            else:
                raise ValueError(f"Unknown ion type: {fragment_type}")
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
def _init_spectrum(
    mz: np.ndarray, intensity: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    mz, intensity = mz.reshape(-1), intensity.reshape(-1)
    order = np.argsort(mz)
    return mz[order], intensity[order]


@nb.njit(cache=True)
def _round(
    mz: np.ndarray, intensity: np.ndarray, decimals: int, combine: str
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    JIT helper function for `MsmsSpectrum.round`.

    Parameters
    ----------
    mz : np.ndarray
        The mass-to-charge ratios of the spectrum fragment peaks.
    intensity : np.ndarray
        The intensities of the corresponding spectrum fragment peaks.
    decimals : int
        Number of decimal places to round the mass-to-charge ratios.
    combine : {'sum', 'max'}
        Method used to combine intensities from merged fragment peaks.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        A tuple consisting of the rounded mass-to-charge ratios and the
        corresponding intensities and peak annotation indexes.
    """
    mz_round = np.round_(mz, decimals, np.empty_like(mz, np.float32))
    mz_unique = np.unique(mz_round)
    # If peaks got merged by rounding the mass-to-charge ratios we need to
    # combine their intensities and annotations as well.
    if len(mz_unique) < len(mz_round):
        intensity_unique = np.zeros_like(mz_unique, np.float32)
        annotations_unique_idx = np.zeros_like(mz_unique, np.int64)
        combine_is_sum = combine == "sum"
        i_orig = 0
        offset = 0
        for i_unique in range(len(mz_unique)):
            # Check whether subsequent mz values got merged.
            while (
                abs(mz_unique[i_unique] - mz_round[i_orig + offset]) <= 1e-06
            ):
                offset += 1
            # Select the annotation of the most intense peak.
            annotations_unique_idx[i_unique] = i_orig + np.argmax(
                intensity[i_orig : i_orig + offset]
            )
            # Combine the corresponding intensities.
            intensity_unique[i_unique] = (
                intensity[i_orig : i_orig + offset].sum()
                if combine_is_sum
                else intensity[annotations_unique_idx[i_unique]]
            )

            i_orig += offset
            offset = 0

        return mz_unique, intensity_unique, annotations_unique_idx
    else:
        return mz_unique, intensity, np.arange(len(mz))


@nb.njit(cache=True)
def _get_mz_range_mask(
    mz: np.ndarray, min_mz: float, max_mz: float
) -> np.ndarray:
    """
    JIT helper function for `MsmsSpectrum.set_mz_range`.

    Parameters
    ----------
    mz : np.ndarray
        The mass-to-charge ratios of the spectrum fragment peaks.
    min_mz : float
        Minimum m/z (inclusive).
    max_mz : float
        Maximum m/z (inclusive).

    Returns
    -------
    np.ndarray
        Index mask specifying which peaks are inside of the given m/z range.
    """
    mask = np.full_like(mz, True, np.bool_)
    mz_i = 0
    while mz_i < len(mz) and mz[mz_i] < min_mz:
        mask[mz_i] = False
        mz_i += 1
    mz_i = len(mz) - 1
    while mz_i > 0 and mz[mz_i] > max_mz:
        mask[mz_i] = False
        mz_i -= 1
    return mask


@nb.njit(cache=True)
def _get_non_precursor_peak_mask(
    mz: np.ndarray,
    pep_mass: float,
    max_charge: int,
    isotope: int,
    fragment_tol_mass: float,
    fragment_tol_mode: str,
) -> np.ndarray:
    """
    JIT helper function for `MsmsSpectrum.remove_precursor_peak`.

    Parameters
    ----------
    mz : np.ndarray
        The mass-to-charge ratios of the spectrum fragment peaks.
    pep_mass : float
        The mono-isotopic mass of the uncharged peptide.
    max_charge : int
        The maximum precursor loss charge.
    isotope : int
        The number of isotopic peaks to be checked.
    fragment_tol_mass : float
            Fragment mass tolerance around the precursor mass to remove the
            precursor peak.
    fragment_tol_mode : {'Da', 'ppm'}
            Fragment mass tolerance unit. Either 'Da' or 'ppm'.

    Returns
    -------
    np.ndarray
        Index mask specifying which peaks are retained after precursor peak
        filtering.
    """
    remove_mz = []
    for charge in range(max_charge, 0, -1):
        for iso in range(isotope + 1):
            remove_mz.append((pep_mass + iso) / charge + 1.0072766)

    fragment_tol_mode_is_da = fragment_tol_mode == "Da"
    mask = np.full_like(mz, True, np.bool_)
    mz_i = remove_i = 0
    while mz_i < len(mz) and remove_i < len(remove_mz):
        md = utils.mass_diff(
            mz[mz_i], remove_mz[remove_i], fragment_tol_mode_is_da
        )
        if md < -fragment_tol_mass:
            mz_i += 1
        elif md > fragment_tol_mass:
            remove_i += 1
        else:
            mask[mz_i] = False
            mz_i += 1

    return mask


@nb.njit(cache=True)
def _get_filter_intensity_mask(
    intensity: np.ndarray, min_intensity: float, max_num_peaks: int
) -> np.ndarray:
    """
    JIT helper function for `MsmsSpectrum.filter_intensity`.

    Parameters
    ----------
    intensity : np.ndarray
        The intensities of the spectrum fragment peaks.
    min_intensity : float
        Remove peaks whose intensity is below `min_intensity` percentage of the
        intensity of the most intense peak.
    max_num_peaks : int
        Only retain the `max_num_peaks` most intense peaks.

    Returns
    -------
    np.ndarray
        Index mask specifying which peaks are retained after filtering the at
        most `max_num_peaks` most intense intensities above the minimum
        intensity threshold.
    """
    intensity_idx = np.argsort(intensity)
    min_intensity *= intensity[intensity_idx[-1]]
    # Discard low-intensity noise peaks.
    start_i = 0
    for start_i, intens in enumerate(intensity[intensity_idx]):
        if intens > min_intensity:
            break
    # Only retain at most the `max_num_peaks` most intense peaks.
    mask = np.full_like(intensity, False, np.bool_)
    mask[
        intensity_idx[max(start_i, len(intensity_idx) - max_num_peaks) :]
    ] = True
    return mask


@nb.njit(nb.float32[:](nb.float32[:], nb.int64), cache=True)
def _get_scaled_intensity_root(
    intensity: np.ndarray, degree: int
) -> np.ndarray:
    """
    JIT helper function for `MsmsSpectrum.scale_intensity`.

    Parameters
    ----------
    intensity : np.ndarray
        The intensities of the spectrum fragment peaks.
    degree : int
        The degree of the root scaling.

    Returns
    -------
    np.ndarray
        The root-scaled intensities.
    """
    return np.power(intensity, 1 / degree).astype(np.float32)


@nb.njit(nb.float32[:](nb.float32[:], nb.float64), cache=True)
def _get_scaled_intensity_log(intensity: np.ndarray, base: int) -> np.ndarray:
    """
    JIT helper function for `MsmsSpectrum.scale_intensity`.

    Parameters
    ----------
    intensity : np.ndarray
        The intensities of the spectrum fragment peaks.
    base : int
        The base of the log scaling.

    Returns
    -------
    np.ndarray
        The log-scaled intensities.
    """
    return (np.log1p(intensity) / np.log(base)).astype(np.float32)


@nb.njit(nb.float32[:](nb.float32[:], nb.int64), cache=True)
def _get_scaled_intensity_rank(
    intensity: np.ndarray, max_rank: int
) -> np.ndarray:
    """
    JIT helper function for `MsmsSpectrum.scale_intensity`.

    Parameters
    ----------
    intensity : np.ndarray
        The intensities of the spectrum fragment peaks.
    max_rank : int
        The maximum rank of the rank scaling.

    Returns
    -------
    np.ndarray
        The rank-scaled intensities.
    """
    return (max_rank - np.argsort(np.argsort(intensity)[::-1])).astype(
        np.float32
    )


@nb.njit(nb.float32[:](nb.float32[:], nb.float32), fastmath=True, cache=True)
def _scale_intensity_max(
    intensity: np.ndarray, max_intensity: float
) -> np.ndarray:
    """
    JIT helper function for `MsmsSpectrum.scale_intensity`.

    Parameters
    ----------
    intensity : np.ndarray
        The intensities of the spectrum fragment peaks.
    max_intensity : float
        Intensity of the most intense peak relative to which the peaks will be
        scaled.

    Returns
    -------
    np.ndarray
        The maximum-scaled intensities.
    """
    return intensity * max_intensity / intensity.max()


@nb.njit(cache=True)
def _get_fragment_annotation_map(
    spectrum_mz: np.ndarray,
    spectrum_intensity: np.ndarray,
    labels_mz: np.ndarray,
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
        The mass-to-charge values of the spectrum fragment peaks.
    spectrum_intensity : np.ndarray
        The intensities of the spectrum fragment peaks.
    labels_mz : np.ndarray
        The mass-to-charge values of the theoretical fragment labels.
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
        A list of (peak index, annotation index) tuples.
    """
    annotation_i_map = []
    peak_i_start = 0
    for fragment_i, fragment_mz in enumerate(labels_mz):
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


@nb.njit(cache=True)
def _get_mz_peak_index(
    spectrum_mz: np.ndarray,
    spectrum_intensity: np.ndarray,
    fragment_mz: float,
    fragment_tol_mass: float,
    fragment_tol_mode: str,
    peak_assignment: str = "most_intense",
) -> Optional[int]:
    """
    Find the best m/z peak in a spectrum relative to the given m/z value.

    Parameters
    ----------
    spectrum_mz : np.ndarray
        The mass-to-charge varlues of the spectrum fragment peaks.
    spectrum_intensity : np.ndarray
        The intensities of the spectrum fragment peaks.
    fragment_mz : float
        The m/z of the requested peak.
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
    int
        None if no peak was found within the given m/z window, else the index
        of the best m/z peak relative to the given m/z value.
    """
    # Binary search to find the best matching peak.
    md = (
        fragment_tol_mass
        if fragment_tol_mode == "Da"
        else fragment_tol_mass / 10 ** 6 * fragment_mz
    )
    peak_i_start, peak_i_stop = np.searchsorted(
        spectrum_mz, [fragment_mz - md, fragment_mz + md]
    )
    if peak_i_start == peak_i_stop:
        return None
    else:
        peak_annotation_i = 0
        if peak_assignment == "nearest_mz":
            peak_annotation_i = np.argmin(
                np.abs(spectrum_mz[peak_i_start:peak_i_stop] - fragment_mz)
            )
        elif peak_assignment == "most_intense":
            peak_annotation_i = np.argmax(
                spectrum_intensity[peak_i_start:peak_i_stop]
            )
        return peak_i_start + peak_annotation_i


class MsmsSpectrum:
    """
    Class representing a tandem mass spectrum.
    """

    def __init__(
        self,
        identifier: str,
        precursor_mz: float,
        precursor_charge: int,
        mz: Union[np.ndarray, Iterable],
        intensity: Union[np.ndarray, Iterable],
        retention_time: Optional[float] = None,
        mz_dtype: npt.DTypeLike = np.float64,
        intensity_dtype: npt.DTypeLike = np.float32,
    ) -> None:
        """
        Instantiate a new `MsmsSpectrum` consisting of fragment peaks.

        Parameters
        ----------
        identifier : str
            Spectrum identifier. It is recommended to use a unique and
            interpretable identifier, such as a `Universal Spectrum Identifier
            (USI) <https://psidev.info/usi>`_ as defined by the Proteomics
            Standards Initiative.
        precursor_mz : float
            Precursor ion mass-to-charge ratio.
        precursor_charge : int
            Precursor ion charge.
        mz : array_like
            Mass-to-charge ratios of the fragment peaks.
        intensity : array_like
            Intensities of the corresponding fragment peaks in `mz`.
        retention_time : Optional[float], optional
            Retention time at which the spectrum was acquired (the default is
            None, which indicates that retention time is unspecified/unknown).
        mz_dtype : npt.DTypeLike
            The dtype of the m/z values (the default is float64).
        intensity_dtype : npt.DTypeLike
            The dtype of the intensity values (the default is float32).
        """
        self.identifier = identifier
        self.precursor_mz = precursor_mz
        self.precursor_charge = precursor_charge
        if len(mz) != len(intensity):
            raise ValueError(
                "The mass-to-charge and intensity arrays should "
                "have equal lengths"
            )
        self._mz, self._intensity = _init_spectrum(
            np.require(mz, mz_dtype, "W"),
            np.require(intensity, intensity_dtype, "W"),
        )
        self.proforma, self._annotation = None, None
        self.retention_time = retention_time

    @classmethod
    def from_usi(
        cls, usi: str, backend: str = "aggregator", **kwargs
    ) -> "MsmsSpectrum":
        """
        Construct a spectrum from a public resource specified by its Universal
        Spectrum Identifier (USI).

        See the `Universal Spectrum Identifier (USI)
        <https://psidev.info/usi>`_ specification by the Proteomics Standards
        Initiative for more information on how to use USIs to refer to spectra.

        See the `Pyteomics documentation
        <https://pyteomics.readthedocs.io/en/latest/api/usi.html#pyteomics.usi.proxi>`_
        for more details on the available PROXI backends.

        Parameters
        ----------
        usi : str
            The USI from which to generate the spectrum.
        backend : str
            PROXI host backend (default: 'aggregator').
        kwargs
            Extra arguments to construct the spectrum that might be missing
            from the PROXI response (e.g. `precursor_mz` or `precursor_charge`)
            and the PROXI backend (see the `pyteomics.usi.proxi`
            documentation).

        Returns
        -------
        MsmsSpectrum
            The spectrum corresponding to the given USI

        Raises
        ------
        ValueError
            If the PROXI response does not contain information about the
            precursor m/z or precursor charge. Explicit precursor m/z and
            precursor charge values can be provided using the `precursor_mz`
            and `precursor_charge` keyword arguments respectively.
        """
        spectrum_dict = pyteomics.usi.proxi(usi, backend, **kwargs)
        for attr in spectrum_dict["attributes"]:
            if attr["accession"] in ("MS:1000827", "MS:1000744", "MS:1002234"):
                kwargs["precursor_mz"] = float(attr["value"])
                break
        else:
            if "precursor_mz" not in kwargs:
                raise ValueError(
                    "Unknown precursor m/z from USI. Specify the "
                    "precursor m/z directly."
                )
        for attr in spectrum_dict["attributes"]:
            if attr["accession"] == "MS:1000041":
                kwargs["precursor_charge"] = int(attr["value"])
                break
        else:
            if "precursor_charge" not in kwargs:
                raise ValueError(
                    "Unknown precursor charge from USI. Specify "
                    "the precursor charge directly."
                )
        mz = spectrum_dict["m/z array"]
        intensity = spectrum_dict["intensity array"]
        return cls(identifier=usi, mz=mz, intensity=intensity, **kwargs)

    @property
    def mz(self) -> np.ndarray:
        """
        The mass-to-charge ratios of the fragment peaks.

        Returns
        -------
        np.ndarray
            The mass-to-charge ratios of the fragment peaks.
        """
        return self._mz

    @property
    def intensity(self) -> np.ndarray:
        """
        TGet or set the intensities of the fragment peaks.

        When setting new intensity values it should be possible to convert the
        given values to a NumPy array and the number of intensity values should
        be equal to the number of m/z (and annotation) values.

        Returns
        -------
        np.ndarray
            The intensity values of the fragment peaks.
        """
        return self._intensity

    @property
    def annotation(self) -> Optional[np.ndarray]:
        """
        The annotated labels of the fragment peaks.

        Returns
        -------
        Optional[np.ndarray]
            The labels of the fragment peaks, or None if the spectrum has not
            been annotated.
        """
        return self._annotation

    def round(self, decimals: int = 0, combine: str = "sum") -> "MsmsSpectrum":
        """
        Round the mass-to-charge ratios of the fragment peaks to the given
        number of decimals.

        Peaks that have the same mass-to-charge ratio after rounding will be
        combined using the specified strategy.
        If multiple peaks are merged into a single peak it will be annotated
        with the annotation of the most intense peak.

        Parameters
        ----------
        decimals : int, optional
            Number of decimal places to round the `mz` to (default: 0).
            If decimals is negative, it specifies the number of positions to
            the left of the decimal point.
        combine : {'sum', 'max'}
            Method used to combine intensities from merged fragment peaks.
            ``sum`` specifies that the intensities of the merged fragment peaks
            will be summed, ``max`` specifies that the maximum intensity of
            the fragment peaks that are merged is used (the default is
            ``sum``).

        Returns
        -------
        MsmsSpectrum
        """
        self._mz, self._intensity, annotation_idx = _round(
            self._mz, self._intensity, decimals, combine
        )
        if self._annotation is not None:
            self._annotation = self._annotation[annotation_idx]
        return self

    def set_mz_range(
        self, min_mz: Optional[float] = None, max_mz: Optional[float] = None
    ) -> "MsmsSpectrum":
        """
        Restrict the mass-to-charge ratios of the fragment peaks to the
        given range.

        Parameters
        ----------
        min_mz : Optional[float], optional
            Minimum m/z (inclusive). If not set no minimal m/z restriction will
            occur.
        max_mz : Optional[float], optional
            Maximum m/z (inclusive). If not set no maximal m/z restriction will
            occur.

        Returns
        -------
        MsmsSpectrum
        """
        if min_mz is None and max_mz is None:
            return self
        else:
            if min_mz is None:
                min_mz = self._mz[0]
            if max_mz is None:
                max_mz = self._mz[-1]
            if max_mz < min_mz:
                min_mz, max_mz = max_mz, min_mz
        mz_range_mask = _get_mz_range_mask(self._mz, min_mz, max_mz)
        self._mz = self._mz[mz_range_mask]
        self._intensity = self._intensity[mz_range_mask]
        if self._annotation is not None:
            self._annotation = self._annotation[mz_range_mask]
        return self

    def remove_precursor_peak(
        self,
        fragment_tol_mass: float,
        fragment_tol_mode: str,
        isotope: int = 0,
    ) -> "MsmsSpectrum":
        """
        Remove fragment peak(s) close to the precursor mass-to-charge ratio.

        Parameters
        ----------
        fragment_tol_mass : float
            Fragment mass tolerance around the precursor mass to remove the
            precursor peak.
        fragment_tol_mode : {'Da', 'ppm'}
            Fragment mass tolerance unit. Either 'Da' or 'ppm'.
        isotope : int
            The number of precursor isotopic peaks to be checked (the default
            is 0 to check only the mono-isotopic peaks).

        Returns
        -------
        MsmsSpectrum
        """
        # FIXME: This assumes M+xH charged ions.
        neutral_mass = (self.precursor_mz - 1.0072766) * self.precursor_charge
        peak_mask = _get_non_precursor_peak_mask(
            self._mz,
            neutral_mass,
            self.precursor_charge,
            isotope,
            fragment_tol_mass,
            fragment_tol_mode,
        )
        self._mz = self.mz[peak_mask]
        self._intensity = self.intensity[peak_mask]
        if self._annotation is not None:
            self._annotation = self._annotation[peak_mask]
        return self

    def filter_intensity(
        self, min_intensity: float = 0.0, max_num_peaks: Optional[int] = None
    ) -> "MsmsSpectrum":
        """
        Remove low-intensity fragment peaks.

        Only the `max_num_peaks` most intense fragment peaks are retained.
        Additionally, noise peaks whose intensity is below `min_intensity`
        percentage of the intensity of the most intense peak are removed.

        Parameters
        ----------
        min_intensity : float, optional
            Remove peaks whose intensity is below `min_intensity` percentage
            of the intensity of the most intense peak (the default is 0, which
            means that no minimum intensity filter will be applied).
        max_num_peaks : Optional[int], optional
            Only retain the `max_num_peaks` most intense peaks (the default is
            None, which retains all peaks).

        Returns
        -------
        MsmsSpectrum
        """
        if max_num_peaks is None:
            max_num_peaks = len(self._intensity)
        intensity_mask = _get_filter_intensity_mask(
            self._intensity, min_intensity, max_num_peaks
        )
        self._mz = self._mz[intensity_mask]
        self._intensity = self._intensity[intensity_mask]
        if self._annotation is not None:
            self._annotation = self._annotation[intensity_mask]
        return self

    def scale_intensity(
        self,
        scaling: Optional[str] = None,
        max_intensity: Optional[float] = None,
        **kwargs,
    ) -> "MsmsSpectrum":
        """
        Scale the intensity of the fragment peaks.

        Two types of scaling can be performed: scaling all peaks using a
        specific transformation and scaling the peaks relative to the most
        intense peak.

        Parameters
        ----------
        scaling : {'root', 'log', 'rank'}, optional
            Method to scale the peak intensities (the default is None, which
            means that no transformation will be performed).
            Potential transformation options are:

            - 'root': Root-transform the peak intensities. The default is a
              square root transformation (`degree` is 2). The degree of the
              root can be specified using the `degree` kwarg.
            - 'log':  Log-transform the peak intensities. The default is a log2
              transformation (`base` is 2) after summing the intensities with 1
              to avoid negative values after the transformation. The base of
              the logarithm can be specified using the `base` kwarg.
            - 'rank': Rank-transform the peak intensities. The maximum rank of
              the most intense peak can be specified using the `max_rank`
              kwarg, by default the number of peaks in the spectrum is used as
              the maximum rank. Note that `max_rank` should be greater than or
              equal to the number of peaks in the spectrum.

        max_intensity : Optional[float], optional
            Intensity of the most intense peak relative to which the peaks will
            be scaled (the default is None, which means that no scaling
            relative to the most intense peak will be performed).

        Returns
        -------
        MsmsSpectrum
        """
        if scaling == "root":
            self._intensity = _get_scaled_intensity_root(
                self._intensity, kwargs.get("degree", 2)
            )
        elif scaling == "log":
            self._intensity = _get_scaled_intensity_log(
                self._intensity, kwargs.get("base", 2)
            )
        elif scaling == "rank":
            max_rank = kwargs.get("max_rank", len(self._intensity))
            if max_rank < len(self._intensity):
                raise ValueError(
                    "`max_rank` should be greater than or equal "
                    "to the number of peaks in the spectrum. See "
                    "`filter_intensity` to reduce the number of "
                    "peaks in the spectrum."
                )
            self._intensity = _get_scaled_intensity_rank(
                self._intensity, max_rank
            )
        if max_intensity is not None:
            self._intensity = _scale_intensity_max(
                self._intensity, max_intensity
            )
        return self

    def annotate_proforma(
        self,
        proforma_str: str,
        fragment_tol_mass: float,
        fragment_tol_mode: str,
        ion_types: str = "by",
        max_ion_charge: Optional[int] = None,
        peak_assignment: str = "most_intense",
        neutral_losses: Optional[Dict[Optional[str], float]] = None,
    ) -> "MsmsSpectrum":
        """
        Assign fragment ion labels to the peaks from a ProForma annotation.

        Parameters
        ----------
        proforma_str : str
            The ProForma spectrum annotation.
        fragment_tol_mass : float
            Fragment mass tolerance to match spectrum peaks against theoretical
            peaks.
        fragment_tol_mode : {'Da', 'ppm'}
            Fragment mass tolerance unit. Either 'Da' or 'ppm'.
        ion_types : str, optional
            Peptide fragment types to annotate. Can be any combination of 'a',
            'b', 'c', 'x', 'y', and 'z' (the default is 'by', which means that
            b-ions and y-ions will be annotated).
        max_ion_charge : Optional[int], optional
            All fragments up to and including the given charge will be
            annotated (by default all fragments with a charge up to the
            precursor minus, or minimum charge one, one will be annotated).
        peak_assignment : {'most_intense', 'nearest_mz'}, optional
            In case multiple peaks occur within the given mass window around a
            theoretical peak, only a single peak will be annotated with the
            fragment type:

            - 'most_intense': The most intense peak will be annotated
              (default).
            - 'nearest_mz': The peak whose m/z is closest to the theoretical
              m/z will be annotated.
        neutral_losses : Dict[Optional[str], float], optional
            Neutral losses to consider, specified as a dictionary with as keys
            the description (molecular formula) and as value the neutral loss.
            Note that the value should typically be negative.

        Returns
        -------
        MsmsSpectrum
        """
        self.proforma = proforma_str
        self._annotation = np.full_like(self.mz, None, object)
        # Be default, peak charges are assumed to be smaller than the precursor
        # charge.
        if max_ion_charge is None:
            max_ion_charge = max(1, self.precursor_charge - 1)
        # Make sure the standard peaks (without a neutral loss) are always
        # considered.
        if neutral_losses is not None and None not in neutral_losses:
            neutral_losses[None] = 0

        # Parse the ProForma string and find peaks that match the theoretical
        # fragments.
        for proteoform in proforma.parse(self.proforma):
            # TODO: Only localized modifications or all of them? Check.
            theoretical_fragments = _get_theoretical_fragments(
                proteoform.sequence,
                proteoform.modifications,
                ion_types,
                max_ion_charge,
                neutral_losses,
            )
            for annotation_i, fragment_i in _get_fragment_annotation_map(
                self.mz,
                self.intensity,
                np.asarray(
                    [fragment.calc_mz for fragment in theoretical_fragments]
                ),
                fragment_tol_mass,
                fragment_tol_mode,
                peak_assignment,
            ):
                # FIXME: Duplicate labels for the same peak are overwritten.
                self._annotation[annotation_i] = theoretical_fragments[fragment_i]

        return self
