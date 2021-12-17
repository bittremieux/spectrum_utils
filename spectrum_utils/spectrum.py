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
    "J": 113.08406,  # leucine / isoleucine
    # "U": 150.95363,     # selenocysteine (in Pyteomics)
    # "O": 237.14772,     # pyrrolysine (in Pyteomics)
    "X": 0,  # any amino acid, gaps (zero mass)
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


def _get_theoretical_fragments(
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
def _get_peak_annotation_indexes(
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


class GnpsBackend(pyteomics.usi._PROXIBackend):

    _url_template = ("https://metabolomics-usi.ucsd.edu/proxi/v{version}/"
                     "spectra?usi1={usi}")

    def __init__(self, **kwargs):
        super(GnpsBackend, self).__init__("GNPS", self._url_template, **kwargs)


pyteomics.usi._proxies["GNPS"] = GnpsBackend
pyteomics.usi.AGGREGATOR = pyteomics.usi.PROXIAggregator()


@nb.experimental.jitclass
class MsmsSpectrumJit:
    """
    Helper MsmsSpectrum class for JIT-compiled fast data processing.
    """

    identifier: nb.types.string
    precursor_mz: nb.float64
    precursor_charge: nb.int8
    _mz: nb.float64[:]
    _intensity: nb.float32[:]
    retention_time: nb.float32

    def __init__(
        self,
        identifier: str,
        precursor_mz: float,
        precursor_charge: int,
        mz: np.ndarray,
        intensity: np.ndarray,
        retention_time: float,
    ) -> None:
        self.identifier = identifier
        self.precursor_mz = precursor_mz
        self.precursor_charge = precursor_charge
        if len(mz) != len(intensity):
            raise ValueError(
                "The m/z and intensity arrays should have equal lengths"
            )
        mz = np.asarray(mz, np.float64).reshape(-1)
        # Make sure the peaks are sorted by m/z.
        intensity = np.asarray(intensity, np.float32).reshape(-1)
        order = np.argsort(mz)
        self._mz, self._intensity = mz[order], intensity[order]
        self.retention_time = retention_time

    @property
    def mz(self) -> np.ndarray:
        return self._mz

    @property
    def intensity(self) -> np.ndarray:
        return self._intensity

    def round(
        self, decimals: int = 0, combine: str = "sum"
    ) -> "MsmsSpectrumJit":
        mz_round = np.round_(self._mz, decimals, np.empty_like(self._mz))
        mz_unique = np.unique(mz_round)
        if len(mz_unique) == len(mz_round):
            self._mz = mz_unique
        # If peaks got merged by rounding the mass-to-charge ratios we need to
        # combine their intensities and annotations as well.
        else:
            intensity_unique = np.zeros_like(mz_unique, np.float32)
            i_orig = offset = 0
            for i_unique in range(len(mz_unique)):
                # Check whether subsequent mz values got merged.
                while (
                    i_orig + offset < len(mz_round)
                    and abs(mz_unique[i_unique] - mz_round[i_orig + offset])
                    <= 1e-06
                ):
                    offset += 1
                # Combine the corresponding intensities.
                intensity_unique[i_unique] = (
                    self._intensity[i_orig : i_orig + offset].sum()
                    if combine == "sum"
                    else self._intensity[i_orig : i_orig + offset].max()
                )
                i_orig += offset
                offset = 0
            self._mz, self._intensity = mz_unique, intensity_unique
        return self

    def set_mz_range(
        self, min_mz: Optional[float] = None, max_mz: Optional[float] = None
    ) -> "MsmsSpectrumJit":
        if min_mz is None and max_mz is None:
            return self
        else:
            if min_mz is None:
                min_mz = self._mz[0]
            if max_mz is None:
                max_mz = self._mz[-1]
            if max_mz < min_mz:
                min_mz, max_mz = max_mz, min_mz
        min_i, max_i = 0, len(self._mz)
        while min_i < len(self._mz) and self._mz[min_i] < min_mz:
            min_i += 1
        while max_i > 0 and self._mz[max_i - 1] > max_mz:
            max_i -= 1
        self._mz = self._mz[min_i:max_i]
        self._intensity = self._intensity[min_i:max_i]
        return self

    def remove_precursor_peak(
        self,
        fragment_tol_mass: float,
        fragment_tol_mode: str,
        isotope: int = 0,
    ) -> "MsmsSpectrumJit":
        # TODO: This assumes [M+H]x charged ions.
        neutral_mass = (self.precursor_mz - 1.0072766) * self.precursor_charge
        remove_mz = [
            (neutral_mass + iso) / charge + 1.0072766
            for charge in range(self.precursor_charge, 0, -1)
            for iso in range(isotope + 1)
        ]
        mask = np.full_like(self._mz, True, np.bool_)
        mz_i = remove_i = 0
        while mz_i < len(self._mz) and remove_i < len(remove_mz):
            md = utils.mass_diff(
                self._mz[mz_i], remove_mz[remove_i], fragment_tol_mode == "Da"
            )
            if md < -fragment_tol_mass:
                mz_i += 1
            elif md > fragment_tol_mass:
                remove_i += 1
            else:
                mask[mz_i] = False
                mz_i += 1
        self._mz, self._intensity = self._mz[mask], self._intensity[mask]
        return self

    def filter_intensity(
        self, min_intensity: float = 0.0, max_num_peaks: Optional[int] = None
    ) -> "MsmsSpectrumJit":
        if max_num_peaks is None:
            max_num_peaks = len(self._intensity)

        intensity_idx = np.argsort(self._intensity)
        min_intensity *= self._intensity[intensity_idx[-1]]
        # Discard low-intensity noise peaks.
        start_i = 0
        for intens in self._intensity[intensity_idx]:
            if intens > min_intensity:
                break
            start_i += 1
        # Only retain at most the `max_num_peaks` most intense peaks.
        mask = np.full_like(self._intensity, False, np.bool_)
        mask[
            intensity_idx[max(start_i, len(intensity_idx) - max_num_peaks) :]
        ] = True
        self._mz, self._intensity = self._mz[mask], self._intensity[mask]
        return self

    def scale_intensity(
        self,
        scaling: Optional[str] = None,
        max_intensity: Optional[float] = None,
        degree: int = 2,
        base: int = 2,
        max_rank: Optional[int] = None,
    ) -> "MsmsSpectrumJit":
        if scaling == "root":
            self._intensity = np.power(
                self._intensity, 1 / degree
            ).astype(np.float32)
        elif scaling == "log":
            self._intensity = (
                np.log1p(self._intensity) / np.log(base)
            ).astype(np.float32)
        elif scaling == "rank":
            if max_rank is None:
                max_rank = len(self._intensity)
            if max_rank < len(self._intensity):
                raise ValueError(
                    "`max_rank` should be greater than or equal to the number "
                    "of peaks in the spectrum. See `filter_intensity` to "
                    "reduce the number of peaks in the spectrum."
                )
            self._intensity = (
                max_rank - np.argsort(np.argsort(self._intensity)[::-1])
            ).astype(np.float32)
        if max_intensity is not None:
            self._intensity = (
                self._intensity * max_intensity / self._intensity.max()
            ).astype(np.float32)
        return self


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
        retention_time: float = np.nan,
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
            Precursor ion m/z.
        precursor_charge : int
            Precursor ion charge.
        mz : array_like
            M/z values of the fragment peaks.
        intensity : array_like
            Intensities of the corresponding fragment peaks in `mz`.
        retention_time : float, optional
            Retention time at which the spectrum was acquired (the default is
            np.nan, which indicates that retention time is
            unspecified/unknown).
        """
        self._inner = MsmsSpectrumJit(
            identifier,
            precursor_mz,
            precursor_charge,
            np.require(mz, np.float64, "W"),
            np.require(intensity, np.float32, "W"),
            retention_time if retention_time is not None else np.nan,
        )
        self.proforma, self._annotation = None, None

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
    def identifier(self) -> str:
        return self._inner.identifier

    @identifier.setter
    def identifier(self, identifier: str):
        self._inner.identifier = identifier

    @property
    def precursor_mz(self) -> float:
        return self._inner.precursor_mz

    @precursor_mz.setter
    def precursor_mz(self, precursor_mz: float):
        self._inner.precursor_mz = precursor_mz

    @property
    def precursor_charge(self) -> int:
        return self._inner.precursor_charge

    @precursor_charge.setter
    def precursor_charge(self, precursor_charge: int):
        self._inner.precursor_charge = precursor_charge

    @property
    def retention_time(self) -> float:
        return self._inner.retention_time

    @retention_time.setter
    def retention_time(self, retention_time: float):
        self._inner.retention_time = retention_time

    @property
    def mz(self) -> np.ndarray:
        """
        The m/z values of the fragment peaks.

        Returns
        -------
        np.ndarray
            The m/z values of the fragment peaks.
        """
        return self._inner.mz

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
        return self._inner.intensity

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
        Round the m/z values of the fragment peaks to the given number of
        decimals.

        Intensities of peaks that have the same m/z value after rounding will
        be combined using the specified strategy.

        Note: This peak manipulation will reset any ProForma annotations that
        were applied previously.

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
        self.proforma, self._annotation = None, None
        self._inner.round(decimals, combine)
        return self

    def set_mz_range(
        self, min_mz: Optional[float] = None, max_mz: Optional[float] = None
    ) -> "MsmsSpectrum":
        """
        Restrict the m/z values of the fragment peaks to the given range.

        Note: This peak manipulation will reset any ProForma annotations that
        were applied previously.

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
        self.proforma, self._annotation = None, None
        self._inner.set_mz_range(min_mz, max_mz)
        return self

    def remove_precursor_peak(
        self,
        fragment_tol_mass: float,
        fragment_tol_mode: str,
        isotope: int = 0,
    ) -> "MsmsSpectrum":
        """
        Remove fragment peak(s) close to the precursor m/z.

        Note: This peak manipulation will reset any ProForma annotations that
        were applied previously.

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
        self.proforma, self._annotation = None, None
        self._inner.remove_precursor_peak(
            fragment_tol_mass, fragment_tol_mode, isotope
        )
        return self

    def filter_intensity(
        self, min_intensity: float = 0.0, max_num_peaks: Optional[int] = None
    ) -> "MsmsSpectrum":
        """
        Remove low-intensity fragment peaks.

        Only the `max_num_peaks` most intense fragment peaks are retained.
        Additionally, noise peaks whose intensity is below `min_intensity`
        percentage of the intensity of the most intense peak are removed.

        Note: This peak manipulation will reset any ProForma annotations that
        were applied previously.

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
        self.proforma, self._annotation = None, None
        self._inner.filter_intensity(min_intensity, max_num_peaks)
        return self

    def scale_intensity(
        self,
        scaling: Optional[str] = None,
        max_intensity: Optional[float] = None,
        *,
        degree: int = 2,
        base: int = 2,
        max_rank: Optional[int] = None,
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
              the maximum rank.

        degree: int
            The degree of the root transformation (the default is 2 for square
            root scaling).
        base: int
            The base of the log transformation (the default is 2 for log2
            scaling).
        max_rank: Optional[int]
            The maximum rank of the rank transformation (the default is None,
            which uses the number of peaks in the spectrum as maximum rank).
            `max_rank` should be greater than or equal to the number of peaks
            in the spectrum.
        max_intensity : Optional[float], optional
            Intensity of the most intense peak relative to which the peaks will
            be scaled (the default is None, which means that no scaling
            relative to the most intense peak will be performed).

        Returns
        -------
        MsmsSpectrum
        """
        self._inner.scale_intensity(
            scaling, max_intensity, degree, base, max_rank
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
            fragments = _get_theoretical_fragments(
                proteoform,
                ion_types,
                max_ion_charge,
                neutral_losses,
            )
            for annotation_i, fragment_i in _get_peak_annotation_indexes(
                self.mz,
                self.intensity,
                np.asarray([f.calc_mz for f in fragments]),
                fragment_tol_mass,
                fragment_tol_mode,
                peak_assignment,
            ):
                # TODO: Support multiple annotations per peak. Currently
                #       duplicate labels for the same peak are overwritten.
                self._annotation[annotation_i] = fragments[fragment_i]

        return self
