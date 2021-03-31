import math
import operator
import warnings
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

import numba as nb
import numpy as np
try:
    from pyteomics import cmass as mass
except ImportError:
    from pyteomics import mass

from spectrum_utils import utils


_aa_mass = mass.std_aa_mass.copy()


def static_modification(amino_acid: str, mass_diff: float) -> None:
    """
    Globally modify the monoisotopic mass of an amino acid to set a static
    modification.

    Parameters
    ----------
    amino_acid : str
        The amino acid whose monoisotopic mass is modified.
    mass_diff : float
        The *mass difference* to be added to the amino acid's original
        monoisotopic mass.
    """
    global _aa_mass
    _aa_mass[amino_acid] += mass_diff


def reset_modifications() -> None:
    """
    Undo all static modifications and reset to the standard amino acid
    monoisotopic masses.
    """
    global _aa_mass
    _aa_mass = mass.std_aa_mass.copy()


class FragmentAnnotation:
    """
    Class representing a general fragment ion annotation.
    """

    def __init__(self, charge: int, calc_mz: float, annotation: str = None)\
            -> None:
        """
        Instantiate a new `FragmentAnnotation`.

        Parameters
        ----------
        charge : int
            The fragment ion charge if known, None otherwise.
        calc_mz : float
            The theoretical m/z value of the fragment.
        annotation : str
            The fragment's annotation string.
        """
        self.ion_type = 'unknown'
        self.charge = charge
        self.calc_mz = calc_mz
        self.annotation = annotation

    def _charge_to_str(self):
        """
        Convert a numeric charge to a string representation.
        """
        if self.charge is None:
            return 'unknown'
        elif self.charge > 0:
            return '+' * self.charge
        elif self.charge < 0:
            return '-' * -self.charge
        else:
            return 'undefined'

    def __repr__(self) -> str:
        return f"FragmentAnnotation(annotation='{self.annotation}', " \
               f"charge={self._charge_to_str()}, mz={self.calc_mz})"

    def __str__(self) -> str:
        return self.annotation

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, FragmentAnnotation):
            return False
        else:
            return self.annotation == other.annotation


class PeptideFragmentAnnotation(FragmentAnnotation):
    """
    Class representing a peptide fragment ion annotation.
    """

    def __init__(self, charge: int, calc_mz: float, ion_type: str,
                 ion_index: int) -> None:
        """
        Instantiate a new `PeptideFragmentAnnotation`.

        Parameters
        ----------
        charge : int
            The peptide fragment ion charge.
        calc_mz : float
            The theoretical m/z value of the peptide fragment.
        ion_type : {'a', 'b', 'c', 'x', 'y', 'z'}
            The peptide fragment ion type.
        ion_index : int
            The peptide fragment ion index.
        """
        super().__init__(charge, calc_mz)
        if ion_type not in 'abcxyz':
            raise ValueError(f'Unknown ion type: {ion_type}')
        self.ion_type = ion_type
        self.ion_index = ion_index
        self.annotation = f'{self.ion_type}{self.ion_index}' \
                          f'{self._charge_to_str()}'

    def __repr__(self) -> str:
        return f"PeptideFragmentAnnotation(ion_type='{self.ion_type}', " \
            f"ion_index={self.ion_index}, charge={self.charge}, " \
            f"mz={self.calc_mz})"

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, PeptideFragmentAnnotation):
            return False
        else:
            return (self.ion_type == other.ion_type and
                    self.ion_index == other.ion_index and
                    self.charge == other.charge and
                    math.isclose(self.calc_mz, other.calc_mz))


class MoleculeFragmentAnnotation(FragmentAnnotation):
    """
    Class representing a molecule fragment ion annotation.
    """

    def __init__(self, charge: int, calc_mz: float, smiles: str) -> None:
        """
        Instantiate a new `MoleculeFragmentAnnotation`.

        Parameters
        ----------
        charge : int
            The molecule fragment ion charge.
        calc_mz : float
            The theoretical m/z value of the molecule fragment.
        smiles : str
            The SMILES representation of the molecule.
        """
        super().__init__(charge, calc_mz)
        self.ion_type = 'molecule'
        self.smiles = smiles
        self.annotation = self.smiles

    def __repr__(self) -> str:
        return f"MoleculeFragmentAnnotation(smiles='{self.smiles}', " \
            f"charge={self.charge}, mz={self.calc_mz})"

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, MoleculeFragmentAnnotation):
            return False
        else:
            return (self.smiles == other.smiles and
                    self.charge == other.charge and
                    math.isclose(self.calc_mz, other.calc_mz))


def _get_theoretical_peptide_fragments(
        peptide: str,
        modifications: Optional[Dict[Union[float, str], float]] = None,
        types: str = 'by', max_charge: int = 1)\
        -> List[PeptideFragmentAnnotation]:
    """
    Get theoretical peptide fragments for the given peptide.

    Parameters
    ----------
    peptide : str
        The peptide sequence for which the fragments will be generated. The
        peptide sequence should only exist of the 20 standard amino acids.
    modifications : Optional[Dict[Union[float, str], float]], optional
        Mapping of modification positions and mass differences. Valid positions
        are any amino acid index in the peptide (0-based), 'N-term', and
        'C-term'.
    types : str, optional
        The fragment type. Can be any combination of 'a', 'b', 'c', 'x', 'y',
        and 'z' (the default is 'by', which means that b-ions and y-ions will
        be generated).
    max_charge : int, optional
        All fragments up to and including the given charge will be generated
        (the default is 1 to only generate singly-charged fragments).

    Returns
    -------
    List[Tuple[PeptideFragmentAnnotation, float]]
        A list of all fragments as (`PeptideFragmentAnnotation`, m/z) tuples
        sorted in ascending m/z order.
    """
    if modifications is not None:
        mods = modifications.copy()
        if 'N-term' in modifications:
            mods[-1] = mods['N-term']
            del mods['N-term']
        if 'C-term' in modifications:
            mods[len(peptide) + 1] = mods['C-term']
            del mods['C-term']
    else:
        mods = {}
    ions = []
    for i in range(1, len(peptide)):
        for ion_type in types:
            if ion_type in 'abc':   # N-terminal fragment.
                ion_index = i
                sequence = peptide[:i]
                mod_mass = sum([md for pos, md in mods.items() if pos < i])
            else:   # C-terminal fragment.
                ion_index = len(peptide) - i
                sequence = peptide[i:]
                mod_mass = sum([md for pos, md in mods.items() if pos >= i])
            for charge in range(1, max_charge + 1):
                ions.append(PeptideFragmentAnnotation(
                    charge, mass.fast_mass(
                        sequence=sequence, ion_type=ion_type,
                        charge=charge, aa_mass=_aa_mass) + mod_mass / charge,
                    ion_type, ion_index))
    return sorted(ions, key=operator.attrgetter('calc_mz'))


@nb.njit(nb.types.Tuple((nb.float32[:], nb.float32[:], nb.int64[:]))
         (nb.float32[::1], nb.float32[::1]))
def _init_spectrum(mz: np.ndarray, intensity: np.ndarray)\
        -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    mz, intensity = mz.reshape(-1), intensity.reshape(-1)
    order = np.argsort(mz)
    return mz[order], intensity[order], order


@nb.njit
def _round(mz: np.ndarray, intensity: np.ndarray, decimals: int, combine: str)\
        -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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
        combine_is_sum = combine == 'sum'
        i_orig = 0
        offset = 0
        for i_unique in range(len(mz_unique)):
            # Check whether subsequent mz values got merged.
            while (abs(mz_unique[i_unique] - mz_round[i_orig + offset])
                   <= 1e-06):
                offset += 1
            # Select the annotation of the most intense peak.
            annotations_unique_idx[i_unique] = i_orig + np.argmax(
                intensity[i_orig: i_orig + offset])
            # Combine the corresponding intensities.
            intensity_unique[i_unique] = (
                intensity[i_orig: i_orig + offset].sum() if combine_is_sum else
                intensity[annotations_unique_idx[i_unique]])

            i_orig += offset
            offset = 0

        return mz_unique, intensity_unique, annotations_unique_idx
    else:
        return mz_unique, intensity, np.arange(len(mz))


@nb.njit
def _get_mz_range_mask(mz: np.ndarray, min_mz: float, max_mz: float)\
        -> np.ndarray:
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


@nb.njit
def _get_non_precursor_peak_mask(mz: np.ndarray, pep_mass: float,
                                 max_charge: int, isotope: int,
                                 fragment_tol_mass: float,
                                 fragment_tol_mode: str)\
        -> np.ndarray:
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

    fragment_tol_mode_is_da = fragment_tol_mode == 'Da'
    mask = np.full_like(mz, True, np.bool_)
    mz_i = remove_i = 0
    while mz_i < len(mz) and remove_i < len(remove_mz):
        md = utils.mass_diff(mz[mz_i], remove_mz[remove_i],
                             fragment_tol_mode_is_da)
        if md < -fragment_tol_mass:
            mz_i += 1
        elif md > fragment_tol_mass:
            remove_i += 1
        else:
            mask[mz_i] = False
            mz_i += 1

    return mask


@nb.njit
def _get_filter_intensity_mask(intensity: np.ndarray, min_intensity: float,
                               max_num_peaks: int) -> np.ndarray:
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
    mask[intensity_idx[max(start_i, len(intensity_idx) - max_num_peaks):]] =\
        True
    return mask


@nb.njit(nb.float32[:](nb.float32[:], nb.int64))
def _get_scaled_intensity_root(intensity: np.ndarray, degree: int)\
        -> np.ndarray:
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


@nb.njit(nb.float32[:](nb.float32[:], nb.float64))
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


@nb.njit(nb.float32[:](nb.float32[:], nb.int64))
def _get_scaled_intensity_rank(intensity: np.ndarray, max_rank: int)\
        -> np.ndarray:
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
    return ((max_rank - np.argsort(np.argsort(intensity)[::-1]))
            .astype(np.float32))


@nb.njit(nb.float32[:](nb.float32[:], nb.float32), fastmath=True)
def _scale_intensity_max(intensity: np.ndarray, max_intensity: float)\
        -> np.ndarray:
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


@nb.njit
def _get_peptide_fragment_annotation_map(
        spectrum_mz: np.ndarray, spectrum_intensity: np.ndarray,
        annotation_mz: List[float], fragment_tol_mass: float,
        fragment_tol_mode: str, peak_assignment: str = 'most_intense')\
        -> List[Tuple[int, int]]:
    """
    JIT helper function for `MsmsSpectrum.annotate_peaks`.

    Parameters
    ----------
    spectrum_mz : np.ndarray
        The mass-to-charge varlues of the spectrum fragment peaks.
    spectrum_intensity : np.ndarray
        The intensities of the spectrum fragment peaks.
    annotation_mz : List[float]
        A list of mass-to-charge values of the peptide fragment annotations.
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
    A list of (peak index, annotation index) tuples.
    """
    annotation_i_map = []
    peak_i_start = 0
    for fragment_i, fragment_mz in enumerate(annotation_mz):
        while (peak_i_start < len(spectrum_mz) and
               utils.mass_diff(spectrum_mz[peak_i_start], fragment_mz,
                               fragment_tol_mode == 'Da')
               < -fragment_tol_mass):
            peak_i_start += 1
        peak_i_stop = peak_i_start
        while (peak_i_stop < len(spectrum_mz) and
               utils.mass_diff(spectrum_mz[peak_i_stop], fragment_mz,
                               fragment_tol_mode == 'Da')
               <= fragment_tol_mass):
            peak_i_stop += 1
        if peak_i_start != peak_i_stop:
            peak_annotation_i = 0
            if peak_assignment == 'nearest_mz':
                peak_annotation_i = np.argmin(np.abs(
                    spectrum_mz[peak_i_start: peak_i_stop] - fragment_mz))
            elif peak_assignment == 'most_intense':
                peak_annotation_i = np.argmax(
                    spectrum_intensity[peak_i_start: peak_i_stop])
            annotation_i_map.append((peak_i_start + peak_annotation_i,
                                     fragment_i))

    return annotation_i_map


@nb.njit
def _get_mz_peak_index(
        spectrum_mz: np.ndarray, spectrum_intensity: np.ndarray,
        fragment_mz: float, fragment_tol_mass: float,
        fragment_tol_mode: str, peak_assignment: str = 'most_intense')\
        -> Optional[int]:
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
    md = (fragment_tol_mass if fragment_tol_mode == 'Da' else
          fragment_tol_mass / 10**6 * fragment_mz)
    peak_i_start, peak_i_stop = np.searchsorted(
        spectrum_mz, [fragment_mz - md, fragment_mz + md])
    if peak_i_start == peak_i_stop:
        return None
    else:
        peak_annotation_i = 0
        if peak_assignment == 'nearest_mz':
            peak_annotation_i = np.argmin(np.abs(
                spectrum_mz[peak_i_start:peak_i_stop] - fragment_mz))
        elif peak_assignment == 'most_intense':
            peak_annotation_i = np.argmax(
                spectrum_intensity[peak_i_start:peak_i_stop])
        return peak_i_start + peak_annotation_i


class MsmsSpectrum:
    """
    Class representing a tandem mass spectrum.
    """

    def __init__(self,
                 identifier: str,
                 precursor_mz: float,
                 precursor_charge: int,
                 mz: Union[np.ndarray, Iterable],
                 intensity: Union[np.ndarray, Iterable],
                 annotation: Optional[Union[np.ndarray, Iterable]] = None,
                 retention_time: Optional[float] = None,
                 peptide: Optional[str] = None,
                 modifications: Optional[Dict[Union[int, str], float]] = None,
                 is_decoy: bool = False) -> None:
        """
        Instantiate a new `MsmsSpectrum` consisting of fragment peaks.

        Parameters
        ----------
        identifier : str
            (Unique) spectrum identifier. See for example the universal
            spectrum identifier specification by the Proteomics Standards
            Initiative.
        precursor_mz : float
            Precursor ion mass-to-charge ratio.
        precursor_charge : int
            Precursor ion charge.
        mz : array_like
            Mass-to-charge ratios of the fragment peaks.
        intensity : array_like
            Intensities of the corresponding fragment peaks in `mz`.
        annotation : Optional[array_like], optional
            Annotations of the corresponding fragment peaks in `mz` (the
            default is None, which indicates that the fragment peaks are not
            annotated).
        retention_time : Optional[float], optional
            Retention time at which the spectrum was acquired (the default is
            None, which indicates that retention time is unspecified/unknown).
        peptide : Optional[str], optional
            The peptide sequence corresponding to the spectrum (the default is
            None, which means that no peptide-spectrum match is specified). The
            peptide sequence should only exist of the 20 standard amino acids.
        modifications : Optional[Dict[Union[int, str], float]], optional
            Mapping of modification positions and mass differences. Valid
            positions are any amino acid index in the peptide (0-based),
            'N-term', and 'C-term'.
        is_decoy : bool, optional
            Flag indicating whether the `peptide` is a target or decoy
            peptide (the default is False, which implies a target peptide).
        """
        self.identifier = identifier

        self.precursor_mz = precursor_mz
        self.precursor_charge = precursor_charge

        if len(mz) != len(intensity):
            raise ValueError('The mass-to-charge and intensity arrays should '
                             'have equal length')

        self._mz, self._intensity, order = _init_spectrum(
            np.require(mz, np.float32, 'W'),
            np.require(intensity, np.float32, 'W'))

        if annotation is not None:
            self._annotation = np.asarray(annotation).reshape(-1)
            if len(self.mz) != len(self.annotation):
                raise ValueError('The mass-to-charge and annotation arrays '
                                 'should have equal length')
            else:
                self._annotation = self.annotation[order]
        else:
            self._annotation = None

        self.retention_time = retention_time
        if peptide is not None:
            self.peptide = peptide.upper()
            for aa in self.peptide:
                if aa not in mass.std_aa_mass:
                    raise ValueError(f'Unknown amino acid: {aa}')
        else:
            self.peptide = None
        if peptide is not None and modifications is not None:
            for mod_pos in modifications.keys():
                if mod_pos not in ('N-term', 'C-term'):
                    if not isinstance(mod_pos, int):
                        raise ValueError(f'Unknown modification position: '
                                         f'{mod_pos}')
                    elif mod_pos < 0 or mod_pos > len(peptide):
                        raise ValueError(f'Modification position exceeds '
                                         f'peptide bounds: {mod_pos}')
        else:
            modifications = None
        self.modifications = modifications
        self.is_decoy = is_decoy

    @property
    def mz(self) -> np.ndarray:
        """
        Get or set the mass-to-charge ratios of the fragment peaks.

        When setting new m/z values it should be possible to convert the given
        values to a NumPy array and the number of m/z values should be equal to
        the number of intensity (and annotation) values.

        Returns
        -------
        np.ndarray
            The mass-to-charge ratios of the fragment peaks.
        """
        return self._mz

    @mz.setter
    def mz(self, mz: Union[np.ndarray, Iterable]) -> None:
        if np.asarray(mz).ndim > 0:
            self._mz = np.asarray(mz)
        else:
            raise ValueError('Invalid m/z values')

    @property
    def intensity(self) -> np.ndarray:
        """
        Get or set the intensities of the fragment peaks.

        When setting new intensity values it should be possible to convert the
        given values to a NumPy array and the number of intensity values should
        be equal to the number of m/z (and annotation) values.

        Returns
        -------
        np.ndarray
            The intensity values of the fragment peaks.
        """
        return self._intensity

    @intensity.setter
    def intensity(self, intensity: Union[np.ndarray, Iterable]) -> None:
        if np.asarray(intensity).ndim > 0:
            self._intensity = np.asarray(intensity)
        else:
            raise ValueError('Invalid intensity values')

    @property
    def annotation(self) -> Optional[np.ndarray]:
        """
        Get or set the annotations of the fragment peaks.

        When setting new annotations it should be possible to convert the given
        values to a NumPy array (or None) and the number of annotations should
        be equal to the number of m/z and intensity values.

        Returns
        -------
        Optional[np.ndarray]
            The annotations of the fragment peaks or None if no annotations
            have been specified.
        """
        return self._annotation

    @annotation.setter
    def annotation(
            self, annotation: Optional[Union[np.ndarray, Iterable]] = None)\
            -> None:
        if annotation is None:
            self._annotation = None
        elif np.asarray(annotation).ndim > 0:
            self._annotation = np.asarray(annotation)
        else:
            raise ValueError('Invalid annotation values')

    def round(self, decimals: int = 0, combine: str = 'sum') -> 'MsmsSpectrum':
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
        self : `MsmsSpectrum`
        """
        self.mz, self.intensity, annotation_idx = _round(
            self.mz, self.intensity, decimals, combine)
        if self.annotation is not None:
            self.annotation = self.annotation[annotation_idx]

        return self

    def set_mz_range(self, min_mz: Optional[float] = None,
                     max_mz: Optional[float] = None) -> 'MsmsSpectrum':
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
        self : `MsmsSpectrum`
        """
        if min_mz is None and max_mz is None:
            return self
        else:
            if min_mz is None:
                min_mz = self.mz[0]
            if max_mz is None:
                max_mz = self.mz[-1]

        mz_range_mask = _get_mz_range_mask(self.mz, min_mz, max_mz)
        self.mz = self.mz[mz_range_mask]
        self.intensity = self.intensity[mz_range_mask]
        if self.annotation is not None:
            self.annotation = self.annotation[mz_range_mask]

        return self

    def remove_precursor_peak(self, fragment_tol_mass: float,
                              fragment_tol_mode: str, isotope: int = 0)\
            -> 'MsmsSpectrum':
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
        self : `MsmsSpectrum`
        """
        pep_mass = (self.precursor_mz - 1.0072766) * self.precursor_charge
        peak_mask = _get_non_precursor_peak_mask(
            self.mz, pep_mass, self.precursor_charge, isotope,
            fragment_tol_mass, fragment_tol_mode)
        self.mz = self.mz[peak_mask]
        self.intensity = self.intensity[peak_mask]
        if self.annotation is not None:
            self.annotation = self.annotation[peak_mask]

        return self

    def filter_intensity(self, min_intensity: float = 0.0,
                         max_num_peaks: Optional[int] = None)\
            -> 'MsmsSpectrum':
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
        self : `MsmsSpectrum`
        """
        if max_num_peaks is None:
            max_num_peaks = len(self.intensity)
        intensity_mask = _get_filter_intensity_mask(
            self.intensity, min_intensity, max_num_peaks)
        self.mz = self.mz[intensity_mask]
        self.intensity = self.intensity[intensity_mask]
        if self.annotation is not None:
            self.annotation = self.annotation[intensity_mask]

        return self

    def scale_intensity(self, scaling: Optional[str] = None,
                        max_intensity: Optional[float] = None,
                        **kwargs) -> 'MsmsSpectrum':
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
        self : `MsmsSpectrum`
        """
        if scaling == 'root':
            self.intensity = _get_scaled_intensity_root(
                self.intensity, kwargs.get('degree', 2))
        elif scaling == 'log':
            self.intensity = _get_scaled_intensity_log(
                self.intensity, kwargs.get('base', 2))
        elif scaling == 'rank':
            max_rank = kwargs.get('max_rank', len(self.intensity))
            if max_rank < len(self.intensity):
                raise ValueError('`max_rank` should be greater than or equal '
                                 'to the number of peaks in the spectrum. See '
                                 '`filter_intensity` to reduce the number of '
                                 'peaks in the spectrum.')
            self.intensity = _get_scaled_intensity_rank(
                self.intensity, max_rank)

        if max_intensity is not None:
            self.intensity = _scale_intensity_max(
                self.intensity, max_intensity)

        return self

    def annotate_peaks(self, *args, **kwargs):
        raise DeprecationWarning('Renamed to annotate_peptide_fragments')

    def annotate_peptide_fragments(self, fragment_tol_mass: float,
                                   fragment_tol_mode: str,
                                   ion_types: str = 'by',
                                   max_ion_charge: Optional[int] = None,
                                   peak_assignment: str = 'most_intense')\
            -> 'MsmsSpectrum':
        """
        Annotate peaks with their corresponding peptide fragment ion
        annotations.

        `self.annotation` will be overwritten and include
        `PeptideFragmentAnnotation` objects for matching peaks.

        Parameters
        ----------
        fragment_tol_mass : float
            Fragment mass tolerance to match spectrum peaks against theoretical
            peaks.
        fragment_tol_mode : {'Da', 'ppm'}
            Fragment mass tolerance unit. Either 'Da' or 'ppm'.
        ion_types : str, optional
            Fragment type to annotate. Can be any combination of 'a', 'b', 'c',
            'x', 'y', and 'z' (the default is 'by', which means that b-ions and
            y-ions will be annotated).
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

        Returns
        -------
        self : `MsmsSpectrum`
        """
        if self.peptide is None:
            raise ValueError('No peptide sequence available for the spectrum')
        if max_ion_charge is None:
            max_ion_charge = max(1, self.precursor_charge - 1)

        theoretical_fragments = _get_theoretical_peptide_fragments(
            self.peptide, self.modifications, ion_types, max_ion_charge)
        self.annotation = np.full_like(self.mz, None, object)
        with warnings.catch_warnings():
            # FIXME: Deprecated reflected list in Numba should be resolved from
            #        version 0.46.0 onwards.
            #  https://numba.pydata.org/numba-doc/latest/reference/deprecation.html#deprecation-of-reflection-for-list-and-set-types
            warnings.simplefilter('ignore', nb.NumbaPendingDeprecationWarning)
            for annotation_i, fragment_i in\
                    _get_peptide_fragment_annotation_map(
                        self.mz, self.intensity,
                        [fragment.calc_mz for fragment
                         in theoretical_fragments],
                        fragment_tol_mass, fragment_tol_mode, peak_assignment):
                self.annotation[annotation_i] =\
                    theoretical_fragments[fragment_i]

        return self

    def annotate_molecule_fragment(self, smiles: str, fragment_mz: float,
                                   fragment_charge: int,
                                   fragment_tol_mass: float,
                                   fragment_tol_mode: str,
                                   peak_assignment: str = 'most_intense')\
            -> 'MsmsSpectrum':
        """
        Annotate a peak (if present) with its corresponding molecule.

        The matching position in `self.annotation` will be overwritten by the
        provided molecule.

        Parameters
        ----------
        smiles : str
            The fragment molecule that will be annotated in SMILES format. Note
            that the SMILES string is not tested for validity.
        fragment_mz : float
            The expected m/z of the molecule.
        fragment_charge : int
            The charge of the molecule.
        fragment_tol_mass : float
            Fragment mass tolerance to match spectrum peaks against the
            theoretical molecule mass.
        fragment_tol_mode : {'Da', 'ppm'}
            Fragment mass tolerance unit. Either 'Da' or 'ppm'.
        peak_assignment : {'most_intense', 'nearest_mz'}, optional
            In case multiple peaks occur within the given mass window around
            the molecule's given mass, only a single peak will be annotated:

            - 'most_intense': The most intense peak will be annotated
              (default).
            - 'nearest_mz': The peak whose m/z is closest to the theoretical
              m/z will be annotated.

        Returns
        -------
        self : `MsmsSpectrum`
        """
        # Find the matching m/z value.
        peak_index = _get_mz_peak_index(self.mz, self.intensity, fragment_mz,
                                        fragment_tol_mass, fragment_tol_mode,
                                        peak_assignment)
        if peak_index is None:
            raise ValueError(f'No matching peak found for molecule "{smiles}"')
        else:
            # Initialize the annotations if they don't exist yet.
            if self.annotation is None:
                self.annotation = np.full_like(self.mz, None, object)
            # Set the molecule's annotation.
            self.annotation[peak_index] =\
                MoleculeFragmentAnnotation(fragment_charge, fragment_mz,
                                           smiles)

        return self

    def annotate_mz_fragment(self, fragment_mz: float, fragment_charge: int,
                             fragment_tol_mass: float, fragment_tol_mode: str,
                             peak_assignment: str = 'most_intense',
                             text: Optional[str] = None) -> 'MsmsSpectrum':
        """
        Annotate a peak (if present) with its m/z value or a custom provided
        string.

        The matching position in `self.annotation` will be overwritten.

        Parameters
        ----------
        fragment_mz : float
            The expected m/z to annotate.
        fragment_charge : int
            The peak charge.
        fragment_tol_mass : float
            Fragment mass tolerance to match spectrum peaks against the given
            m/z.
        fragment_tol_mode : {'Da', 'ppm'}
            Fragment mass tolerance unit. Either 'Da' or 'ppm'.
        peak_assignment : {'most_intense', 'nearest_mz'}, optional
            In case multiple peaks occur within the given mass window around
            the given m/z, only a single peak will be annotated:

            - 'most_intense': The most intense peak will be annotated
              (default).
            - 'nearest_mz': The peak whose m/z is closest to the given m/z will
              be annotated.
        text : Optional[str], optional
            The text to annotate the peak with. If None, its m/z value will be
            used.

        Returns
        -------
        self : `MsmsSpectrum`
        """
        # Find the matching m/z value.
        peak_index = _get_mz_peak_index(self.mz, self.intensity, fragment_mz,
                                        fragment_tol_mass, fragment_tol_mode,
                                        peak_assignment)
        if peak_index is None:
            raise ValueError(f'No matching peak found for {fragment_mz} m/z')
        else:
            # Initialize the annotations if they don't exist yet.
            if self.annotation is None:
                self.annotation = np.full_like(self.mz, None, object)
            # Set the peak annotation.
            self.annotation[peak_index] =\
                FragmentAnnotation(
                    fragment_charge, fragment_mz,
                    text if text is not None else str(fragment_mz))

        return self
