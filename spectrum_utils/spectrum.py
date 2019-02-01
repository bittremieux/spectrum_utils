import operator
from typing import Any
from typing import Iterable
from typing import List
from typing import Tuple
from typing import Union

import numba as nb
import numpy as np
from pyteomics import mass
from pyteomics import parser

from spectrum_utils import utils


class FragmentAnnotation:
    """
    Class representing a fragment io annotation.
    """

    def __init__(self, ion_type: str, ion_index: int, charge: int) -> None:
        """
        Instantiate a new `FragmentAnnotation`.

        Parameters
        ----------
        ion_type : {'a', 'b', 'c', 'x', 'y', 'z'}
            The fragment ion type.
        ion_index : int
            The fragment ion index.
        charge : int
            The fragment ion charge.
        """
        self.ion_type = ion_type
        self.ion_index = ion_index
        self.charge = charge

    def __repr__(self) -> str:
        return f"FragmentAnnotation(ion_type='{self.ion_type}', " \
            f"ion_index={self.ion_index}, charge={self.charge})"

    def __str__(self) -> str:
        return f'{self.ion_type}{self.ion_index}{"+" * self.charge}'

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, FragmentAnnotation):
            return False
        else:
            return (self.ion_type == other.ion_type and
                    self.ion_index == other.ion_index and
                    self.charge == other.charge)


def _get_theoretical_peptide_fragments(peptide: str, types: str = 'by',
                                       max_charge: int = 1)\
        -> List[Tuple[FragmentAnnotation, float]]:
    """
    Get theoretical fragments for the given peptide.

    Parameters
    ----------
    peptide : str
        The peptide sequence for which the fragments will be generated.
    types : str, optional
        The fragment type. Can be any combination of 'a', 'b', 'c', 'x', 'y',
        and 'z' (the default is 'by', which means that b-ions and y-ions will
        be generated).
    max_charge : int, optional
        All fragments up to and including the given charge will be generated
        (the default is 1 to only generate singly-charged fragments).

    Returns
    -------
        A list of all fragments as (`FragmentAnnotation`, m/z) tuples sorted in
        ascending m/z order.
    """
    ions = []
    amino_acids = parser.parse(peptide)
    for i in range(1, len(amino_acids)):
        for ion_type in types:
            for charge in range(1, max_charge + 1):
                if ion_type in 'abc':
                    ions.append((
                        FragmentAnnotation(ion_type, i, charge),
                        mass.calculate_mass(sequence=''.join(amino_acids[:i]),
                                            ion_type=ion_type,
                                            charge=charge)))
                else:
                    ions.append((
                        FragmentAnnotation(ion_type, len(peptide) - i, charge),
                        mass.calculate_mass(sequence=''.join(amino_acids[i:]),
                                            ion_type=ion_type,
                                            charge=charge)))
    return sorted(ions, key=operator.itemgetter(1))


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
        return mz, intensity, np.arange(len(mz))


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
    remove_mz = np.array([(pep_mass + iso) / charge + 1.0072766
                          for charge in range(max_charge, 0, -1)
                          for iso in range(isotope + 1)], np.float32)

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


@nb.njit
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
    return np.power(intensity, 1 / degree)


@nb.njit
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
    return np.log1p(intensity) / np.log(base)


@nb.njit
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


@nb.njit(fastmath=True)
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
                 annotation: Union[np.ndarray, Iterable] = None,
                 retention_time: float = None,
                 peptide: str = None,
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
        annotation : array_like, optional
            Annotations of the corresponding fragment peaks in `mz` (the
            default is None, which indicates that the fragment peaks are not
            annotated).
        retention_time : float, optional
            Retention time at which the spectrum was acquired (the default is
            None, which indicates that retention time is unspecified/unknown).
        peptide : str, optional
            The peptide sequence corresponding to the spectrum (the default is
            None, which means that no peptide-spectrum match is specified).
        is_decoy : bool, optional
            Flag indicating whether the `peptide` is a target or decoy
            peptide (the default is False, which implies a target peptide).
        """
        self.identifier = identifier

        self.precursor_mz = precursor_mz
        self.precursor_charge = precursor_charge

        self.mz = np.asarray(mz, np.float32).reshape(-1)
        self.intensity = np.asarray(intensity, np.float32).reshape(-1)
        if len(self.mz) != len(self.intensity):
            raise ValueError('The mass-to-charge and intensity arrays should '
                             'have equal length')
        self.annotation = (np.asarray(annotation).reshape(-1)
                           if annotation is not None else
                           np.full_like(self.mz, None, object))
        if len(self.mz) != len(self.annotation):
            raise ValueError('The mass-to-charge and annotation arrays should '
                             'have equal length')
        # Make sure the fragment peaks are sorted in ascending mass-to-charge
        # order.
        order = np.argsort(self.mz)
        self.mz = self.mz[order]
        self.intensity = self.intensity[order]
        self.annotation = self.annotation[order]

        self.retention_time = retention_time
        self.peptide = peptide
        self.is_decoy = is_decoy

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
        self.annotation = self.annotation[annotation_idx]

        return self

    def set_mz_range(self, min_mz: float, max_mz: float) -> 'MsmsSpectrum':
        """
        Restrict the mass-to-charge ratios of the fragment peaks to the
        given range.

        Parameters
        ----------
        min_mz : float
            Minimum m/z (inclusive).
        max_mz : float
            Maximum m/z (inclusive).

        Returns
        -------
        self : `MsmsSpectrum`
        """
        mz_range_mask = _get_mz_range_mask(self.mz, min_mz, max_mz)
        self.mz = self.mz[mz_range_mask]
        self.intensity = self.intensity[mz_range_mask]
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
        self.annotation = self.annotation[peak_mask]

        return self

    def filter_intensity(self, min_intensity: float = 0.0,
                         max_num_peaks: int = None) -> 'MsmsSpectrum':
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
        max_num_peaks : int, optional
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
        self.annotation = self.annotation[intensity_mask]

        return self

    def scale_intensity(self, scaling: str = None, max_intensity: float = None,
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
                      square root transformation (`degree` is 2). The degree
                      of the root can be specified using the `degree` kwarg.
            - 'log':  Log-transform the peak intensities. The default is a log2
                      transformation (`base` is 2) after summing the
                      intensities with 1 to avoid negative values after the
                      transformation. The base of the logarithm can be
                      specified using the `base` kwarg.
            - 'rank': Rank-transform the peak intensities. The maximum rank of
                      the most intense peak can be specified using the
                      `max_rank` kwarg, by default the number of peaks in the
                      spectrum is used as the maximum rank. Note that
                      `max_rank` should be greater than or equal to the number
                      of peaks in the spectrum.
        max_intensity : float, optional
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

    def annotate_peaks(self, fragment_tol_mass: float, fragment_tol_mode: str,
                       ion_types: str = 'by', max_ion_charge: int = None,
                       peak_assignment: str = 'most_intense')\
            -> 'MsmsSpectrum':
        """

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
        max_ion_charge : int, optional
            All fragments up to and including the given charge will be
            annotated (by default all fragments with a charge up to the
            precursor minus one will be annotated).
        peak_assignment : {'most_intense', 'nearest_mz'}, optional
            In case multiple peaks occur within the given mass window around a
            theoretical peak, only a single peak will be annotated with the
            fragment type:
            - 'most_intense': The most intense peak will be annotated
                              (default).
            - 'nearest_mz':   The peak whose m/z is closest to the theoretical
                              m/z will be annotated.

        Returns
        -------
        self : `MsmsSpectrum`
        """
        if self.peptide is None:
            raise ValueError('No peptide sequence available for the spectrum')
        if max_ion_charge is None:
            max_ion_charge = self.precursor_charge - 1

        theoretical_fragments = _get_theoretical_peptide_fragments(
            self.peptide, ion_types, max_ion_charge)
        self.annotation = np.empty(len(self.mz), object)
        peak_i_start = 0
        for fragment_annotation, fragment_mz in theoretical_fragments:
            while (peak_i_start < len(self.mz) and
                   utils.mass_diff(self.mz[peak_i_start], fragment_mz,
                                   fragment_tol_mode == 'Da')
                   < -fragment_tol_mass):
                peak_i_start += 1
            peak_i_stop = peak_i_start
            annotation_candidates_i = []
            while (peak_i_stop < len(self.mz) and
                   utils.mass_diff(self.mz[peak_i_stop], fragment_mz,
                                   fragment_tol_mode == 'Da')
                   <= fragment_tol_mass):
                annotation_candidates_i.append(peak_i_stop)
                peak_i_stop += 1
            if len(annotation_candidates_i) > 0:
                peak_annotation_i = 0
                if peak_assignment == 'nearest_mz':
                    peak_annotation_i = np.argmin(np.abs(
                        self.mz[peak_i_start: peak_i_stop] - fragment_mz))
                elif peak_assignment == 'most_intense':
                    peak_annotation_i = np.argmax(
                        self.intensity[peak_i_start: peak_i_stop])
                self.annotation[peak_i_start + peak_annotation_i] =\
                    fragment_annotation

        return self
