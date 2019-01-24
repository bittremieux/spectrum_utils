import operator

import numpy as np
from pyteomics import mass
from pyteomics import parser

from spectrum_utils import utils


def _get_theoretical_peptide_fragments(peptide: str, types: str = 'by',
                                       max_charge: int = 1):
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


class FragmentAnnotation:
    """
    Class representing a fragment io annotation.
    """

    def __init__(self, ion_type: str, ion_index: int, charge: int):
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

    def __repr__(self):
        return f"FragmentAnnotation(ion_type='{self.ion_type}', " \
            f"ion_index={self.ion_index}, charge={self.charge})"

    def __str__(self):
        return f'{self.ion_type}{self.ion_index}{"+" * self.charge}'

    def __eq__(self, other):
        if not isinstance(other, FragmentAnnotation):
            return False
        else:
            return (self.ion_type == other.ion_type and
                    self.ion_index == other.ion_index and
                    self.charge == other.charge)



class MsmsSpectrum:
    """
    Class representing a tandem mass spectrum.
    """

    def __init__(self,
                 identifier: str,
                 precursor_mz: float,
                 precursor_charge: int,
                 mz,
                 intensity,
                 annotation=None,
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
        annotation : optional, array_like
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

    def round(self, decimals: int = 0, combine: str = 'sum'):
        """
        Round the mass-to-charge ratios of the fragment peaks to the given
        number of decimals.

        Peaks have the same mass-to-charge ratio after rounding will be
        combined using the specified strategy.

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
        self.mz, unique_indices, unique_inverse = np.unique(
            np.around(self.mz, decimals), True, True)
        # ``around`` sorts the masses. Although ``masses`` is sorted in the
        # constructor, just make sure the intensities and annotations are
        # sorted similarly in case some change occurred in the mean time.
        if len(self.mz) == len(self.intensity) == len(self.annotation):
            self.intensity = self.intensity[unique_indices]
            self.annotation = self.annotation[unique_indices]
        # If peaks got merged by rounding the mass-to-charge ratios we need
        # to combine their intensities and annotations as well.
        else:
            duplicate_indices = np.setdiff1d(
                    np.arange(len(unique_inverse)), unique_indices, True)
            replaced_indices = unique_inverse[duplicate_indices]
            # Combine intensity.
            intensity = self.intensity[unique_indices]
            if combine == 'sum':
                for replaced_i, duplicate_i in zip(replaced_indices,
                                                   duplicate_indices):
                    intensity[replaced_i] += self.intensity[duplicate_i]
            elif combine == 'max':
                for replaced_i, duplicate_i in zip(replaced_indices,
                                                   duplicate_indices):
                    intensity[replaced_i] = max(intensity[replaced_i],
                                                self.intensity[duplicate_i])
            self.intensity = intensity
            # Combine annotation.
            annotation = self.annotation[unique_indices]
            for replaced_i, duplicate_i in zip(replaced_indices,
                                               duplicate_indices):
                if (annotation[replaced_i] is None and
                        self.annotation[duplicate_i] is not None):
                    annotation[replaced_i] = annotation[duplicate_i]
                elif (annotation[replaced_i] is not None and
                        self.annotation[duplicate_i] is not None):
                    annotation[replaced_i] = f'{annotation[replaced_i]} / '\
                                             f'{self.annotation[duplicate_i]}'
            self.annotation = annotation

        return self

    def set_mz_range(self, min_mz: float, max_mz: float):
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
        mz_range = np.where(np.logical_and(min_mz <= self.mz,
                                           self.mz <= max_mz))[0]
        self.mz = self.mz[mz_range]
        self.intensity = self.intensity[mz_range]
        self.annotation = self.annotation[mz_range]

        return self

    def remove_precursor_peak(self, fragment_tol_mass: float,
                              fragment_tol_mode: str):
        """
        Remove fragment peak(s) close to the precursor mass-to-charge ratio.

        Parameters
        ----------
        fragment_tol_mass : float
            Fragment mass tolerance around the precursor mass to remove the
            precursor peak.
        fragment_tol_mode : {'Da', 'ppm'}
            Fragment mass tolerance unit. Either 'Da' or 'ppm'.

        Returns
        -------
        self : `MsmsSpectrum`
        """
        mass_diff = utils.mass_diff(self.mz, self.precursor_mz,
                                    fragment_tol_mode)
        peak_mask = np.where(np.abs(mass_diff) > fragment_tol_mass)[0]
        self.mz = self.mz[peak_mask]
        self.intensity = self.intensity[peak_mask]
        self.annotation = self.annotation[peak_mask]

        return self

    def filter_intensity(self, min_intensity: float = 0.0,
                         max_num_peaks: int = None):
        """
        Remove low-intensity fragment peaks.

        Only the `max_num_peaks` most intense fragment peaks are retained.
        Additionally, noise peaks whose intensity is below `min_intensity`
        percentage of the intensity of the most intense peak are removed.

        Parameters
        ----------
        min_intensity : int, optional
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
        # Only retain the `max_num_peaks` most intense peaks.
        if max_num_peaks is None:
            max_num_peaks = len(self.intensity)
        filter_num_peaks = np.argsort(self.intensity)[::-1][:max_num_peaks]
        # Discard further low-intensity noise peaks.
        max_intensity = self.intensity[filter_num_peaks[0]]
        filter_noise = np.where(
                self.intensity >= min_intensity * max_intensity)[0]

        # Indices get sorted in `intersect1d`, retaining the mz order.
        filter_intensity = np.intersect1d(filter_num_peaks, filter_noise, True)
        self.mz = self.mz[filter_intensity]
        self.intensity = self.intensity[filter_intensity]
        self.annotation = self.annotation[filter_intensity]

        return self

    def scale_intensity(self, scaling: str = None, max_intensity: float = None,
                        **kwargs):
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
            self.intensity = np.power(self.intensity,
                                      1 / kwargs.get('degree', 2))
        elif scaling == 'log':
            self.intensity = (np.log1p(self.intensity) /
                              np.log(kwargs.get('base', 2)))
        elif scaling == 'rank':
            max_rank = kwargs.get('max_rank', len(self.intensity))
            if max_rank < len(self.intensity):
                raise ValueError('`max_rank` should be greater than or equal '
                                 'to the number of peaks in the spectrum. See '
                                 '`filter_intensity` to reduce the number of '
                                 'peaks in the spectrum.')
            self.intensity = ((max_rank -
                               np.argsort(np.argsort(self.intensity)[::-1]))
                              .astype(np.float32))

        if max_intensity is not None:
            self.intensity *= max_intensity / self.intensity.max()

        return self

    def annotate_peaks(self, fragment_tol_mass: float, fragment_tol_mode: str,
                       ion_types: str = 'by', max_ion_charge: int = None,
                       peak_assignment: str = 'most_intense'):
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
                                  fragment_tol_mode) < -fragment_tol_mass):
                peak_i_start += 1
            peak_i_stop = peak_i_start
            annotation_candidates_i = []
            while (peak_i_stop < len(self.mz) and
                   utils.mass_diff(self.mz[peak_i_stop], fragment_mz,
                                  fragment_tol_mode) <= fragment_tol_mass):
                annotation_candidates_i.append(peak_i_stop)
                peak_i_stop += 1
            if len(annotation_candidates_i) > 0:
                if peak_assignment == 'nearest_mz':
                    peak_annotation_i = np.argmin(np.abs(
                        self.mz[peak_i_start: peak_i_stop] - fragment_mz))
                elif peak_assignment == 'most_intense':
                    peak_annotation_i = np.argmax(
                        self.intensity[peak_i_start: peak_i_stop])
                self.annotation[peak_i_start + peak_annotation_i] =\
                    fragment_annotation

        return self
