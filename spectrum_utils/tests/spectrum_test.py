import operator

import numpy as np
import pytest
from pyteomics import mass
from pyteomics import parser

from spectrum_utils import spectrum


np.random.seed(13)


def test_mz_intensity_len():
    mz = np.random.uniform(100, 1400, 150)
    intensity = np.random.exponential(1, 100)
    with pytest.raises(ValueError):
        spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity)


def test_mz_annotation_len():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    annotation = [str(this_mz) for this_mz in mz[:100]]
    with pytest.raises(ValueError):
        spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                              annotation)


def test_init_mz_sorted():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity)
    for mz1, mz2 in zip(spec.mz[:-1], spec.mz[1:]):
        assert mz1 <= mz2


def test_init_intensity_order():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    mz_intensity_tuples = sorted(zip(mz, intensity),
                                 key=operator.itemgetter(0))
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity)
    for this_mz, this_intensity, mz_intensity_tuple in zip(
            spec.mz, spec.intensity, mz_intensity_tuples):
        assert (this_mz, this_intensity) == pytest.approx(mz_intensity_tuple)


def test_init_annotation_order():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    annotation = [str(this_mz) for this_mz in mz]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 annotation)
    for this_mz, this_annotation in zip(spec.mz, spec.annotation):
        assert this_mz == pytest.approx(float(this_annotation))


def test_round_no_merge():
    num_peaks = 150
    mz = np.arange(1, num_peaks + 1) + np.random.uniform(-0.49, 0.5, num_peaks)
    intensity = np.random.exponential(1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz.copy(),
                                 intensity.copy())
    spec.round(0)
    assert len(spec.mz) == num_peaks
    assert len(spec.intensity) == num_peaks
    assert len(spec.annotation) == num_peaks
    np.testing.assert_allclose(spec.mz, mz)
    np.testing.assert_allclose(spec.intensity, intensity)


def test_round_merge_len():
    num_peaks = 10
    mz = np.arange(1, num_peaks + 1) + np.random.uniform(-0.2, 0.2, num_peaks)
    mz[4] = mz[3] + 0.002
    mz[5] = mz[3] + 0.005
    mz[7] = mz[8] - 0.0037
    intensity = np.random.exponential(1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity)
    spec.round(1)
    assert len(spec.mz) == len(mz) - 3
    assert len(spec.mz) == len(spec.intensity)
    assert len(spec.mz) == len(spec.annotation)


def test_round_merge_sum():
    num_peaks = 10
    mz = np.arange(1, num_peaks + 1) + np.random.uniform(-0.2, 0.2, num_peaks)
    mz[4] = mz[3] + 0.002
    mz[5] = mz[3] + 0.005
    mz[7] = mz[8] - 0.0037
    intensity = np.random.exponential(1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity.copy())
    spec.round(1, 'sum')
    assert np.sum(spec.intensity) == pytest.approx(np.sum(intensity))


def test_round_merge_max():
    num_peaks = 10
    mz = np.arange(1, num_peaks + 1) + np.random.uniform(-0.2, 0.2, num_peaks)
    mz[4] = mz[3] + 0.002
    mz[5] = mz[3] + 0.005
    mz[7] = mz[8] - 0.0037
    intensity = np.arange(1, 11)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity.copy())
    spec.round(1, 'max')
    np.testing.assert_allclose(spec.intensity, [1, 2, 3, 6, 7, 9, 10])


def test_round_merge_annotation():
    num_peaks = 10
    mz = np.arange(1, num_peaks + 1) + np.random.uniform(-0.2, 0.2, num_peaks)
    mz[4] = mz[3] + 0.002
    mz[5] = mz[3] + 0.005
    mz[7] = mz[8] - 0.0037
    intensity = np.arange(1, 11)
    annotation = [None, None, None, spectrum.FragmentAnnotation('b', 4, 1),
                  None, spectrum.FragmentAnnotation('b', 6, 1),
                  spectrum.FragmentAnnotation('b', 7, 1), None, None,
                  spectrum.FragmentAnnotation('b', 10, 1)]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 annotation.copy())
    spec.round(1, 'max')
    np.testing.assert_array_equal(spec.annotation, [
        None, None, None, spectrum.FragmentAnnotation('b', 6, 1),
        spectrum.FragmentAnnotation('b', 7, 1), None,
        spectrum.FragmentAnnotation('b', 10, 1)])


def test_set_mz_range_keep_all():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity)
    min_mz, max_mz = 0, 1500
    spec.set_mz_range(min_mz, max_mz)
    assert len(spec.mz) == num_peaks
    assert len(spec.intensity) == num_peaks
    assert len(spec.annotation) == num_peaks


def test_set_mz_range_truncate():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity)
    min_mz, max_mz = 400, 1200
    assert spec.mz.min() < min_mz
    assert spec.mz.max() > max_mz
    spec.set_mz_range(min_mz, max_mz)
    assert len(spec.mz) < num_peaks
    assert len(spec.intensity) < num_peaks
    assert len(spec.annotation) < num_peaks
    assert spec.mz.min() >= min_mz
    assert spec.mz.max() <= max_mz


def test_remove_precursor_peak():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    fragment_tol_mass = np.random.uniform(0, 0.5)
    fragment_tol_mode = 'Da'
    precursor_mz = mz[np.random.randint(0, num_peaks)] + fragment_tol_mass / 2
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', precursor_mz, 2, mz,
                                 intensity)
    spec.remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
    assert np.abs(precursor_mz - spec.mz).all() > fragment_tol_mass
    assert len(spec.mz) <= num_peaks - 1
    assert len(spec.intensity) <= num_peaks - 1
    assert len(spec.annotation) <= num_peaks - 1


def test_remove_precursor_peak_none():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    fragment_tol_mass = np.random.uniform(0, 0.5)
    fragment_tol_mode = 'Da'
    precursor_mz = mz[np.random.randint(0, num_peaks)] + fragment_tol_mass * 2
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', precursor_mz, 2, mz,
                                 intensity)
    spec.remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
    assert len(spec.mz) == num_peaks
    assert len(spec.intensity) == num_peaks
    assert len(spec.annotation) == num_peaks
    assert np.abs(precursor_mz - spec.mz).all() > fragment_tol_mass


def test_remove_precursor_peak_charge():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    fragment_tol_mass = np.random.uniform(0, 0.5)
    fragment_tol_mode = 'Da'
    precursor_mz = mz[np.random.randint(0, num_peaks)] + fragment_tol_mass / 2
    precursor_charge = 3
    mz[-1] = ((precursor_mz - 1.0072766) * precursor_charge) / 2 + 1.0072766
    mz[-2] = ((precursor_mz - 1.0072766) * precursor_charge) + 1.0072766
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', precursor_mz,
                                 precursor_charge, mz, intensity)
    spec.remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
    assert np.abs(precursor_mz - spec.mz).all() > fragment_tol_mass
    assert len(spec.mz) <= num_peaks - 3
    assert len(spec.intensity) <= num_peaks - 3
    assert len(spec.annotation) <= num_peaks - 3


def test_remove_precursor_peak_isotope():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    fragment_tol_mass = np.random.uniform(0, 0.5)
    fragment_tol_mode = 'Da'
    precursor_mz = mz[np.random.randint(0, num_peaks)] + fragment_tol_mass / 2
    precursor_charge = 3
    mz[-1] = precursor_mz + 1 / precursor_charge
    mz[-2] = precursor_mz + 2 / precursor_charge
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', precursor_mz,
                                 precursor_charge, mz, intensity)
    spec.remove_precursor_peak(fragment_tol_mass, fragment_tol_mode, 2)
    assert np.abs(precursor_mz - spec.mz).all() > fragment_tol_mass
    assert len(spec.mz) <= num_peaks - 3
    assert len(spec.intensity) <= num_peaks - 3
    assert len(spec.annotation) <= num_peaks - 3


def test_filter_intensity_keep_all():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity)
    spec.filter_intensity()
    assert len(spec.mz) == num_peaks
    assert len(spec.intensity) == num_peaks
    assert len(spec.annotation) == num_peaks


def test_filter_intensity_remove_low_intensity():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity)
    min_intensity = 0.05
    assert spec.intensity.min() < min_intensity * spec.intensity.max()
    spec.filter_intensity(min_intensity=min_intensity)
    assert len(spec.mz) < num_peaks
    assert len(spec.intensity) < num_peaks
    assert len(spec.annotation) < num_peaks
    assert spec.intensity.min() >= min_intensity * spec.intensity.max()


def test_filter_intensity_max_num_peaks():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity)
    max_num_peaks = 50
    spec.filter_intensity(max_num_peaks=max_num_peaks)
    assert len(spec.mz) == max_num_peaks
    assert len(spec.intensity) == max_num_peaks
    assert len(spec.annotation) == max_num_peaks


def test_filter_intensity_remove_low_intensity_max_num_peaks():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity)
    min_intensity = 0.05
    assert spec.intensity.min() < min_intensity * spec.intensity.max()
    max_num_peaks = 50
    spec.filter_intensity(min_intensity=min_intensity,
                          max_num_peaks=max_num_peaks)
    assert len(spec.mz) <= max_num_peaks
    assert len(spec.intensity) <= max_num_peaks
    assert len(spec.annotation) <= max_num_peaks
    assert spec.intensity.min() >= min_intensity * spec.intensity.max()


def test_scale_intensity_root():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    for degree in [2, 4, 10]:
        spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz,
                                     intensity.copy())
        intensity_unscaled = spec.intensity.copy()
        spec.scale_intensity(scaling='root', degree=degree)
        np.testing.assert_allclose(spec.intensity ** degree,
                                   intensity_unscaled, rtol=1e-5)


def test_scale_intensity_log():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    for base in [2, np.e, 10]:
        spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz,
                                     intensity.copy())
        intensity_unscaled = spec.intensity.copy()
        spec.scale_intensity(scaling='log', base=base)
        np.testing.assert_allclose(base ** spec.intensity - 1,
                                   intensity_unscaled, rtol=1e-5)


def test_scale_intensity_rank():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity)
    spec.scale_intensity(scaling='rank')
    np.testing.assert_allclose(np.sort(spec.intensity),
                               np.arange(1, num_peaks + 1))


def test_scale_intensity_rank_less_peaks():
    num_peaks = 50
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity)
    max_rank = num_peaks + 50
    spec.scale_intensity(scaling='rank', max_rank=max_rank)
    np.testing.assert_allclose(np.sort(spec.intensity),
                               np.arange(num_peaks + 1, max_rank + 1))


def test_scale_intensity_rank_more_peaks():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity)
    with pytest.raises(ValueError):
        spec.scale_intensity(scaling='rank', max_rank=num_peaks - 50)


def test_scale_intensity_max():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity)
    intensity_copy, max_intensity = spec.intensity.copy(), spec.intensity.max()
    spec.scale_intensity(max_intensity=1.)
    assert spec.intensity.max() == pytest.approx(1.)
    np.testing.assert_allclose(spec.intensity * max_intensity, intensity_copy,
                               rtol=1e-5)


def test_annotate_peaks():
    fragment_tol_mass = 0.02
    fragment_tol_mode = 'Da'
    peptides = ['SYELPDGQVITIGNER', 'MFLSFPTTK', 'DLYANTVLSGGTTMYPGIADR',
                'YLYEIAR', 'VAPEEHPVLLTEAPLNPK']
    for peptide in peptides:
        fragment_mz = np.asarray([mz for _, mz in
                                  spectrum._get_theoretical_peptide_fragments(
                                      peptide)])
        fragment_mz += np.random.uniform(
            -fragment_tol_mass, fragment_tol_mass, len(fragment_mz))
        num_peaks = 150
        mz = np.random.uniform(100, 1400, num_peaks)
        mz[: len(fragment_mz)] = fragment_mz
        intensity = np.random.lognormal(0, 1, num_peaks)
        charge = 2
        spec = spectrum.MsmsSpectrum(
            'test_spectrum', mass.calculate_mass(sequence=peptide,
                                                 charge=charge),
            charge, mz, intensity, peptide=peptide)
        spec.annotate_peaks(fragment_tol_mass, fragment_tol_mode)
        assert np.count_nonzero(spec.annotation) == len(fragment_mz)


def test_annotate_peaks_nearest_mz():
    fragment_tol_mass = 0.02
    fragment_tol_mode = 'Da'
    peptide = 'YLYEIAR'
    fragment_mz = np.asarray([mz for _, mz in
                              spectrum._get_theoretical_peptide_fragments(
                                  peptide)])
    mz = np.asarray([fragment_mz[0] - 0.005, fragment_mz[0] + 0.015])
    intensity = np.asarray([10, 20])
    charge = 2
    spec = spectrum.MsmsSpectrum(
        'test_spectrum', mass.calculate_mass(sequence=peptide,
                                             charge=charge),
        charge, mz, intensity, peptide=peptide)
    spec.annotate_peaks(fragment_tol_mass, fragment_tol_mode,
                        peak_assignment='nearest_mz')
    assert spec.annotation[0] is not None
    assert spec.annotation[1] is None


def test_annotate_peaks_most_intense():
    fragment_tol_mass = 0.02
    fragment_tol_mode = 'Da'
    peptide = 'YLYEIAR'
    fragment_mz = np.asarray([mz for _, mz in
                              spectrum._get_theoretical_peptide_fragments(
                                peptide)])
    mz = np.asarray([fragment_mz[0] - 0.01, fragment_mz[0] + 0.01])
    intensity = np.asarray([10, 20])
    charge = 2
    spec = spectrum.MsmsSpectrum(
        'test_spectrum', mass.calculate_mass(sequence=peptide,
                                             charge=charge),
        charge, mz, intensity, peptide=peptide)
    spec.annotate_peaks(fragment_tol_mass, fragment_tol_mode,
                        peak_assignment='most_intense')
    assert spec.annotation[0] is None
    assert spec.annotation[1] is not None
