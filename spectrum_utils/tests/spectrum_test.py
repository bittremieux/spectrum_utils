import operator

import numpy as np
import pytest
from pyteomics import mass

from spectrum_utils import spectrum


np.random.seed(13)


def test_get_theoretical_fragments():
    peptide = 'HPYLEDR'
    fragments = {'b1_1': 138.066147, 'b2_1': 235.118912, 'b3_1': 398.182220,
                 'b4_1': 511.266266, 'b5_1': 640.308899, 'b6_1': 755.335815,
                 'y1_1': 175.118912, 'y2_1': 290.145844, 'y3_1': 419.188446,
                 'y4_1': 532.272522, 'y5_1': 695.335815, 'y6_1': 792.388550,
                 'b1_2': 69.536731,  'b2_2': 118.063111, 'b3_2': 199.594776,
                 'b4_2': 256.136806, 'b5_2': 320.658101, 'b6_2': 378.171571,
                 'y1_2': 88.063114,  'y2_2': 145.576584, 'y3_2': 210.097879,
                 'y4_2': 266.639909, 'y5_2': 348.171574, 'y6_2': 396.697954,
                 'b1_3': 46.693580,  'b2_3': 79.044500,  'b3_3': 133.398943,
                 'b4_3': 171.093630, 'b5_3': 214.107826, 'b6_3': 252.450140,
                 'y1_3': 59.044501,  'y2_3': 97.386815,  'y3_3': 140.401011,
                 'y4_3': 178.095698, 'y5_3': 232.450141, 'y6_3': 264.801061}
    for fragment in spectrum._get_theoretical_peptide_fragments(
            peptide, max_charge=3):
        fragment_mz = fragments[f'{fragment.ion_type}{fragment.ion_index}'
                                f'_{fragment.charge}']
        assert fragment.calc_mz == pytest.approx(fragment_mz)


def test_get_theoretical_fragments_static_mod():
    spectrum.static_modification('Y', 79.96633)
    peptide = 'HPYLEDR'
    fragments = {'b1_1': 138.066147, 'b2_1': 235.118912, 'b3_1': 478.148590,
                 'b4_1': 591.232666, 'b5_1': 720.275269, 'b6_1': 835.302185,
                 'y1_1': 175.118912, 'y2_1': 290.145844, 'y3_1': 419.188446,
                 'y4_1': 532.272522, 'y5_1': 775.302185, 'y6_1': 872.354980,
                 'b1_2': 69.536731,  'b2_2': 118.063111, 'b3_2': 239.577941,
                 'b4_2': 296.119971, 'b5_2': 360.641266, 'b6_2': 418.154736,
                 'y1_2': 88.063114,  'y2_2': 145.576584, 'y3_2': 210.097879,
                 'y4_2': 266.639909, 'y5_2': 388.154739, 'y6_2': 436.681119,
                 'b1_3': 46.693580,  'b2_3': 79.044500,  'b3_3': 160.054386,
                 'b4_3': 197.749073, 'b5_3': 240.763270, 'b6_3': 279.105583,
                 'y1_3': 59.044501,  'y2_3': 97.386815,  'y3_3': 140.401011,
                 'y4_3': 178.095698, 'y5_3': 259.105585, 'y6_3': 291.456505}
    for fragment in spectrum._get_theoretical_peptide_fragments(
            peptide, None, max_charge=3):
        fragment_mz = fragments[f'{fragment.ion_type}{fragment.ion_index}'
                                f'_{fragment.charge}']
        assert fragment.calc_mz == pytest.approx(fragment_mz)
    assert spectrum._aa_mass['Y'] == pytest.approx(163.06333 + 79.96633)
    spectrum.reset_modifications()
    assert spectrum._aa_mass['Y'] == pytest.approx(163.06333)


def test_get_theoretical_fragments_mod():
    peptide = 'HPYLEDR'
    modifications = {2: 79.96633}
    fragments = {'b1_1': 138.066147, 'b2_1': 235.118912, 'b3_1': 478.148590,
                 'b4_1': 591.232666, 'b5_1': 720.275269, 'b6_1': 835.302185,
                 'y1_1': 175.118912, 'y2_1': 290.145844, 'y3_1': 419.188446,
                 'y4_1': 532.272522, 'y5_1': 775.302185, 'y6_1': 872.354980,
                 'b1_2': 69.536731,  'b2_2': 118.063111, 'b3_2': 239.577941,
                 'b4_2': 296.119971, 'b5_2': 360.641266, 'b6_2': 418.154736,
                 'y1_2': 88.063114,  'y2_2': 145.576584, 'y3_2': 210.097879,
                 'y4_2': 266.639909, 'y5_2': 388.154739, 'y6_2': 436.681119,
                 'b1_3': 46.693580,  'b2_3': 79.044500,  'b3_3': 160.054386,
                 'b4_3': 197.749073, 'b5_3': 240.763270, 'b6_3': 279.105583,
                 'y1_3': 59.044501,  'y2_3': 97.386815,  'y3_3': 140.401011,
                 'y4_3': 178.095698, 'y5_3': 259.105585, 'y6_3': 291.456505}
    for fragment in spectrum._get_theoretical_peptide_fragments(
            peptide, modifications, max_charge=3):
        fragment_mz = fragments[f'{fragment.ion_type}{fragment.ion_index}'
                                f'_{fragment.charge}']
        assert fragment.calc_mz == pytest.approx(fragment_mz)


def test_get_theoretical_fragments_mod_term():
    peptide = 'HPYLEDR'
    modifications = {'N-term': 42.01056}
    fragments = {'b1': 180.076706, 'b2': 277.129486, 'b3': 440.192810,
                 'b4': 553.276917, 'b5': 682.319519, 'b6': 797.346436,
                 'y1': 175.118912, 'y2': 290.145844, 'y3': 419.188446,
                 'y4': 532.272522, 'y5': 695.335815, 'y6': 792.388550}
    for fragment in spectrum._get_theoretical_peptide_fragments(
            peptide, modifications):
        fragment_mz = fragments[f'{fragment.ion_type}{fragment.ion_index}']
        assert fragment.calc_mz == pytest.approx(fragment_mz)


def test_get_theoretical_fragments_mod_multiple():
    peptide = 'HPYLEDR'
    modifications = {'N-term': 42.01056, 2: 79.96633}
    fragments = {'b1': 180.076706, 'b2': 277.129486, 'b3': 520.159180,
                 'b4': 633.243225, 'b5': 762.285828, 'b6': 877.312744,
                 'y1': 175.118912, 'y2': 290.145844, 'y3': 419.188446,
                 'y4': 532.272522, 'y5': 775.302185, 'y6': 872.354980}
    for fragment in spectrum._get_theoretical_peptide_fragments(
            peptide, modifications):
        fragment_mz = fragments[f'{fragment.ion_type}{fragment.ion_index}']
        assert fragment.calc_mz == pytest.approx(fragment_mz)


def test_mz_intensity_len():
    mz = np.random.uniform(100, 1400, 150)
    intensity = np.random.exponential(1, 100)
    with pytest.raises(ValueError):
        spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity)


def test_mz_annotation_len():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    annotation = [spectrum.FragmentAnnotation(None, this_mz, str(this_mz))
                  for this_mz in mz[:100]]
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
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 annotation)
    for this_mz, this_annotation in zip(spec.mz, spec.annotation):
        assert this_mz == pytest.approx(this_annotation.calc_mz)


def test_mz_array():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks).tolist()
    intensity = np.random.lognormal(0, 1, num_peaks)
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 annotation)
    assert type(spec.mz) == np.ndarray
    spec.mz = np.random.uniform(100, 1400, num_peaks).tolist()
    assert type(spec.mz) == np.ndarray
    with pytest.raises(ValueError):
        spec.mz = 14


def test_intensity_array():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks).tolist()
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 annotation)
    assert type(spec.intensity) == np.ndarray
    spec.intensity = np.random.lognormal(0, 1, num_peaks).tolist()
    assert type(spec.intensity) == np.ndarray
    with pytest.raises(ValueError):
        spec.intensity = 14


def test_annotation_array():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 annotation)
    assert type(spec.annotation) == np.ndarray
    spec.annotation = annotation
    assert type(spec.annotation) == np.ndarray
    with pytest.raises(ValueError):
        spec.annotation = 14


def test_annotation_none():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity)
    assert spec.annotation is None
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec.annotation = annotation
    assert type(spec.annotation) == np.ndarray
    spec.annotation = None
    assert spec.annotation is None


def test_init_peptide():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    with pytest.raises(ValueError):
        spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                              peptide='BJOUXZ')


def test_init_modification_none():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    modifications = {0: +1, 4: -56, 'N-term': +16, 'C-term': -42}
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 peptide=None, modifications=modifications)
    assert spec.peptide is None
    assert spec.modifications is None


def test_init_modification_index():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    with pytest.raises(ValueError):
        spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                              peptide='PEPTIDER', modifications={12: 16})


def test_init_modification_str():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    with pytest.raises(ValueError):
        spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                              peptide='PEPTIDER', modifications={'bla': 42})


def test_init_modification():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    modifications = {0: +1, 4: -56, 'N-term': +16, 'C-term': -42}
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 peptide='PEPTIDER',
                                 modifications=modifications)
    assert len(spec.modifications) == 4


def test_round_no_merge():
    num_peaks = 150
    mz = np.arange(1, num_peaks + 1) + np.random.uniform(-0.49, 0.5, num_peaks)
    intensity = np.random.exponential(1, num_peaks)
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz.copy(),
                                 intensity.copy(), annotation.copy())
    decimals = 0
    spec.round(decimals)
    assert len(spec.mz) == num_peaks
    assert len(spec.intensity) == num_peaks
    assert len(spec.annotation) == num_peaks
    np.testing.assert_allclose(spec.mz, np.around(mz, decimals))
    np.testing.assert_allclose(spec.intensity, intensity)


def test_round_merge_len():
    num_peaks = 10
    mz = np.arange(1, num_peaks + 1) + np.random.uniform(-0.2, 0.2, num_peaks)
    mz[4] = mz[3] + 0.0002
    mz[5] = mz[3] + 0.0005
    mz[7] = mz[8] - 0.00037
    intensity = np.random.exponential(1, num_peaks)
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 annotation)
    spec.round(1)
    assert len(spec.mz) == len(mz) - 3
    assert len(spec.mz) == len(spec.intensity)
    assert len(spec.mz) == len(spec.annotation)


def test_round_merge_sum():
    num_peaks = 10
    mz = np.arange(1, num_peaks + 1) + np.random.uniform(-0.2, 0.2, num_peaks)
    mz[4] = mz[3] + 0.0002
    mz[5] = mz[3] + 0.0005
    mz[7] = mz[8] - 0.00037
    intensity = np.random.exponential(1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity.copy())
    spec.round(1, 'sum')
    assert np.sum(spec.intensity) == pytest.approx(np.sum(intensity))


def test_round_merge_max():
    num_peaks = 10
    mz = np.arange(1, num_peaks + 1) + np.random.uniform(-0.2, 0.2, num_peaks)
    mz[4] = mz[3] + 0.0002
    mz[5] = mz[3] + 0.0005
    mz[7] = mz[8] - 0.00037
    intensity = np.arange(1, 11)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity.copy())
    spec.round(1, 'max')
    np.testing.assert_allclose(spec.intensity, [1, 2, 3, 6, 7, 9, 10])


def test_round_merge_annotation():
    num_peaks = 10
    mz = np.arange(1, num_peaks + 1) + np.random.uniform(-0.2, 0.2, num_peaks)
    mz[4] = mz[3] + 0.0002
    mz[5] = mz[3] + 0.0005
    mz[7] = mz[8] - 0.00037
    intensity = np.arange(1, 11)
    annotation = [None, None, None,
                  spectrum.PeptideFragmentAnnotation(1, mz[3], 'b', 4), None,
                  spectrum.PeptideFragmentAnnotation(1, mz[5], 'b', 6),
                  spectrum.PeptideFragmentAnnotation(1, mz[6], 'b', 7), None,
                  None, spectrum.PeptideFragmentAnnotation(1, mz[9], 'b', 10)]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz.copy(), intensity,
                                 annotation.copy())
    spec.round(1, 'max')
    np.testing.assert_array_equal(spec.annotation, [
        None, None, None, spectrum.PeptideFragmentAnnotation(1, mz[5], 'b', 6),
        spectrum.PeptideFragmentAnnotation(1, mz[6], 'b', 7), None,
        spectrum.PeptideFragmentAnnotation(1, mz[9], 'b', 10)])


def test_set_mz_range_keep_all():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 annotation)
    min_mz, max_mz = 0, 1500
    spec.set_mz_range(min_mz, max_mz)
    assert len(spec.mz) == num_peaks
    assert len(spec.intensity) == num_peaks
    assert len(spec.annotation) == num_peaks


def test_set_mz_range_truncate():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 annotation)
    min_mz, max_mz = 400, 1200
    assert spec.mz.min() < min_mz
    assert spec.mz.max() > max_mz
    spec.set_mz_range(min_mz, max_mz)
    assert len(spec.mz) < num_peaks
    assert len(spec.intensity) < num_peaks
    assert len(spec.annotation) < num_peaks
    assert spec.mz.min() >= min_mz
    assert spec.mz.max() <= max_mz


def test_set_mz_range_truncate_left():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 annotation)
    min_mz, max_mz = 400, 1500
    assert spec.mz.min() < min_mz
    spec.set_mz_range(min_mz, max_mz)
    assert len(spec.mz) < num_peaks
    assert len(spec.intensity) < num_peaks
    assert len(spec.annotation) < num_peaks
    assert spec.mz.min() >= min_mz


def test_set_mz_range_truncate_right():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 annotation)
    min_mz, max_mz = 0, 1200
    assert spec.mz.max() > max_mz
    spec.set_mz_range(min_mz, max_mz)
    assert len(spec.mz) < num_peaks
    assert len(spec.intensity) < num_peaks
    assert len(spec.annotation) < num_peaks
    assert spec.mz.max() <= max_mz


def test_set_mz_range_none():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 annotation)
    spec.set_mz_range(None, None)
    assert len(spec.mz) == num_peaks
    assert len(spec.intensity) == num_peaks
    assert len(spec.annotation) == num_peaks


def test_remove_precursor_peak():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    fragment_tol_mass = np.random.uniform(0, 0.5)
    fragment_tol_mode = 'Da'
    precursor_mz = mz[np.random.randint(0, num_peaks)] + fragment_tol_mass / 2
    intensity = np.random.lognormal(0, 1, num_peaks)
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', precursor_mz, 2, mz,
                                 intensity, annotation)
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
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', precursor_mz, 2, mz,
                                 intensity, annotation)
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
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', precursor_mz,
                                 precursor_charge, mz, intensity, annotation)
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
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', precursor_mz,
                                 precursor_charge, mz, intensity, annotation)
    spec.remove_precursor_peak(fragment_tol_mass, fragment_tol_mode, 2)
    assert np.abs(precursor_mz - spec.mz).all() > fragment_tol_mass
    assert len(spec.mz) <= num_peaks - 3
    assert len(spec.intensity) <= num_peaks - 3
    assert len(spec.annotation) <= num_peaks - 3


def test_filter_intensity_keep_all():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 annotation)
    spec.filter_intensity()
    assert len(spec.mz) == num_peaks
    assert len(spec.intensity) == num_peaks
    assert len(spec.annotation) == num_peaks


def test_filter_intensity_remove_low_intensity():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    max_intensity = intensity.max()
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 annotation)
    min_intensity = 0.05
    assert spec.intensity.min() < min_intensity * spec.intensity.max()
    spec.filter_intensity(min_intensity=min_intensity)
    assert len(spec.mz) < num_peaks
    assert len(spec.intensity) < num_peaks
    assert len(spec.annotation) < num_peaks
    assert spec.intensity.max() == pytest.approx(max_intensity)
    assert spec.intensity.min() >= min_intensity * max_intensity


def test_filter_intensity_max_num_peaks():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    max_intensity = intensity.max()
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 annotation)
    max_num_peaks = 50
    spec.filter_intensity(max_num_peaks=max_num_peaks)
    assert len(spec.mz) == max_num_peaks
    assert len(spec.intensity) == max_num_peaks
    assert len(spec.annotation) == max_num_peaks
    assert spec.intensity.max() == pytest.approx(max_intensity)


def test_filter_intensity_remove_low_intensity_max_num_peaks():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    max_intensity = intensity.max()
    annotation = [spectrum.PeptideFragmentAnnotation(1, this_mz, 'y', i)
                  for i, this_mz in enumerate(mz)]
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 2, mz, intensity,
                                 annotation)
    min_intensity = 0.05
    assert spec.intensity.min() < min_intensity * max_intensity
    max_num_peaks = 50
    spec.filter_intensity(min_intensity=min_intensity,
                          max_num_peaks=max_num_peaks)
    assert len(spec.mz) <= max_num_peaks
    assert len(spec.intensity) <= max_num_peaks
    assert len(spec.annotation) <= max_num_peaks
    assert spec.intensity.max() == pytest.approx(max_intensity)
    assert spec.intensity.min() >= min_intensity * max_intensity


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


def test_annotate_peptide_fragments():
    fragment_tol_mass = 0.02
    fragment_tol_mode = 'Da'
    peptides = ['SYELPDGQVITIGNER', 'MFLSFPTTK', 'DLYANTVLSGGTTMYPGIADR',
                'YLYEIAR', 'VAPEEHPVLLTEAPLNPK']
    for charge, peptide in enumerate(peptides, 1):
        fragment_mz = np.asarray([fragment.calc_mz for fragment in
                                  spectrum._get_theoretical_peptide_fragments(
                                      peptide, max_charge=charge - 1)])
        fragment_mz += np.random.uniform(
            -0.9 * fragment_tol_mass, 0.9 * fragment_tol_mass,
            len(fragment_mz))
        num_peaks = 150
        mz = np.random.uniform(100, 1400, num_peaks)
        mz[: len(fragment_mz)] = fragment_mz
        intensity = np.random.lognormal(0, 1, num_peaks)
        spec = spectrum.MsmsSpectrum(
            'test_spectrum', mass.calculate_mass(sequence=peptide,
                                                 charge=charge),
            charge, mz, intensity, peptide=peptide)
        spec.annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode)
        assert np.count_nonzero(spec.annotation) == len(fragment_mz)


def test_annotate_peptide_fragments_nearest_mz():
    fragment_tol_mass = 0.02
    fragment_tol_mode = 'Da'
    peptide = 'YLYEIAR'
    fragment_mz = np.asarray([fragment.calc_mz for fragment in
                              spectrum._get_theoretical_peptide_fragments(
                                  peptide)])
    mz = np.asarray([fragment_mz[0] - 0.005, fragment_mz[0] + 0.015])
    intensity = np.asarray([10, 20])
    charge = 2
    spec = spectrum.MsmsSpectrum(
        'test_spectrum', mass.calculate_mass(sequence=peptide,
                                             charge=charge),
        charge, mz, intensity, peptide=peptide)
    spec.annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode,
                                    peak_assignment='nearest_mz')
    assert spec.annotation[0] == spectrum.PeptideFragmentAnnotation(
        1, fragment_mz[0], 'b', 1)
    assert spec.annotation[1] is None


def test_annotate_peptide_fragments_most_intense():
    fragment_tol_mass = 0.02
    fragment_tol_mode = 'Da'
    peptide = 'YLYEIAR'
    fragment_mz = np.asarray([fragment.calc_mz for fragment in
                              spectrum._get_theoretical_peptide_fragments(
                                  peptide)])
    mz = np.asarray([fragment_mz[0] - 0.01, fragment_mz[0] + 0.01])
    intensity = np.asarray([10, 20])
    charge = 2
    spec = spectrum.MsmsSpectrum(
        'test_spectrum', mass.calculate_mass(sequence=peptide,
                                             charge=charge),
        charge, mz, intensity, peptide=peptide)
    spec.annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode,
                                    peak_assignment='most_intense')
    assert spec.annotation[0] is None
    assert spec.annotation[1] == spectrum.PeptideFragmentAnnotation(
        1, fragment_mz[0], 'b', 1)


def test_annotate_molecule_fragments():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 1, mz, intensity)
    fragment_tol_mass = 0.02
    fragment_tol_mode = 'Da'
    smiles = ['[Cu+2].[O-]S(=O)(=O)[O-] ', 'O=Cc1ccc(O)c(OC)c1',
              'CC(=O)NCCC1=CNc2c1cc(OC)cc2', 'CN1CCC[C@H]1c2cccnc2',
              'CC(C)[C@@]12C[C@@H]1[C@@H](C)C(=O)C2']
    for sm, fragment_mz in zip(smiles, np.random.choice(spec.mz, len(smiles),
                                                        False)):
        spec.annotate_molecule_fragment(sm, fragment_mz, 1, fragment_tol_mass,
                                        fragment_tol_mode)
    assert np.count_nonzero(spec.annotation) == len(smiles)


def test_annotate_molecule_fragments_invalid_mz():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 1, mz, intensity)
    fragment_tol_mass = 0.02
    fragment_tol_mode = 'Da'
    with pytest.raises(ValueError):
        spec.annotate_molecule_fragment('CCCCCCCC', 1600, 1,
                                        fragment_tol_mass, fragment_tol_mode)


def test_annotate_molecule_fragments_nearest_mz():
    fragment_tol_mass = 1
    fragment_tol_mode = 'Da'
    mz, intensity = np.asarray([200, 200.5, 201]), np.asarray([10, 20, 30])
    spec = spectrum.MsmsSpectrum('test_spectrum', 100, 1, mz, intensity)
    fragment_mz = 200.15
    charge = 1
    spec.annotate_molecule_fragment('CCCCCCCC', fragment_mz, charge,
                                    fragment_tol_mass, fragment_tol_mode,
                                    peak_assignment='nearest_mz')
    assert spec.annotation[0] == spectrum.MoleculeFragmentAnnotation(
        charge, fragment_mz, 'CCCCCCCC')
    assert spec.annotation[1] is None
    assert spec.annotation[2] is None


def test_annotate_molecule_fragments_most_intense():
    fragment_tol_mass = 1
    fragment_tol_mode = 'Da'
    mz, intensity = np.asarray([200, 200.5, 201.5]), np.asarray([10, 20, 30])
    spec = spectrum.MsmsSpectrum('test_spectrum', 100, 1, mz, intensity)
    fragment_mz = 200.15
    charge = 1
    spec.annotate_molecule_fragment('CCCCCCCC', fragment_mz, charge,
                                    fragment_tol_mass, fragment_tol_mode,
                                    peak_assignment='most_intense')
    assert spec.annotation[0] is None
    assert spec.annotation[1] == spectrum.MoleculeFragmentAnnotation(
        charge, fragment_mz, 'CCCCCCCC')
    assert spec.annotation[2] is None


def test_annotate_mz_fragments():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 1, mz, intensity)
    fragment_tol_mass = 0.02
    fragment_tol_mode = 'Da'
    num_annotation = 4
    for fragment_mz in np.random.choice(spec.mz, num_annotation, False):
        spec.annotate_mz_fragment(fragment_mz, 1, fragment_tol_mass,
                                  fragment_tol_mode)
    assert np.count_nonzero(spec.annotation) == num_annotation


def test_annotate_mz_fragments_text():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 1, mz, intensity)
    fragment_tol_mass = 0.02
    fragment_tol_mode = 'Da'
    annotation = ['one', 'two', 'three', 'four']
    for a, fragment_mz in zip(annotation, np.random.choice(
            spec.mz, len(annotation), False)):
        spec.annotate_mz_fragment(fragment_mz, 1, fragment_tol_mass,
                                  fragment_tol_mode, a)
    assert np.count_nonzero(spec.annotation) == len(annotation)


def test_annotate_mz_fragments_invalid_mz():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum('test_spectrum', 500, 1, mz, intensity)
    fragment_tol_mass = 0.02
    fragment_tol_mode = 'Da'
    with pytest.raises(ValueError):
        spec.annotate_mz_fragment(1600, 1, fragment_tol_mass,
                                  fragment_tol_mode)


def test_annotate_mz_fragments_nearest_mz():
    fragment_tol_mass = 1
    fragment_tol_mode = 'Da'
    mz, intensity = np.asarray([200, 200.5, 201]), np.asarray([10, 20, 30])
    spec = spectrum.MsmsSpectrum('test_spectrum', 100, 1, mz, intensity)
    fragment_mz = 200.15
    charge = 1
    spec.annotate_mz_fragment(fragment_mz, charge, fragment_tol_mass,
                              fragment_tol_mode, 'nearest_mz', 'hello')
    assert spec.annotation[0] == spectrum.FragmentAnnotation(
        charge, fragment_mz, 'hello')
    assert spec.annotation[1] is None
    assert spec.annotation[2] is None


def test_annotate_mz_fragments_most_intense():
    fragment_tol_mass = 1
    fragment_tol_mode = 'Da'
    mz, intensity = np.asarray([200, 200.5, 201.5]), np.asarray([10, 20, 30])
    spec = spectrum.MsmsSpectrum('test_spectrum', 100, 1, mz, intensity)
    fragment_mz = 200.15
    charge = 1
    spec.annotate_mz_fragment(fragment_mz, charge, fragment_tol_mass,
                              fragment_tol_mode, 'most_intense', 'hello')
    assert spec.annotation[0] is None
    assert spec.annotation[1] == spectrum.FragmentAnnotation(
        charge, fragment_mz, 'hello')
    assert spec.annotation[2] is None
