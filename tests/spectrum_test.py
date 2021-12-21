import operator

import numpy as np
import pytest
from pyteomics import mass

from spectrum_utils import proforma, spectrum


@pytest.fixture(autouse=True)
def set_random_seed():
    np.random.seed(13)


def test_mz_intensity_len():
    mz = np.random.uniform(100, 1400, 150)
    intensity = np.random.exponential(1, 100)
    with pytest.raises(ValueError):
        spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)


def test_init_mz_sorted():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    for mz1, mz2 in zip(spec.mz[:-1], spec.mz[1:]):
        assert mz1 <= mz2


def test_init_intensity_order():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    mz_intensity_tuples = sorted(
        zip(mz, intensity), key=operator.itemgetter(0)
    )
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    for this_mz, this_intensity, mz_intensity_tuple in zip(
        spec.mz, spec.intensity, mz_intensity_tuples
    ):
        assert (this_mz, this_intensity) == pytest.approx(mz_intensity_tuple)


def test_mz_array():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks).tolist()
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    assert type(spec.mz) == np.ndarray
    with pytest.raises(AttributeError):
        spec.mz = np.random.uniform(100, 1400, num_peaks)


def test_intensity_array():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks).tolist()
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    assert type(spec.intensity) == np.ndarray
    with pytest.raises(AttributeError):
        spec.intensity = np.random.lognormal(0, 1, num_peaks)


def test_from_usi():
    usi = ("mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:"
           "VLHPLEGAVVIIFK/2")
    spec = spectrum.MsmsSpectrum.from_usi(usi)
    assert spec.identifier == usi
    assert spec.precursor_mz == 767.9700
    assert spec.precursor_charge == 2
    assert len(spec.mz) == len(spec.intensity)
    usi = "mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555"
    spec = spectrum.MsmsSpectrum.from_usi(usi, precursor_charge=2)
    assert spec.identifier == usi
    assert spec.precursor_mz == 767.9700
    assert spec.precursor_charge == 2
    assert len(spec.mz) == len(spec.intensity)
    usi = "mzspec:PXD003534:DY_HS_Exp7-Ad1:scan:30372"
    precursor_mz, precursor_charge = 507.7484, 2
    spec = spectrum.MsmsSpectrum.from_usi(
        usi, precursor_mz=precursor_mz, precursor_charge=precursor_charge
    )
    assert spec.identifier == usi
    assert spec.precursor_mz == precursor_mz
    assert spec.precursor_charge == precursor_charge
    assert len(spec.mz) == len(spec.intensity)
    with pytest.raises(ValueError):
        spectrum.MsmsSpectrum.from_usi(usi, precursor_mz=precursor_mz)
    with pytest.raises(ValueError):
        spectrum.MsmsSpectrum.from_usi(usi, precursor_charge=precursor_charge)


def test_round_no_merge():
    num_peaks = 150
    mz = np.arange(1, num_peaks + 1) + np.random.uniform(-0.49, 0.5, num_peaks)
    intensity = np.random.exponential(1, num_peaks)
    spec = spectrum.MsmsSpectrum(
        "test_spectrum",
        500,
        2,
        mz.copy(),
        intensity.copy(),
    )
    decimals = 0
    spec.round(decimals)
    assert len(spec.mz) == num_peaks
    assert len(spec.intensity) == num_peaks
    np.testing.assert_allclose(spec.mz, np.around(mz, decimals))
    np.testing.assert_allclose(spec.intensity, intensity)


def test_round_merge_len():
    num_peaks = 10
    mz = np.arange(1, num_peaks + 1) + np.random.uniform(-0.2, 0.2, num_peaks)
    mz[4] = mz[3] + 0.0002
    mz[5] = mz[3] + 0.0005
    mz[7] = mz[8] - 0.00037
    intensity = np.random.exponential(1, num_peaks)
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    spec.annotate_proforma(f"X[+{mz[3]}]", 10, "ppm")
    spec.round(1)
    assert len(spec.mz) == len(mz) - 3
    assert len(spec.mz) == len(spec.intensity)
    assert spec.annotation is None


def test_round_merge_sum():
    num_peaks = 10
    mz = np.arange(1, num_peaks + 1) + np.random.uniform(-0.2, 0.2, num_peaks)
    mz[4] = mz[3] + 0.0002
    mz[5] = mz[3] + 0.0005
    mz[7] = mz[8] - 0.00037
    intensity = np.random.exponential(1, num_peaks)
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity.copy())
    spec.round(1, "sum")
    assert np.sum(spec.intensity) == pytest.approx(np.sum(intensity))


def test_round_merge_max():
    num_peaks = 10
    mz = np.arange(1, num_peaks + 1) + np.random.uniform(-0.2, 0.2, num_peaks)
    mz[4] = mz[3] + 0.0002
    mz[5] = mz[3] + 0.0005
    mz[7] = mz[8] - 0.00037
    intensity = np.arange(1, 11)
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity.copy())
    spec.round(1, "max")
    np.testing.assert_allclose(spec.intensity, [1, 2, 3, 6, 7, 9, 10])


def test_set_mz_range_keep_all():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    min_mz, max_mz = 0, 1500
    spec.set_mz_range(min_mz, max_mz)
    assert len(spec.mz) == num_peaks
    assert len(spec.intensity) == num_peaks


def test_set_mz_range_truncate():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    spec.annotate_proforma(f"X[+{mz[75]}]", 10, "ppm")
    min_mz, max_mz = 400, 1200
    assert spec.mz.min() < min_mz
    assert spec.mz.max() > max_mz
    spec.set_mz_range(min_mz, max_mz)
    assert len(spec.mz) < num_peaks
    assert len(spec.intensity) < num_peaks
    assert spec.annotation is None
    assert spec.mz.min() >= min_mz
    assert spec.mz.max() <= max_mz


def test_set_mz_range_truncate_left():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    min_mz, max_mz = 400, 1500
    assert spec.mz.min() < min_mz
    spec.set_mz_range(min_mz, max_mz)
    assert len(spec.mz) < num_peaks
    assert len(spec.intensity) < num_peaks
    assert spec.mz.min() >= min_mz


def test_set_mz_range_truncate_right():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    min_mz, max_mz = 0, 1200
    assert spec.mz.max() > max_mz
    spec.set_mz_range(min_mz, max_mz)
    assert len(spec.mz) < num_peaks
    assert len(spec.intensity) < num_peaks
    assert spec.mz.max() <= max_mz


def test_set_mz_range_none():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    spec.set_mz_range(None, None)
    assert len(spec.mz) == num_peaks
    assert len(spec.intensity) == num_peaks
    spec.set_mz_range(None, 1500)
    assert len(spec.mz) == num_peaks
    assert len(spec.intensity) == num_peaks
    spec.set_mz_range(0, None)
    assert len(spec.mz) == num_peaks
    assert len(spec.intensity) == num_peaks


def test_set_mz_range_reversed():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    min_mz, max_mz = 400, 1200
    assert spec.mz.min() < min_mz
    assert spec.mz.max() > max_mz
    spec.set_mz_range(max_mz, min_mz)
    assert len(spec.mz) < num_peaks
    assert len(spec.intensity) < num_peaks
    assert spec.mz.min() >= min_mz
    assert spec.mz.max() <= max_mz


def test_remove_precursor_peak():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    fragment_tol_mass = np.random.uniform(0, 0.5)
    fragment_tol_mode = "Da"
    precursor_mz = mz[np.random.randint(0, num_peaks)] + fragment_tol_mass / 2
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum(
        "test_spectrum", precursor_mz, 2, mz, intensity
    )
    spec.annotate_proforma(f"X[+{mz[75]}]", 10, "ppm")
    spec.remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
    assert np.abs(precursor_mz - spec.mz).all() > fragment_tol_mass
    assert len(spec.mz) <= num_peaks - 1
    assert len(spec.intensity) <= num_peaks - 1
    assert spec.annotation is None


def test_remove_precursor_peak_none():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    fragment_tol_mass = np.random.uniform(0, 0.5)
    fragment_tol_mode = "Da"
    precursor_mz = mz[np.random.randint(0, num_peaks)] + fragment_tol_mass * 2
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum(
        "test_spectrum", precursor_mz, 2, mz, intensity
    )
    spec.remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
    assert len(spec.mz) == num_peaks
    assert len(spec.intensity) == num_peaks
    assert np.abs(precursor_mz - spec.mz).all() > fragment_tol_mass


def test_remove_precursor_peak_charge():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    fragment_tol_mass = np.random.uniform(0, 0.5)
    fragment_tol_mode = "Da"
    precursor_mz = mz[np.random.randint(0, num_peaks)] + fragment_tol_mass / 2
    precursor_charge = 3
    mz[-1] = ((precursor_mz - 1.0072766) * precursor_charge) / 2 + 1.0072766
    mz[-2] = ((precursor_mz - 1.0072766) * precursor_charge) + 1.0072766
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum(
        "test_spectrum", precursor_mz, precursor_charge, mz, intensity
    )
    spec.remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
    assert np.abs(precursor_mz - spec.mz).all() > fragment_tol_mass
    assert len(spec.mz) <= num_peaks - 3
    assert len(spec.intensity) <= num_peaks - 3


def test_remove_precursor_peak_isotope():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    fragment_tol_mass = np.random.uniform(0, 0.5)
    fragment_tol_mode = "Da"
    precursor_mz = mz[np.random.randint(0, num_peaks)] + fragment_tol_mass / 2
    precursor_charge = 3
    mz[-1] = precursor_mz + 1 / precursor_charge
    mz[-2] = precursor_mz + 2 / precursor_charge
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum(
        "test_spectrum", precursor_mz, precursor_charge, mz, intensity
    )
    spec.remove_precursor_peak(fragment_tol_mass, fragment_tol_mode, 2)
    assert np.abs(precursor_mz - spec.mz).all() > fragment_tol_mass
    assert len(spec.mz) <= num_peaks - 3
    assert len(spec.intensity) <= num_peaks - 3


def test_filter_intensity_keep_all():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    spec.filter_intensity()
    assert len(spec.mz) == num_peaks
    assert len(spec.intensity) == num_peaks


def test_filter_intensity_remove_low_intensity():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    max_intensity = intensity.max()
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    spec.annotate_proforma(f"X[+{mz[75]}]", 10, "ppm")
    min_intensity = 0.05
    assert spec.intensity.min() < min_intensity * spec.intensity.max()
    spec.filter_intensity(min_intensity=min_intensity)
    assert len(spec.mz) < num_peaks
    assert len(spec.intensity) < num_peaks
    assert spec.annotation is None
    assert spec.intensity.max() == pytest.approx(max_intensity)
    assert spec.intensity.min() >= min_intensity * max_intensity


def test_filter_intensity_max_num_peaks():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    max_intensity = intensity.max()
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    max_num_peaks = 50
    spec.filter_intensity(max_num_peaks=max_num_peaks)
    assert len(spec.mz) == max_num_peaks
    assert len(spec.intensity) == max_num_peaks
    assert spec.intensity.max() == pytest.approx(max_intensity)


def test_filter_intensity_remove_low_intensity_max_num_peaks():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    max_intensity = intensity.max()
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    min_intensity = 0.05
    assert spec.intensity.min() < min_intensity * max_intensity
    max_num_peaks = 50
    spec.filter_intensity(
        min_intensity=min_intensity, max_num_peaks=max_num_peaks
    )
    assert len(spec.mz) <= max_num_peaks
    assert len(spec.intensity) <= max_num_peaks
    assert spec.intensity.max() == pytest.approx(max_intensity)
    assert spec.intensity.min() >= min_intensity * max_intensity


def test_scale_intensity_root():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    for degree in [2, 4, 10]:
        spec = spectrum.MsmsSpectrum(
            "test_spectrum", 500, 2, mz, intensity.copy()
        )
        intensity_unscaled = spec.intensity.copy()
        spec.scale_intensity(scaling="root", degree=degree)
        np.testing.assert_allclose(
            spec.intensity ** degree, intensity_unscaled, rtol=1e-5
        )


def test_scale_intensity_log():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    for base in [2, np.e, 10]:
        spec = spectrum.MsmsSpectrum(
            "test_spectrum", 500, 2, mz, intensity.copy()
        )
        intensity_unscaled = spec.intensity.copy()
        spec.scale_intensity(scaling="log", base=base)
        np.testing.assert_allclose(
            base ** spec.intensity - 1, intensity_unscaled, rtol=1e-5
        )


def test_scale_intensity_rank():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    spec.scale_intensity(scaling="rank")
    np.testing.assert_allclose(
        np.sort(spec.intensity), np.arange(1, num_peaks + 1)
    )


def test_scale_intensity_rank_less_peaks():
    num_peaks = 50
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    max_rank = num_peaks + 50
    spec.scale_intensity(scaling="rank", max_rank=max_rank)
    np.testing.assert_allclose(
        np.sort(spec.intensity), np.arange(num_peaks + 1, max_rank + 1)
    )


def test_scale_intensity_rank_more_peaks():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    with pytest.raises(ValueError):
        spec.scale_intensity(scaling="rank", max_rank=num_peaks - 50)


def test_scale_intensity_max():
    num_peaks = 150
    mz = np.random.uniform(100, 1400, num_peaks)
    intensity = np.random.lognormal(0, 1, num_peaks)
    spec = spectrum.MsmsSpectrum("test_spectrum", 500, 2, mz, intensity)
    intensity_copy, max_intensity = spec.intensity.copy(), spec.intensity.max()
    spec.scale_intensity(max_intensity=1.0)
    assert spec.intensity.max() == pytest.approx(1.0)
    np.testing.assert_allclose(
        spec.intensity * max_intensity, intensity_copy, rtol=1e-5
    )


def test_annotate_proforma():
    fragment_tol_mass, fragment_tol_mode = 0.02, "Da"
    peptides = [
        "SYELPDGQVITIGNER",
        "MFLSFPTTK",
        "DLYANTVLSGGTTMYPGIADR",
        "YLYEIAR",
        "VAPEEHPVLLTEAPLNPK",
    ]
    for charge, peptide in enumerate(peptides, 1):
        fragment_mz = np.asarray(
            [
                fragment.calc_mz
                for fragment in spectrum._get_theoretical_fragments(
                    proforma.parse(peptide)[0], max_charge=charge - 1
                )
            ]
        )
        fragment_mz += np.random.uniform(
            -0.9 * fragment_tol_mass, 0.9 * fragment_tol_mass, len(fragment_mz)
        )
        num_peaks = 150
        mz = np.random.uniform(100, 1400, num_peaks)
        mz[: len(fragment_mz)] = fragment_mz
        intensity = np.random.lognormal(0, 1, num_peaks)
        spec = spectrum.MsmsSpectrum(
            "test_spectrum",
            mass.calculate_mass(sequence=peptide, charge=charge),
            charge,
            mz,
            intensity,
        )
        spec.annotate_proforma(peptide, fragment_tol_mass, fragment_tol_mode)
        assert np.count_nonzero(spec.annotation) >= len(fragment_mz)


def test_annotate_proforma_nearest_mz():
    fragment_tol_mass = 0.02
    fragment_tol_mode = "Da"
    peptide = "YLYEIAR"
    fragment_mz = np.asarray(
        [
            fragment.calc_mz
            for fragment in spectrum._get_theoretical_fragments(
                proforma.parse(peptide)[0]
            )
        ]
    )
    mz = np.asarray([fragment_mz[0] - 0.005, fragment_mz[0] + 0.015])
    intensity = np.asarray([10, 20])
    charge = 2
    spec = spectrum.MsmsSpectrum(
        "test_spectrum",
        mass.calculate_mass(sequence=peptide, charge=charge),
        charge,
        mz,
        intensity,
    )
    spec.annotate_proforma(
        peptide,
        fragment_tol_mass,
        fragment_tol_mode,
        peak_assignment="nearest_mz",
    )
    assert spec.annotation[0] == spectrum.FragmentAnnotation(
        "b1", charge=1, calc_mz=fragment_mz[0]
    )
    assert spec.annotation[1] is None


def test_annotate_proforma_most_intense():
    fragment_tol_mass = 0.02
    fragment_tol_mode = "Da"
    peptide = "YLYEIAR"
    fragment_mz = np.asarray(
        [
            fragment.calc_mz
            for fragment in spectrum._get_theoretical_fragments(
                proforma.parse(peptide)[0]
            )
        ]
    )
    mz = np.asarray([fragment_mz[0] - 0.01, fragment_mz[0] + 0.01])
    intensity = np.asarray([10, 20])
    charge = 2
    spec = spectrum.MsmsSpectrum(
        "test_spectrum",
        mass.calculate_mass(sequence=peptide, charge=charge),
        charge,
        mz,
        intensity,
    )
    spec.annotate_proforma(
        peptide,
        fragment_tol_mass,
        fragment_tol_mode,
        peak_assignment="most_intense",
    )
    assert spec.annotation[0] is None
    assert spec.annotation[1] == spectrum.FragmentAnnotation(
        "b1", charge=1, calc_mz=fragment_mz[0]
    )


def test_annotate_proforma_neutral_loss():
    fragment_tol_mass, fragment_tol_mode = 0.02, "Da"
    neutral_loss = "H2O", 18.010565  # water
    n_peaks = 150
    peptides = [
        "SYELPDGQVITIGNER",
        "MFLSFPTTK",
        "DLYANTVLSGGTTMYPGIADR",
        "YLYEIAR",
        "VAPEEHPVLLTEAPLNPK",
    ]
    for charge, peptide in enumerate(peptides, 2):
        fragment_mz = np.asarray(
            [
                fragment.calc_mz
                for fragment in spectrum._get_theoretical_fragments(
                    proforma.parse(peptide)[0],
                    neutral_losses={
                        None: 0,
                        neutral_loss[0]: -neutral_loss[1],
                    },
                )
            ]
        )
        fragment_mz += np.random.uniform(
            -0.9 * fragment_tol_mass, 0.9 * fragment_tol_mass, len(fragment_mz)
        )
        mz = np.random.uniform(100, 1400, n_peaks)
        mz[: len(fragment_mz)] = fragment_mz
        intensity = np.random.lognormal(0, 1, n_peaks)
        spec = spectrum.MsmsSpectrum(
            "test_spectrum",
            mass.calculate_mass(sequence=peptide, charge=charge),
            charge,
            mz,
            intensity,
        )
        spec.annotate_proforma(
            peptide,
            fragment_tol_mass,
            fragment_tol_mode,
            neutral_losses={neutral_loss[0]: -neutral_loss[1]},
        )
        n_fragments = (
            len(fragment_mz)
            - (
                len(fragment_mz)
                - (
                    np.partition(
                        np.abs(
                            fragment_mz.reshape(-1, 1)
                            - fragment_mz.reshape(1, -1)
                        ),
                        1,
                        axis=1,
                    )[:, 1]
                    >= fragment_tol_mass
                ).sum()
            )
            // 2
        )
        assert np.count_nonzero(spec.annotation) >= n_fragments
