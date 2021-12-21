import numpy as np
import pytest

from spectrum_utils import proforma, spectrum


@pytest.fixture(autouse=True)
def set_random_seed():
    np.random.seed(13)


def test_fragmentannotation_unknown():
    spectrum.FragmentAnnotation("?")
    with pytest.raises(ValueError):
        spectrum.FragmentAnnotation("?", neutral_loss="-H2O")
    with pytest.raises(ValueError):
        spectrum.FragmentAnnotation("?", isotope=1)
    with pytest.raises(ValueError):
        spectrum.FragmentAnnotation("?", charge=1)
    with pytest.raises(ValueError):
        spectrum.FragmentAnnotation("?", adduct="[M+H]")


def test_fragment_annotation_primary():
    spectrum.FragmentAnnotation(
        "b5", neutral_loss="-H2O", isotope=1, charge=1, adduct="[M+H]"
    )
    with pytest.raises(ValueError):
        spectrum.FragmentAnnotation("b5", charge=0)
    with pytest.raises(ValueError):
        spectrum.FragmentAnnotation("b5", charge=-2)


def test_get_theoretical_fragments():
    peptide = proforma.parse("HPYLEDR")[0]
    fragments = {
        "b1_1": 138.066147,
        "b2_1": 235.118912,
        "b3_1": 398.182220,
        "b4_1": 511.266266,
        "b5_1": 640.308899,
        "b6_1": 755.335815,
        "y1_1": 175.118912,
        "y2_1": 290.145844,
        "y3_1": 419.188446,
        "y4_1": 532.272522,
        "y5_1": 695.335815,
        "y6_1": 792.388550,
        "b1_2": 69.536731,
        "b2_2": 118.063111,
        "b3_2": 199.594776,
        "b4_2": 256.136806,
        "b5_2": 320.658101,
        "b6_2": 378.171571,
        "y1_2": 88.063114,
        "y2_2": 145.576584,
        "y3_2": 210.097879,
        "y4_2": 266.639909,
        "y5_2": 348.171574,
        "y6_2": 396.697954,
        "b1_3": 46.693580,
        "b2_3": 79.044500,
        "b3_3": 133.398943,
        "b4_3": 171.093630,
        "b5_3": 214.107826,
        "b6_3": 252.450140,
        "y1_3": 59.044501,
        "y2_3": 97.386815,
        "y3_3": 140.401011,
        "y4_3": 178.095698,
        "y5_3": 232.450141,
        "y6_3": 264.801061,
    }
    for fragment in spectrum._get_theoretical_fragments(peptide, max_charge=3):
        fragment_mz = fragments[f"{fragment.ion_type}_{fragment.charge}"]
        assert fragment.calc_mz == pytest.approx(fragment_mz)


def test_get_theoretical_fragments_static_mod():
    peptide = proforma.parse("<[+79.96633]@Y>HPYLEDR")[0]
    fragments = {
        "b1_1": 138.066147,
        "b2_1": 235.118912,
        "b3_1": 478.148590,
        "b4_1": 591.232666,
        "b5_1": 720.275269,
        "b6_1": 835.302185,
        "y1_1": 175.118912,
        "y2_1": 290.145844,
        "y3_1": 419.188446,
        "y4_1": 532.272522,
        "y5_1": 775.302185,
        "y6_1": 872.354980,
        "b1_2": 69.536731,
        "b2_2": 118.063111,
        "b3_2": 239.577941,
        "b4_2": 296.119971,
        "b5_2": 360.641266,
        "b6_2": 418.154736,
        "y1_2": 88.063114,
        "y2_2": 145.576584,
        "y3_2": 210.097879,
        "y4_2": 266.639909,
        "y5_2": 388.154739,
        "y6_2": 436.681119,
        "b1_3": 46.693580,
        "b2_3": 79.044500,
        "b3_3": 160.054386,
        "b4_3": 197.749073,
        "b5_3": 240.763270,
        "b6_3": 279.105583,
        "y1_3": 59.044501,
        "y2_3": 97.386815,
        "y3_3": 140.401011,
        "y4_3": 178.095698,
        "y5_3": 259.105585,
        "y6_3": 291.456505,
    }
    for fragment in spectrum._get_theoretical_fragments(peptide, max_charge=3):
        fragment_mz = fragments[f"{fragment.ion_type}_{fragment.charge}"]
        assert fragment.calc_mz == pytest.approx(fragment_mz)


def test_get_theoretical_fragments_mod():
    peptide = proforma.parse("HPY[+79.96633]LEDR")[0]
    fragments = {
        "b1_1": 138.066147,
        "b2_1": 235.118912,
        "b3_1": 478.148590,
        "b4_1": 591.232666,
        "b5_1": 720.275269,
        "b6_1": 835.302185,
        "y1_1": 175.118912,
        "y2_1": 290.145844,
        "y3_1": 419.188446,
        "y4_1": 532.272522,
        "y5_1": 775.302185,
        "y6_1": 872.354980,
        "b1_2": 69.536731,
        "b2_2": 118.063111,
        "b3_2": 239.577941,
        "b4_2": 296.119971,
        "b5_2": 360.641266,
        "b6_2": 418.154736,
        "y1_2": 88.063114,
        "y2_2": 145.576584,
        "y3_2": 210.097879,
        "y4_2": 266.639909,
        "y5_2": 388.154739,
        "y6_2": 436.681119,
        "b1_3": 46.693580,
        "b2_3": 79.044500,
        "b3_3": 160.054386,
        "b4_3": 197.749073,
        "b5_3": 240.763270,
        "b6_3": 279.105583,
        "y1_3": 59.044501,
        "y2_3": 97.386815,
        "y3_3": 140.401011,
        "y4_3": 178.095698,
        "y5_3": 259.105585,
        "y6_3": 291.456505,
    }
    for fragment in spectrum._get_theoretical_fragments(peptide, max_charge=3):
        fragment_mz = fragments[f"{fragment.ion_type}_{fragment.charge}"]
        assert fragment.calc_mz == pytest.approx(fragment_mz)


def test_get_theoretical_fragments_mod_term():
    peptide = proforma.parse("[+42.01056]-HPYLEDR")[0]
    fragments = {
        "b1": 180.076706,
        "b2": 277.129486,
        "b3": 440.192810,
        "b4": 553.276917,
        "b5": 682.319519,
        "b6": 797.346436,
        "y1": 175.118912,
        "y2": 290.145844,
        "y3": 419.188446,
        "y4": 532.272522,
        "y5": 695.335815,
        "y6": 792.388550,
    }
    for fragment in spectrum._get_theoretical_fragments(peptide):
        fragment_mz = fragments[f"{fragment.ion_type}"]
        assert fragment.calc_mz == pytest.approx(fragment_mz)


def test_get_theoretical_fragments_mod_multiple():
    peptide = proforma.parse("[+42.01056]-HPY[+79.96633]LEDR")[0]
    fragments = {
        "b1": 180.076706,
        "b2": 277.129486,
        "b3": 520.159180,
        "b4": 633.243225,
        "b5": 762.285828,
        "b6": 877.312744,
        "y1": 175.118912,
        "y2": 290.145844,
        "y3": 419.188446,
        "y4": 532.272522,
        "y5": 775.302185,
        "y6": 872.354980,
    }
    for fragment in spectrum._get_theoretical_fragments(peptide):
        fragment_mz = fragments[f"{fragment.ion_type}"]
        assert fragment.calc_mz == pytest.approx(fragment_mz)


def test_get_theoretical_fragments_neutral_loss():
    peptide = proforma.parse("HPYLEDR")[0]
    fragments = {
        "b1^1": 138.066147,
        "b2^1": 235.118912,
        "b3^1": 398.182220,
        "b4^1": 511.266266,
        "b5^1": 640.308899,
        "b6^1": 755.335815,
        "y1^1": 175.118912,
        "y2^1": 290.145844,
        "y3^1": 419.188446,
        "y4^1": 532.272522,
        "y5^1": 695.335815,
        "y6^1": 792.388550,
        "b1^2": 69.536731,
        "b2^2": 118.063111,
        "b3^2": 199.594776,
        "b4^2": 256.136806,
        "b5^2": 320.658101,
        "b6^2": 378.171571,
        "y1^2": 88.063114,
        "y2^2": 145.576584,
        "y3^2": 210.097879,
        "y4^2": 266.639909,
        "y5^2": 348.171574,
        "y6^2": 396.697954,
        "b1^3": 46.693580,
        "b2^3": 79.044500,
        "b3^3": 133.398943,
        "b4^3": 171.093630,
        "b5^3": 214.107826,
        "b6^3": 252.450140,
        "y1^3": 59.044501,
        "y2^3": 97.386815,
        "y3^3": 140.401011,
        "y4^3": 178.095698,
        "y5^3": 232.450141,
        "y6^3": 264.801061,
    }
    neutral_loss = "H2O", 18.010565  # water
    neutral_loss_fragments = {}
    for fragment, mz in fragments.items():
        charge = int(fragment.split("^")[1])
        fragment = f"{fragment}-{neutral_loss[0]}"
        neutral_loss_fragments[fragment] = mz - (neutral_loss[1] / charge)
    fragments = {**fragments, **neutral_loss_fragments}
    for fragment in spectrum._get_theoretical_fragments(
        peptide,
        max_charge=3,
        neutral_losses={None: 0, neutral_loss[0]: -neutral_loss[1]},
    ):
        fragment_mz = fragments[
            f"""{fragment.ion_type}^{fragment.charge}{fragment.neutral_loss
            if fragment.neutral_loss is not None else ''}"""
        ]
        assert fragment_mz == pytest.approx(fragment.calc_mz), repr(fragment)


def test_get_theoretical_fragments_mod_neutral_loss():
    peptide = proforma.parse("HPY[+79.96633]LEDR")[0]
    fragments = {
        "b1^1": 138.066147,
        "b2^1": 235.118912,
        "b3^1": 478.148590,
        "b4^1": 591.232666,
        "b5^1": 720.275269,
        "b6^1": 835.302185,
        "y1^1": 175.118912,
        "y2^1": 290.145844,
        "y3^1": 419.188446,
        "y4^1": 532.272522,
        "y5^1": 775.302185,
        "y6^1": 872.354980,
        "b1^2": 69.536731,
        "b2^2": 118.063111,
        "b3^2": 239.577941,
        "b4^2": 296.119971,
        "b5^2": 360.641266,
        "b6^2": 418.154736,
        "y1^2": 88.063114,
        "y2^2": 145.576584,
        "y3^2": 210.097879,
        "y4^2": 266.639909,
        "y5^2": 388.154739,
        "y6^2": 436.681119,
        "b1^3": 46.693580,
        "b2^3": 79.044500,
        "b3^3": 160.054386,
        "b4^3": 197.749073,
        "b5^3": 240.763270,
        "b6^3": 279.105583,
        "y1^3": 59.044501,
        "y2^3": 97.386815,
        "y3^3": 140.401011,
        "y4^3": 178.095698,
        "y5^3": 259.105585,
        "y6^3": 291.456505,
    }
    neutral_loss = "H2O", 18.010565  # water
    neutral_loss_fragments = {}
    for fragment, mz in fragments.items():
        charge = int(fragment.split("^")[1])
        fragment = f"{fragment}-{neutral_loss[0]}"
        neutral_loss_fragments[fragment] = mz - (neutral_loss[1] / charge)
    fragments = {**fragments, **neutral_loss_fragments}
    for fragment in spectrum._get_theoretical_fragments(
        peptide,
        max_charge=3,
        neutral_losses={None: 0, neutral_loss[0]: -neutral_loss[1]},
    ):
        fragment_mz = fragments[
            f"""{fragment.ion_type}^{fragment.charge}{fragment.neutral_loss
            if fragment.neutral_loss is not None else ''}"""
        ]
        assert fragment_mz == pytest.approx(fragment.calc_mz), repr(fragment)


def test_get_theoretical_fragments_ambiguous():
    with pytest.raises(ValueError):
        spectrum._get_theoretical_fragments(proforma.parse("HPYLEBDR")[0])
    with pytest.raises(ValueError):
        spectrum._get_theoretical_fragments(proforma.parse("HPZYLEDR")[0])


def test_get_theoretical_fragments_unsupported_fragment_type():
    with pytest.raises(ValueError):
        spectrum._get_theoretical_fragments(proforma.parse("HPYLEDR")[0], "l")

