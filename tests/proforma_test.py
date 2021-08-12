import math
import os
import shutil

import pytest

from spectrum_utils import proforma


def test_proforma_canonical():
    assert proforma.parse('PEPTIDE') == ('PEPTIDE', {})
    assert proforma.parse('ACDEFGHIKL') == ('ACDEFGHIKL', {})
    assert proforma.parse('MNPQRSTVWY') == ('MNPQRSTVWY', {})
    assert proforma.parse('PEPTJDE') == ('PEPTJDE', {})
    assert proforma.parse('PEPTIOE') == ('PEPTIOE', {})
    assert proforma.parse('PEUTIDE') == ('PEUTIDE', {})
    assert proforma.parse('PEPTIXDE') == ('PEPTIXDE', {})


def test_proforma_name():
    # Unimod without prefix.
    assert (proforma.parse('EM[Oxidation]EVEES[Phospho]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 79.966331}))
    # Unimod with prefix.
    assert (proforma.parse('EM[U:Oxidation]EVEES[U:Phospho]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 79.966331}))
    # Unimod with internal brackets in the name.
    assert (proforma.parse('EM[Oxidation]EVE[Cation:Mg[II]]ES[Phospho]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 4: 21.969392, 6: 79.966331}))
    # PSI-MOD without prefix.
    assert (proforma.parse('EM[L-methionine sulfoxide]EVEE'
                           'S[O-phospho-L-serine]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 79.966331}))
    # PSI-MOD with prefix.
    assert (proforma.parse('EM[M:L-methionine sulfoxide]EVEE'
                           'S[M:O-phospho-L-serine]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 79.966331}))
    # RESID (mandatory prefix).
    assert (proforma.parse('EM[R:L-methionine sulfone]EVEE'
                           'S[R:O-phospho-L-serine]PEK') ==
            ('EMEVEESPEK', {1: 31.989829, 6: 79.966331}))
    # XL-MOD (mandatory prefix).
    assert (proforma.parse('EMEVTK[X:DSS#XL1]SESPEK') ==
            ('EMEVTKSESPEK', {5: 138.06807961}))
    # GNO (mandatory prefix).
    assert (proforma.parse('NEEYN[G:G59626AS]K') ==
            ('NEEYNK', {4: 1931.69}))


def test_proforma_accession():
    # Valid accessions.
    assert (proforma.parse('EM[UNIMOD:35]EVEES[UNIMOD:56]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 45.029395}))
    assert (proforma.parse('EM[MOD:00719]EVEES[MOD:00046]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 79.966331}))
    assert (proforma.parse('EM[RESID:AA0581]EVEES[RESID:AA0037]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 79.966331}))
    assert (proforma.parse('NEEYN[GNO:G59626AS]K') ==
            ('NEEYNK', {4: 1931.69}))
    assert (proforma.parse('YPVLN[GNO:G62765YT]VTMPN[GNO:G02815KT]NSNGKFDK') ==
            ('YPVLNVTMPNNSNGKFDK', {4: 1720.59, 9: 1234.43}))
    # Invalid accessions.
    with pytest.raises(KeyError):
        proforma.parse('EM[U:35]EVEES[U:56]PEK')
    with pytest.raises(KeyError):
        proforma.parse('EM[M:00719]EVEES[M:00046]PEK')
    with pytest.raises(KeyError):
        proforma.parse('EM[R:AA0581]EVEES[R:AA0037]PEK')


def test_proforma_xlink():
    # Crosslink notation.
    assert (proforma.parse('EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1]') ==
            ('EMEVTKSESPEK', {5: 138.06807961}))
    assert (proforma.parse('EMK[XLMOD:02000#XL1]EVTKSE[XLMOD:02010#XL2]S'
                           'K[#XL1]PEK[#XL2]AR') ==
            ('EMKEVTKSESKPEKAR', {2: 138.06807961, 8: -18.01056027}))
    assert (proforma.parse('EMEVTK[XLMOD:02001#XL1]SESPEK') ==
            ('EMEVTKSESPEK', {5: 138.06807961}))
    assert (proforma.parse('EMEVTK[XLMOD:02001]SESPEK') ==
            ('EMEVTKSESPEK', {5: 138.06807961}))
    # Inter-chain crosslinks.
    with pytest.raises(NotImplementedError):
        proforma.parse(
            'SEK[XLMOD:02001#XL1]UENCE\\EMEVTK[XLMOD:02001#XL1]SESPEK')
    with pytest.raises(NotImplementedError):
        proforma.parse('SEK[XLMOD:02001#XL1]UENCE\\EMEVTK[#XL1]SESPEK')
    # Disulfide linkages.
    assert (proforma.parse('EVTSEKC[MOD:00034#XL1]LEMSC[#XL1]EFD') ==
            ('EVTSEKCLEMSCEFD', {6: -2.015650}))
    assert (proforma.parse('EVTSEKC[L-cystine (cross-link)#XL1]LEMSC[#XL1]EFD')
            == ('EVTSEKCLEMSCEFD', {6: -2.015650}))
    with pytest.raises(NotImplementedError):
        proforma.parse('FVNQHLC[MOD:00034#XL1]GSHLVEALYLVC[MOD:00034#XL2]GERGF'
                       'FYTPKA\\GIVEQC[MOD:00034#XL3]C[#XL1]TSIC[#XL3]SLYQLENY'
                       'C[#XL2]N')
    assert (proforma.parse('EVTSEKC[XLMOD:02009#XL1]LEMSC[#XL1]EFD') ==
            ('EVTSEKCLEMSCEFD', {6: -2.01565007}))
    assert (proforma.parse('EVTSEKC[X:Disulfide#XL1]LEMSC[#XL1]EFD') ==
            ('EVTSEKCLEMSCEFD', {6: -2.01565007}))
    assert (proforma.parse('EVTSEKC[half cystine]LEMSC[half cystine]EFD') ==
            ('EVTSEKCLEMSCEFD', {6: -1.007825, 11: -1.007825}))
    assert (proforma.parse('EVTSEKC[MOD:00798]LEMSC[MOD:00798]EFDEVTSEK'
                           'C[MOD:00798]LEMSC[MOD:00798]EFD') ==
            ('EVTSEKCLEMSCEFDEVTSEKCLEMSCEFD', {6: -1.007825, 11: -1.007825,
                                                21: -1.007825, 26: -1.007825}))


def test_proforma_delta_mass():
    # No prefix.
    assert (proforma.parse('EM[+15.9949]EVEES[+79.9663]PEK') ==
            ('EMEVEESPEK', {1: 15.9949, 6: 79.9663}))
    assert (proforma.parse('EM[+15.995]EVEES[-18.01]PEK') ==
            ('EMEVEESPEK', {1: 15.995, 6: -18.01}))
    # CV prefix.
    assert (proforma.parse('EM[U:+15.9949]EVEES[U:+79.9663]PEK') ==
            ('EMEVEESPEK', {1: 15.9949, 6: 79.9663}))
    assert (proforma.parse('EM[U:+15.995]EVEES[U:+79.966]PEK') ==
            ('EMEVEESPEK', {1: 15.995, 6: 79.966}))
    # Experimentally observed prefix.
    assert (proforma.parse('EM[U:+15.995]EVEES[Obs:+79.978]PEK') ==
            ('EMEVEESPEK', {1: 15.995, 6: 79.978}))


def test_proforma_gap():
    assert proforma.parse('RTAAX[+367.0537]WT') == ('RTAAXWT', {4: 367.0537})


def test_proforma_formula():
    sequence, mass_diffs = proforma.parse('SEQUEN[Formula:C12H20O2]CE')
    assert sequence == 'SEQUENCE'
    assert len(mass_diffs) == 1
    assert math.isclose(mass_diffs[5], 196.1463, abs_tol=0.00005)
    sequence, mass_diffs = proforma.parse('SEQUEN[Formula:C12 H20 O2]CE')
    assert sequence == 'SEQUENCE'
    assert len(mass_diffs) == 1
    assert math.isclose(mass_diffs[5], 196.1463, abs_tol=0.00005)


def test_proforma_glycan():
    assert (proforma.parse('SEQUEN[Glycan:HexNAc1Hex2]CE') ==
            ('SEQUENCE', {5: 527.185019356}))
    assert (proforma.parse('SEQUEN[Glycan:HexNAc1 Hex2]CE') ==
            ('SEQUENCE', {5: 527.185019356}))


def test_proforma_special():
    # N-terminal and C-terminal modifications.
    assert (proforma.parse('[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK') ==
            ('EMEVNESPEK', {'N-term': 144.102063, 1: 15.994915, 6: 79.966331}))
    assert (proforma.parse('[iTRAQ4plex]-EM[U:Oxidation]EVNES[Phospho]PE'
                           'K[iTRAQ4plex]-[Methyl]') ==
            ('EMEVNESPEK', {'N-term': 144.102063, 1: 15.994915, 6: 79.966331,
                            9: 144.102063, 'C-term': 14.01565}))
    # Labile modifications.
    assert (proforma.parse('{Glycan:Hex}EM[U:Oxidation]EVNES[Phospho]PE'
                           'K[iTRAQ4plex]') ==
            ('EMEVNESPEK', {1: 15.994915, 6: 79.966331, 9: 144.102063}))
    assert (proforma.parse('{Glycan:Hex}[iTRAQ4plex]-EM[Oxidation]EVNE'
                           'S[Phospho]PEK[iTRAQ4plex]') ==
            ('EMEVNESPEK', {'N-term': 144.102063, 1: 15.994915, 6: 79.966331,
                            9: 144.102063}))
    assert (proforma.parse('{Glycan:Hex}[iTRAQ4plex]-EM[Oxidation]EVNE'
                           'S[Phospho]PEK[iTRAQ4plex]-[Methyl]') ==
            ('EMEVNESPEK', {'N-term': 144.102063, 1: 15.994915, 6: 79.966331,
                            9: 144.102063, 'C-term': 14.01565}))
    assert (proforma.parse('{Glycan:Hex}{Glycan:NeuAc}EMEVNESPEK') ==
            ('EMEVNESPEK', {}))


def test_proforma_ambiguous_position():
    # Unknown modification position.
    assert (proforma.parse('[Phospho]?EM[Oxidation]EVTSESPEK') ==
            ('EMEVTSESPEK', {1: 15.994915}))
    assert (proforma.parse(
                '[Phospho][Phospho]?[Acetyl]-EM[Oxidation]EVTSESPEK') ==
            ('EMEVTSESPEK', {'N-term': 42.010565, 1: 15.994915}))
    assert (proforma.parse(
                '[Phospho]^2[Methyl]?[Acetyl]-EM[Oxidation]EVTSESPEK') ==
            ('EMEVTSESPEK', {'N-term': 42.010565, 1: 15.994915}))
    assert (proforma.parse('[Phospho]^2?[Acetyl]-EM[Oxidation]EVTSESPEK') ==
            ('EMEVTSESPEK', {'N-term': 42.010565, 1: 15.994915}))
    # Multiple possible modification positions.
    assert (proforma.parse('EM[Oxidation]EVT[#g1]S[#g1]ES[Phospho#g1]PEK') ==
            ('EMEVTSESPEK', {1: 15.994915, 7: 79.966331}))
    with pytest.raises(ValueError):
        proforma.parse('EM[Oxidation]EVT[#g1]S[Phospho#g1]ES[Phospho#g1]PEK')
    # Ranges of modification positions.
    with pytest.raises(NotImplementedError):
        proforma.parse('PROT(EOSFORMS)[+19.0523]ISK')
    with pytest.raises(NotImplementedError):
        proforma.parse('PROT(EOC[Carbamidomethyl]FORMS)[+19.0523]ISK')
    with pytest.raises(NotImplementedError):
        proforma.parse('P(ROT(EOSFORMS)[+19.0523]IS)[+19.0523]K')
    # Modification localization score.
    assert (proforma.parse('EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]E'
                           'S[Phospho#g1(0.90)]PEK') ==
            ('EMEVTSESPEK', {1: 15.994915, 7: 79.966331}))
    assert (proforma.parse('[Phospho#s1]?EM[Oxidation]EVT[#s1(0.01)]'
                           'S[#s1(0.09)]ES[#s1(0.90)]PEK') ==
            ('EMEVTSESPEK', {1: 15.994915}))


def test_proforma_global_modification():
    assert (proforma.parse('<[S-carboxamidomethyl-L-cysteine]@C>'
                           'ATPEILTCNSIGCLK') ==
            ('ATPEILTCNSIGCLK', {7: 57.021464, 12: 57.021464}))
    assert (proforma.parse('<[MOD:01090]@C>ATPEILTCNSIGCLK') ==
            ('ATPEILTCNSIGCLK', {7: 57.021464, 12: 57.021464}))
    assert (proforma.parse('<[Oxidation]@C,M>MTPEILTCNSIGCLK') ==
            ('MTPEILTCNSIGCLK', {0: 15.994915, 7: 15.994915, 12: 15.994915}))


def test_proforma_info():
    assert proforma.parse('ELV[INFO:AnyString]IS') == ('ELVIS', {})
    assert proforma.parse('ELV[info:AnyString]IS') == ('ELVIS', {})
    assert (proforma.parse('ELVIS[Phospho|INFO:newly discovered]K')
            == ('ELVISK', {4: 79.966331}))
    assert (proforma.parse('ELVIS[Phospho|INFO:newly discovered|'
                           'INFO:really awesome]K') ==
            ('ELVISK', {4: 79.966331}))


def test_proforma_pipe():
    assert (proforma.parse('ELVIS[U:Phospho|+79.966331]K') ==
            ('ELVISK', {4: 79.966331}))
    assert (proforma.parse('ELVIS[U:Phospho|Obs:+79.978]K') ==
            ('ELVISK', {4: 79.966331}))
    assert (proforma.parse('ELVIS[Phospho|O-phospho-L-serine]K') ==
            ('ELVISK', {4: 79.966331}))
    assert (proforma.parse('ELVIS[UNIMOD:21|MOD:00046]K') ==
            ('ELVISK', {4: 79.966331}))
    assert (proforma.parse('ELVIS[UNIMOD:21|Phospho]K') ==
            ('ELVISK', {4: 79.966331}))
    assert (proforma.parse('ELVIS[Phospho|O-phospho-L-serine|Obs:+79.966]K') ==
            ('ELVISK', {4: 79.966331}))
    assert (proforma.parse('ELVIS[Obs:+79.966|Phospho|Sulfo]K') ==
            ('ELVISK', {4: 79.966}))


def test_proforma_cache():
    cache_dir = '.cache_test'
    shutil.rmtree(cache_dir, ignore_errors=True)
    # Disable cache.
    proforma.cache_dir = None
    assert not os.path.isfile(os.path.join(cache_dir, 'UNIMOD.pkl'))
    assert (proforma.parse('EM[U:Oxidation]EVEES[U:Phospho]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 79.966331}))
    assert not os.path.isfile(os.path.join(cache_dir, 'UNIMOD.pkl'))
    assert (proforma.parse('EM[U:Oxidation]EVEES[U:Phospho]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 79.966331}))
    assert not os.path.isfile(os.path.join(cache_dir, 'UNIMOD.pkl'))
    # Enable cache.
    proforma.cache_dir = cache_dir
    assert not os.path.isfile(os.path.join(cache_dir, 'UNIMOD.pkl'))
    assert (proforma.parse('EM[U:Oxidation]EVEES[U:Phospho]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 79.966331}))
    assert os.path.isfile(os.path.join(cache_dir, 'UNIMOD.pkl'))
    assert (proforma.parse('EM[U:Oxidation]EVEES[U:Phospho]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 79.966331}))
    assert os.path.isfile(os.path.join(cache_dir, 'UNIMOD.pkl'))
    # Clear cache.
    proforma.clear_cache()
    assert not os.path.isfile(os.path.join(cache_dir, 'UNIMOD.pkl'))
    assert (proforma.parse('EM[U:Oxidation]EVEES[U:Phospho]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 79.966331}))
    assert os.path.isfile(os.path.join(cache_dir, 'UNIMOD.pkl'))
    # Clean-up.
    shutil.rmtree(cache_dir, ignore_errors=True)
