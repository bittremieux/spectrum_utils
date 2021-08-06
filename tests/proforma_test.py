import math

import pytest

import spectrum_utils.spectrum as sus


def test_proforma_canonical():
    assert sus._parse_proforma('PEPTIDE') == ('PEPTIDE', {})
    assert sus._parse_proforma('ACDEFGHIKL') == ('ACDEFGHIKL', {})
    assert sus._parse_proforma('MNPQRSTVWY') == ('MNPQRSTVWY', {})
    assert sus._parse_proforma('PEPTJDE') == ('PEPTJDE', {})
    assert sus._parse_proforma('PEPTIOE') == ('PEPTIOE', {})
    assert sus._parse_proforma('PEUTIDE') == ('PEUTIDE', {})
    assert sus._parse_proforma('PEPTIXDE') == ('PEPTIXDE', {})


def test_proforma_name():
    # Unimod without prefix.
    assert (sus._parse_proforma('EM[Oxidation]EVEES[Phospho]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 79.966331}))
    # Unimod with prefix.
    assert (sus._parse_proforma('EM[U:Oxidation]EVEES[U:Phospho]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 79.966331}))
    # Unimod with internal brackets in the name.
    assert (sus._parse_proforma('EM[Oxidation]EVE[Cation:Mg[II]]E'
                                'S[Phospho]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 4: 21.969392, 6: 79.966331}))
    # PSI-MOD without prefix.
    assert (sus._parse_proforma('EM[L-methionine sulfoxide]EVEE'
                                'S[O-phospho-L-serine]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 79.966331}))
    # PSI-MOD with prefix.
    assert (sus._parse_proforma('EM[M:L-methionine sulfoxide]EVEE'
                                'S[M:O-phospho-L-serine]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 79.966331}))
    # RESID (mandatory prefix).
    assert (sus._parse_proforma('EM[R:L-methionine sulfone]EVEE'
                                'S[R:O-phospho-L-serine]PEK') ==
            ('EMEVEESPEK', {1: 31.989829, 6: 79.966331}))
    # XL-MOD (mandatory prefix).
    assert (sus._parse_proforma('EMEVTK[X:DSS#XL1]SESPEK') ==
            ('EMEVTKSESPEK', {5: 138.06807961}))
    # GNO (mandatory prefix).
    assert (sus._parse_proforma('NEEYN[G:G59626AS]K') ==
            ('NEEYNK', {4: 1931.69}))


def test_proforma_accession():
    # Valid accessions.
    assert (sus._parse_proforma('EM[UNIMOD:35]EVEES[UNIMOD:56]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 45.029395}))
    assert (sus._parse_proforma('EM[MOD:00719]EVEES[MOD:00046]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 79.966331}))
    assert (sus._parse_proforma('EM[RESID:AA0581]EVEES[RESID:AA0037]PEK') ==
            ('EMEVEESPEK', {1: 15.994915, 6: 79.966331}))
    assert (sus._parse_proforma('NEEYN[GNO:G59626AS]K') ==
            ('NEEYNK', {4: 1931.69}))
    assert (sus._parse_proforma('YPVLN[GNO:G62765YT]VTMP'
                                'N[GNO:G02815KT]NSNGKFDK') ==
            ('YPVLNVTMPNNSNGKFDK', {4: 1720.59, 9: 1234.43}))
    # Invalid accessions.
    with pytest.raises(KeyError):
        sus._parse_proforma('EM[U:35]EVEES[U:56]PEK')
    with pytest.raises(KeyError):
        sus._parse_proforma('EM[M:00719]EVEES[M:00046]PEK')
    with pytest.raises(KeyError):
        sus._parse_proforma('EM[R:AA0581]EVEES[R:AA0037]PEK')


def test_proforma_xlink():
    # Crosslink notation.
    assert (sus._parse_proforma('EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1]') ==
            ('EMEVTKSESPEK', {5: 138.06807961}))
    assert (sus._parse_proforma('EMK[XLMOD:02000#XL1]EVTKSE[XLMOD:02010#XL2]S'
                                'K[#XL1]PEK[#XL2]AR') ==
            ('EMKEVTKSESKPEKAR', {2: 138.06807961, 8: -18.01056027}))
    assert (sus._parse_proforma('EMEVTK[XLMOD:02001#XL1]SESPEK') ==
            ('EMEVTKSESPEK', {5: 138.06807961}))
    assert (sus._parse_proforma('EMEVTK[XLMOD:02001]SESPEK') ==
            ('EMEVTKSESPEK', {5: 138.06807961}))
    # Inter-chain crosslinks.
    with pytest.raises(NotImplementedError):
        sus._parse_proforma('SEK[XLMOD:02001#XL1]UENCE\\EMEVT'
                            'K[XLMOD:02001#XL1]SESPEK')
    with pytest.raises(NotImplementedError):
        sus._parse_proforma('SEK[XLMOD:02001#XL1]UENCE\\EMEVTK[#XL1]SESPEK')
    # Disulfide linkages.
    assert (sus._parse_proforma('EVTSEKC[MOD:00034#XL1]LEMSC[#XL1]EFD') ==
            ('EVTSEKCLEMSCEFD', {6: -2.015650}))
    assert (sus._parse_proforma('EVTSEKC[L-cystine (cross-link)#XL1]LEMS'
                                'C[#XL1]EFD') ==
            ('EVTSEKCLEMSCEFD', {6: -2.015650}))
    with pytest.raises(NotImplementedError):
        sus._parse_proforma('FVNQHLC[MOD:00034#XL1]GSHLVEALYLV'
                            'C[MOD:00034#XL2]GERGFFYTPKA\\GIVEQ'
                            'C[MOD:00034#XL3]C[#XL1]TSIC[#XL3]SLYQLENY'
                            'C[#XL2]N')
    assert (sus._parse_proforma('EVTSEKC[XLMOD:02009#XL1]LEMSC[#XL1]EFD') ==
            ('EVTSEKCLEMSCEFD', {6: -2.01565007}))
    assert (sus._parse_proforma('EVTSEKC[X:Disulfide#XL1]LEMSC[#XL1]EFD') ==
            ('EVTSEKCLEMSCEFD', {6: -2.01565007}))
    assert (sus._parse_proforma('EVTSEKC[half cystine]LEMSC[half cystine]EFD')
            == ('EVTSEKCLEMSCEFD', {6: -1.007825, 11: -1.007825}))
    assert (sus._parse_proforma('EVTSEKC[MOD:00798]LEMSC[MOD:00798]EFDEVTSEK'
                                'C[MOD:00798]LEMSC[MOD:00798]EFD') ==
            ('EVTSEKCLEMSCEFDEVTSEKCLEMSCEFD', {6: -1.007825, 11: -1.007825,
                                                21: -1.007825, 26: -1.007825}))


def test_proforma_delta_mass():
    # No prefix.
    assert (sus._parse_proforma('EM[+15.9949]EVEES[+79.9663]PEK') ==
            ('EMEVEESPEK', {1: 15.9949, 6: 79.9663}))
    assert (sus._parse_proforma('EM[+15.995]EVEES[-18.01]PEK') ==
            ('EMEVEESPEK', {1: 15.995, 6: -18.01}))
    # CV prefix.
    assert (sus._parse_proforma('EM[U:+15.9949]EVEES[U:+79.9663]PEK') ==
            ('EMEVEESPEK', {1: 15.9949, 6: 79.9663}))
    assert (sus._parse_proforma('EM[U:+15.995]EVEES[U:+79.966]PEK') ==
            ('EMEVEESPEK', {1: 15.995, 6: 79.966}))
    # Experimentally observed prefix.
    assert (sus._parse_proforma('EM[U:+15.995]EVEES[Obs:+79.978]PEK') ==
            ('EMEVEESPEK', {1: 15.995, 6: 79.978}))


def test_proforma_gap():
    assert (sus._parse_proforma('RTAAX[+367.0537]WT') ==
            ('RTAAXWT', {4: 367.0537}))


def test_proforma_formula():
    sequence, mass_diffs = sus._parse_proforma('SEQUEN[Formula:C12H20O2]CE')
    assert sequence == 'SEQUENCE'
    assert len(mass_diffs) == 1
    assert math.isclose(mass_diffs[5], 196.1463, abs_tol=0.00005)
    sequence, mass_diffs = sus._parse_proforma('SEQUEN[Formula:C12 H20 O2]CE')
    assert sequence == 'SEQUENCE'
    assert len(mass_diffs) == 1
    assert math.isclose(mass_diffs[5], 196.1463, abs_tol=0.00005)


def test_proforma_glycan():
    assert (sus._parse_proforma('SEQUEN[Glycan:HexNAc1Hex2]CE') ==
            ('SEQUENCE', {5: 527.185019356}))
    assert (sus._parse_proforma('SEQUEN[Glycan:HexNAc1 Hex2]CE') ==
            ('SEQUENCE', {5: 527.185019356}))


def test_proforma_special():
    # N-terminal and C-terminal modifications.
    assert (sus._parse_proforma('[iTRAQ4plex]-EM[Oxidation]EVNE'
                                'S[Phospho]PEK') ==
            ('EMEVNESPEK', {'N-term': 144.102063, 1: 15.994915, 6: 79.966331}))
    assert (sus._parse_proforma('[iTRAQ4plex]-EM[U:Oxidation]EVNE'
                                'S[Phospho]PEK[iTRAQ4plex]-[Methyl]') ==
            ('EMEVNESPEK', {'N-term': 144.102063, 1: 15.994915, 6: 79.966331,
                            9: 144.102063, 'C-term': 14.01565}))
    # Labile modifications.
    assert (sus._parse_proforma('{Glycan:Hex}EM[U:Oxidation]EVNES[Phospho]PE'
                                'K[iTRAQ4plex]') ==
            ('EMEVNESPEK', {1: 15.994915, 6: 79.966331, 9: 144.102063}))
    assert (sus._parse_proforma('{Glycan:Hex}[iTRAQ4plex]-EM[Oxidation]EVNE'
                                'S[Phospho]PEK[iTRAQ4plex]') ==
            ('EMEVNESPEK', {'N-term': 144.102063, 1: 15.994915, 6: 79.966331,
                            9: 144.102063}))
    assert (sus._parse_proforma('{Glycan:Hex}[iTRAQ4plex]-EM[Oxidation]EVNE'
                                'S[Phospho]PEK[iTRAQ4plex]-[Methyl]') ==
            ('EMEVNESPEK', {'N-term': 144.102063, 1: 15.994915, 6: 79.966331,
                            9: 144.102063, 'C-term': 14.01565}))
    assert (sus._parse_proforma('{Glycan:Hex}{Glycan:NeuAc}EMEVNESPEK') ==
            ('EMEVNESPEK', {}))


def test_proforma_ambiguous_position():
    # Unknown modification position.
    assert (sus._parse_proforma('[Phospho]?EM[Oxidation]EVTSESPEK') ==
            ('EMEVTSESPEK', {1: 15.994915}))
    assert (sus._parse_proforma('[Phospho][Phospho]?[Acetyl]-EM[Oxidation]'
                                'EVTSESPEK') ==
            ('EMEVTSESPEK', {'N-term': 42.010565, 1: 15.994915}))
    assert (sus._parse_proforma('[Phospho]^2[Methyl]?[Acetyl]-EM[Oxidation]'
                                'EVTSESPEK') ==
            ('EMEVTSESPEK', {'N-term': 42.010565, 1: 15.994915}))
    assert (sus._parse_proforma('[Phospho]^2?[Acetyl]-EM[Oxidation]EVTSESPEK')
            == ('EMEVTSESPEK', {'N-term': 42.010565, 1: 15.994915}))
    # Multiple possible modification positions.
    assert (sus._parse_proforma('EM[Oxidation]EVT[#g1]S[#g1]ES[Phospho#g1]PEK')
            == ('EMEVTSESPEK', {1: 15.994915, 7: 79.966331}))
    with pytest.raises(ValueError):
        sus._parse_proforma('EM[Oxidation]EVT[#g1]S[Phospho#g1]E'
                            'S[Phospho#g1]PEK')
    # Ranges of modification positions.
    with pytest.raises(NotImplementedError):
        sus._parse_proforma('PROT(EOSFORMS)[+19.0523]ISK')
    with pytest.raises(NotImplementedError):
        sus._parse_proforma('PROT(EOC[Carbamidomethyl]FORMS)[+19.0523]ISK')
    with pytest.raises(NotImplementedError):
        sus._parse_proforma('P(ROT(EOSFORMS)[+19.0523]IS)[+19.0523]K')
    # Modification localization score.
    assert (sus._parse_proforma('EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]E'
                                'S[Phospho#g1(0.90)]PEK') ==
            ('EMEVTSESPEK', {1: 15.994915, 7: 79.966331}))
    assert (sus._parse_proforma('[Phospho#s1]?EM[Oxidation]EVT[#s1(0.01)]'
                                'S[#s1(0.09)]ES[#s1(0.90)]PEK') ==
            ('EMEVTSESPEK', {1: 15.994915}))


def test_proforma_global_modification():
    assert (sus._parse_proforma('<[S-carboxamidomethyl-L-cysteine]@C>'
                                'ATPEILTCNSIGCLK') ==
            ('ATPEILTCNSIGCLK', {7: 57.021464, 12: 57.021464}))
    assert (sus._parse_proforma('<[MOD:01090]@C>ATPEILTCNSIGCLK') ==
            ('ATPEILTCNSIGCLK', {7: 57.021464, 12: 57.021464}))
    assert (sus._parse_proforma('<[Oxidation]@C,M>MTPEILTCNSIGCLK') ==
            ('MTPEILTCNSIGCLK', {0: 15.994915, 7: 15.994915, 12: 15.994915}))


def test_proforma_info():
    assert sus._parse_proforma('ELV[INFO:AnyString]IS') == ('ELVIS', {})
    assert sus._parse_proforma('ELV[info:AnyString]IS') == ('ELVIS', {})
    assert (sus._parse_proforma('ELVIS[Phospho|INFO:newly discovered]K') ==
            ('ELVISK', {4: 79.966331}))
    assert (sus._parse_proforma('ELVIS[Phospho|INFO:newly discovered|'
                                'INFO:really awesome]K') ==
            ('ELVISK', {4: 79.966331}))


def test_proforma_pipe():
    assert (sus._parse_proforma('ELVIS[U:Phospho|+79.966331]K') ==
            ('ELVISK', {4: 79.966331}))
    assert (sus._parse_proforma('ELVIS[U:Phospho|Obs:+79.978]K') ==
            ('ELVISK', {4: 79.966331}))
    assert (sus._parse_proforma('ELVIS[Phospho|O-phospho-L-serine]K') ==
            ('ELVISK', {4: 79.966331}))
    assert (sus._parse_proforma('ELVIS[UNIMOD:21|MOD:00046]K') ==
            ('ELVISK', {4: 79.966331}))
    assert (sus._parse_proforma('ELVIS[UNIMOD:21|Phospho]K') ==
            ('ELVISK', {4: 79.966331}))
    assert (sus._parse_proforma('ELVIS[Phospho|O-phospho-L-serine|'
                                'Obs:+79.966]K') ==
            ('ELVISK', {4: 79.966331}))
    assert (sus._parse_proforma('ELVIS[Obs:+79.966|Phospho|Sulfo]K') ==
            ('ELVISK', {4: 79.966}))
