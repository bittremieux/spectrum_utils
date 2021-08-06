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
    # FIXME
    # assert (sus._parse_proforma('EMEVTK[X:DSS#XL1]SESPEK') ==
    #         ('EMEVTKSESPEK', {5: 138.06807961}))
    # GNO (mandatory prefix).
    assert (sus._parse_proforma('NEEYN[G:G59626AS]K') ==
            ('NEEYNK', {4: 1931.69}))
