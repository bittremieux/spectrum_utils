import pytest

from spectrum_utils import proforma


def test_proforma_canonical():
    assert proforma.parse("PEPTIDE") == [
        proforma.Proteoform(sequence="PEPTIDE", modifications=[])
    ]
    assert proforma.parse("ACDEFGHIKL") == [
        proforma.Proteoform(sequence="ACDEFGHIKL", modifications=[])
    ]
    assert proforma.parse("MNPQRSTVWY") == [
        proforma.Proteoform(sequence="MNPQRSTVWY", modifications=[])
    ]
    # Ambiguous amino acids.
    assert proforma.parse("PEPTJDE") == [
        proforma.Proteoform(sequence="PEPTJDE", modifications=[])
    ]
    assert proforma.parse("PEPTIOE") == [
        proforma.Proteoform(sequence="PEPTIOE", modifications=[])
    ]
    assert proforma.parse("PEUTIDE") == [
        proforma.Proteoform(sequence="PEUTIDE", modifications=[])
    ]
    assert proforma.parse("PEPTIXDE") == [
        proforma.Proteoform(sequence="PEPTIXDE", modifications=[])
    ]
    assert proforma.parse("BEPTIDE") == [
        proforma.Proteoform(sequence="BEPTIDE", modifications=[])
    ]
    assert proforma.parse("PEPTIZE") == [
        proforma.Proteoform(sequence="PEPTIZE", modifications=[])
    ]


def test_proforma_name():
    # Unimod without prefix.
    assert proforma.parse("EM[Oxidation]EVEES[Phospho]PEK") == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Oxidation"
                        )
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("EM[Oxidation]EVEES[Phospho]PEK", True) == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.966331,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                ),
            ],
        )
    ]
    # Unimod with prefix.
    assert proforma.parse("EM[U:Oxidation]EVEES[U:Phospho]PEK") == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD", name="Oxidation"
                        )
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD", name="Phospho"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("EM[U:Oxidation]EVEES[U:Phospho]PEK", True) == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.966331,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                ),
            ],
        )
    ]
    # Unimod with internal brackets in the name.
    assert proforma.parse("EM[Oxidation]EVE[Cation:Mg[II]]ES[Phospho]PEK") == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Oxidation"
                        )
                    ],
                ),
                proforma.Modification(
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Cation:Mg[II]"
                        )
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "EM[Oxidation]EVE[Cation:Mg[II]]ES[Phospho]PEK", True
    ) == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=21.969392,
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:956",
                            name="Cation:Mg[II]",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.966331,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                ),
            ],
        )
    ]
    # PSI-MOD without prefix.
    assert proforma.parse(
        "EM[L-methionine sulfoxide]EVEE" "S[O-phospho-L-serine]PEK"
    ) == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None,
                            name="L-methionine sulfoxide",
                        )
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None,
                            name="O-phospho-L-serine",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "EM[L-methionine sulfoxide]EVEE" "S[O-phospho-L-serine]PEK", True
    ) == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00719",
                            name="L-methionine sulfoxide",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.966331,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00046",
                            name="O-phospho-L-serine",
                        )
                    ],
                ),
            ],
        )
    ]
    # PSI-MOD with prefix.
    assert proforma.parse(
        "EM[M:L-methionine sulfoxide]EVEE" "S[M:O-phospho-L-serine]PEK"
    ) == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            name="L-methionine sulfoxide",
                        )
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            name="O-phospho-L-serine",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "EM[L-methionine sulfoxide]EVEE" "S[O-phospho-L-serine]PEK", True
    ) == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00719",
                            name="L-methionine sulfoxide",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.966331,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00046",
                            name="O-phospho-L-serine",
                        )
                    ],
                ),
            ],
        )
    ]
    # RESID (mandatory prefix).
    assert proforma.parse(
        "EM[R:L-methionine sulfone]EVEE" "S[R:O-phospho-L-serine]PEK"
    ) == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="RESID",
                            name="L-methionine sulfone",
                        )
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="RESID",
                            name="O-phospho-L-serine",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "EM[R:L-methionine sulfone]EVEE" "S[R:O-phospho-L-serine]PEK", True
    ) == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    mass=31.989829,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="RESID",
                            accession="RESID:AA0251",
                            name="L-methionine sulfone",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.966331,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="RESID",
                            accession="RESID:AA0037",
                            name="O-phospho-L-serine",
                        )
                    ],
                ),
            ],
        )
    ]
    # XL-MOD (mandatory prefix).
    assert proforma.parse("EMEVTK[X:DSS#XL1]SESPEK") == [
        proforma.Proteoform(
            sequence="EMEVTKSESPEK",
            modifications=[
                proforma.Modification(
                    position=5,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD", name="DSS"
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse("EMEVTK[X:DSS#XL1]SESPEK", True) == [
        proforma.Proteoform(
            sequence="EMEVTKSESPEK",
            modifications=[
                proforma.Modification(
                    mass=138.06807961,
                    position=5,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02001",
                            name="DSS",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        )
    ]
    # GNO (mandatory prefix).
    assert proforma.parse("NEEYN[G:G59626AS]K") == [
        proforma.Proteoform(
            sequence="NEEYNK",
            modifications=[
                proforma.Modification(
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="GNO", name="G59626AS"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("NEEYN[G:G59626AS]K", True) == [
        proforma.Proteoform(
            sequence="NEEYNK",
            modifications=[
                proforma.Modification(
                    mass=1931.69,
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="GNO",
                            accession="GNO:G59626AS",
                            name="G59626AS",
                        )
                    ],
                ),
            ],
        )
    ]
    # Resolving incorrect terms.
    # Non-existing terms without prefix.
    with pytest.raises(KeyError):
        proforma.parse("EM[WeirdMod]EVEES[Phospho]PEK", True)
    # Non-existing terms with prefix.
    with pytest.raises(KeyError):
        proforma.parse("EM[U:WeirdMod]EVEES[U:Phospho]PEK", True)
    # Non-existing prefix (will be considered as no prefix).
    with pytest.raises(KeyError):
        proforma.parse("EM[RandomPrefix:Oxidation]EVEES[Phospho]PEK", True)


def test_proforma_accession():
    # Unimod.
    assert proforma.parse("EM[UNIMOD:35]EVEES[UNIMOD:56]PEK") == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                        )
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:56",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("EM[UNIMOD:35]EVEES[UNIMOD:56]PEK", True) == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=45.029395,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:56",
                            name="Acetyl:2H(3)",
                        )
                    ],
                ),
            ],
        )
    ]
    # PSI-MOD.
    assert proforma.parse("EM[MOD:00719]EVEES[MOD:00046]PEK") == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD", accession="MOD:00719"
                        )
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD", accession="MOD:00046"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("EM[MOD:00719]EVEES[MOD:00046]PEK", True) == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00719",
                            name="L-methionine sulfoxide",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.966331,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00046",
                            name="O-phospho-L-serine",
                        )
                    ],
                ),
            ],
        )
    ]
    # RESID.
    assert proforma.parse("EM[RESID:AA0581]EVEES[RESID:AA0037]PEK") == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="RESID",
                            accession="RESID:AA0581",
                        )
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="RESID",
                            accession="RESID:AA0037",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("EM[RESID:AA0581]EVEES[RESID:AA0037]PEK", True) == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="RESID",
                            accession="RESID:AA0581",
                            name="L-methionine (R)-sulfoxide",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.966331,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="RESID",
                            accession="RESID:AA0037",
                            name="O-phospho-L-serine",
                        )
                    ],
                ),
            ],
        )
    ]
    # GNO.
    assert proforma.parse("NEEYN[GNO:G59626AS]K") == [
        proforma.Proteoform(
            sequence="NEEYNK",
            modifications=[
                proforma.Modification(
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="GNO",
                            accession="GNO:G59626AS",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("NEEYN[GNO:G59626AS]K", True) == [
        proforma.Proteoform(
            sequence="NEEYNK",
            modifications=[
                proforma.Modification(
                    mass=1931.69,
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="GNO",
                            accession="GNO:G59626AS",
                            name="G59626AS",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "YPVLN[GNO:G62765YT]VTMP" "N[GNO:G02815KT]NSNGKFDK"
    ) == [
        proforma.Proteoform(
            sequence="YPVLNVTMPNNSNGKFDK",
            modifications=[
                proforma.Modification(
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="GNO",
                            accession="GNO:G62765YT",
                        )
                    ],
                ),
                proforma.Modification(
                    position=9,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="GNO",
                            accession="GNO:G02815KT",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "YPVLN[GNO:G62765YT]VTMP" "N[GNO:G02815KT]NSNGKFDK", True
    ) == [
        proforma.Proteoform(
            sequence="YPVLNVTMPNNSNGKFDK",
            modifications=[
                proforma.Modification(
                    mass=1720.59,
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="GNO",
                            accession="GNO:G62765YT",
                            name="G62765YT",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=1234.43,
                    position=9,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="GNO",
                            accession="GNO:G02815KT",
                            name="G02815KT",
                        )
                    ],
                ),
            ],
        )
    ]
    # Invalid accessions.
    with pytest.raises(KeyError):
        proforma.parse("EM[UNIMOD:one]EVEES[UNIMOD:two]PEK", True)
    with pytest.raises(KeyError):
        proforma.parse("EM[MOD:one]EVEES[MOD:two]PEK", True)
    with pytest.raises(KeyError):
        proforma.parse("EM[RESID:one]EVEES[RESID:two]PEK", True)
    with pytest.raises(KeyError):
        proforma.parse("YPVLN[GNO:one]VTMPN[GNO:two]NSNGKFDK", True)


def test_proforma_xlink():
    # Crosslink notation.
    assert proforma.parse("EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1]") == [
        proforma.Proteoform(
            sequence="EMEVTKSESPEK",
            modifications=[
                proforma.Modification(
                    position=5,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02001",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    position=11,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse("EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1]", True) == [
        proforma.Proteoform(
            sequence="EMEVTKSESPEK",
            modifications=[
                proforma.Modification(
                    mass=138.06807961,
                    position=5,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02001",
                            name="DSS",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    position=11,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse(
        "EMK[XLMOD:02000#XL1]EVTKSE[XLMOD:02010#XL2]S" "K[#XL1]PEK[#XL2]AR"
    ) == [
        proforma.Proteoform(
            sequence="EMKEVTKSESKPEKAR",
            modifications=[
                proforma.Modification(
                    position=2,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02000",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    position=8,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02010",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL2"
                    ),
                ),
                proforma.Modification(
                    position=10,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    position=13,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL2"
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse(
        "EMK[XLMOD:02000#XL1]EVTKSE[XLMOD:02010#XL2]S" "K[#XL1]PEK[#XL2]AR",
        True,
    ) == [
        proforma.Proteoform(
            sequence="EMKEVTKSESKPEKAR",
            modifications=[
                proforma.Modification(
                    mass=138.06807961,
                    position=2,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02000",
                            name="BS3",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    mass=-18.01056027,
                    position=8,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02010",
                            name="1-ethyl-3-(3-Dimethylaminopropyl)"
                            "carbodiimide hydrochloride",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL2"
                    ),
                ),
                proforma.Modification(
                    position=10,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    position=13,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL2"
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse("EMEVTK[XLMOD:02001#XL1]SESPEK") == [
        proforma.Proteoform(
            sequence="EMEVTKSESPEK",
            modifications=[
                proforma.Modification(
                    position=5,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02001",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse("EMEVTK[XLMOD:02001#XL1]SESPEK", True) == [
        proforma.Proteoform(
            sequence="EMEVTKSESPEK",
            modifications=[
                proforma.Modification(
                    mass=138.06807961,
                    position=5,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02001",
                            name="DSS",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse("EMEVTK[XLMOD:02001]SESPEK") == [
        proforma.Proteoform(
            sequence="EMEVTKSESPEK",
            modifications=[
                proforma.Modification(
                    position=5,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02001",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("EMEVTK[XLMOD:02001]SESPEK", True) == [
        proforma.Proteoform(
            sequence="EMEVTKSESPEK",
            modifications=[
                proforma.Modification(
                    mass=138.06807961,
                    position=5,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02001",
                            name="DSS",
                        )
                    ],
                ),
            ],
        )
    ]
    # Inter-chain crosslinks.
    assert proforma.parse(
        "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[XLMOD:02001#XL1]SESPEK"
    ) == [
        proforma.Proteoform(
            sequence="SEKUENCE",
            modifications=[
                proforma.Modification(
                    position=2,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02001",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        ),
        proforma.Proteoform(
            sequence="EMEVTKSESPEK",
            modifications=[
                proforma.Modification(
                    position=5,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02001",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        ),
    ]
    assert proforma.parse(
        "SEK[XLMOD:02001#XL1]UENCE//" "EMEVTK[XLMOD:02001#XL1]SESPEK", True
    ) == [
        proforma.Proteoform(
            sequence="SEKUENCE",
            modifications=[
                proforma.Modification(
                    mass=138.06807961,
                    position=2,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02001",
                            name="DSS",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        ),
        proforma.Proteoform(
            sequence="EMEVTKSESPEK",
            modifications=[
                proforma.Modification(
                    mass=138.06807961,
                    position=5,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02001",
                            name="DSS",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        ),
    ]
    assert proforma.parse("SEK[XLMOD:02001#XL1]UENCE//EMEVTK[#XL1]SESPEK") == [
        proforma.Proteoform(
            sequence="SEKUENCE",
            modifications=[
                proforma.Modification(
                    position=2,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02001",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        ),
        proforma.Proteoform(
            sequence="EMEVTKSESPEK",
            modifications=[
                proforma.Modification(
                    position=5,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        ),
    ]
    assert proforma.parse(
        "SEK[XLMOD:02001#XL1]UENCE//" "EMEVTK[XLMOD:02001#XL1]SESPEK", True
    ) == [
        proforma.Proteoform(
            sequence="SEKUENCE",
            modifications=[
                proforma.Modification(
                    mass=138.06807961,
                    position=2,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02001",
                            name="DSS",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        ),
        proforma.Proteoform(
            sequence="EMEVTKSESPEK",
            modifications=[
                proforma.Modification(
                    mass=138.06807961,
                    position=5,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02001",
                            name="DSS",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        ),
    ]
    # Disulfide linkages.
    assert proforma.parse("EVTSEKC[MOD:00034#XL1]LEMSC[#XL1]EFD") == [
        proforma.Proteoform(
            sequence="EVTSEKCLEMSCEFD",
            modifications=[
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD", accession="MOD:00034"
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    position=11,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse("EVTSEKC[MOD:00034#XL1]LEMSC[#XL1]EFD", True) == [
        proforma.Proteoform(
            sequence="EVTSEKCLEMSCEFD",
            modifications=[
                proforma.Modification(
                    mass=-2.015650,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00034",
                            name="L-cystine (cross-link)",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    position=11,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse(
        "EVTSEKC[L-cystine (cross-link)#XL1]LEMS" "C[#XL1]EFD"
    ) == [
        proforma.Proteoform(
            sequence="EVTSEKCLEMSCEFD",
            modifications=[
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None,
                            name="L-cystine (cross-link)",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    position=11,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse(
        "EVTSEKC[L-cystine (cross-link)#XL1]LEMS" "C[#XL1]EFD", True
    ) == [
        proforma.Proteoform(
            sequence="EVTSEKCLEMSCEFD",
            modifications=[
                proforma.Modification(
                    mass=-2.015650,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00034",
                            name="L-cystine (cross-link)",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    position=11,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse(
        "FVNQHLC[MOD:00034#XL1]GSHLVEALYLVC[MOD:00034#XL2]G"
        "ERGFFYTPKA//"
        "GIVEQC[MOD:00034#XL3]C[#XL1]TSIC[#XL3]SLYQLENY"
        "C[#XL2]N"
    ) == [
        proforma.Proteoform(
            sequence="FVNQHLCGSHLVEALYLVCGERGFFYTPKA",
            modifications=[
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD", accession="MOD:00034"
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    position=18,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD", accession="MOD:00034"
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL2"
                    ),
                ),
            ],
        ),
        proforma.Proteoform(
            sequence="GIVEQCCTSICSLYQLENYCN",
            modifications=[
                proforma.Modification(
                    position=5,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD", accession="MOD:00034"
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL3"
                    ),
                ),
                proforma.Modification(
                    position=6,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    position=10,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL3"
                    ),
                ),
                proforma.Modification(
                    position=19,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL2"
                    ),
                ),
            ],
        ),
    ]
    assert proforma.parse(
        "FVNQHLC[MOD:00034#XL1]GSHLVEALYLVC[MOD:00034#XL2]G"
        "ERGFFYTPKA//"
        "GIVEQC[MOD:00034#XL3]C[#XL1]TSIC[#XL3]SLYQLENY"
        "C[#XL2]N",
        True,
    ) == [
        proforma.Proteoform(
            sequence="FVNQHLCGSHLVEALYLVCGERGFFYTPKA",
            modifications=[
                proforma.Modification(
                    mass=-2.015650,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00034",
                            name="L-cystine (cross-link)",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    mass=-2.015650,
                    position=18,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00034",
                            name="L-cystine (cross-link)",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL2"
                    ),
                ),
            ],
        ),
        proforma.Proteoform(
            sequence="GIVEQCCTSICSLYQLENYCN",
            modifications=[
                proforma.Modification(
                    mass=-2.015650,
                    position=5,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00034",
                            name="L-cystine (cross-link)",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL3"
                    ),
                ),
                proforma.Modification(
                    position=6,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    position=10,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL3"
                    ),
                ),
                proforma.Modification(
                    position=19,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL2"
                    ),
                ),
            ],
        ),
    ]
    assert proforma.parse("EVTSEKC[XLMOD:02009#XL1]LEMSC[#XL1]EFD") == [
        proforma.Proteoform(
            sequence="EVTSEKCLEMSCEFD",
            modifications=[
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02009",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    position=11,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse("EVTSEKC[XLMOD:02009#XL1]LEMSC[#XL1]EFD", True) == [
        proforma.Proteoform(
            sequence="EVTSEKCLEMSCEFD",
            modifications=[
                proforma.Modification(
                    mass=-2.01565007,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02009",
                            name="Disulfide",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    position=11,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse("EVTSEKC[X:Disulfide#XL1]LEMSC[#XL1]EFD") == [
        proforma.Proteoform(
            sequence="EVTSEKCLEMSCEFD",
            modifications=[
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD", name="Disulfide"
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    position=11,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse("EVTSEKC[X:Disulfide#XL1]LEMSC[#XL1]EFD", True) == [
        proforma.Proteoform(
            sequence="EVTSEKCLEMSCEFD",
            modifications=[
                proforma.Modification(
                    mass=-2.01565007,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="XLMOD",
                            accession="XLMOD:02009",
                            name="Disulfide",
                        )
                    ],
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
                proforma.Modification(
                    position=11,
                    label=proforma.Label(
                        type=proforma.LabelType.XL, label="XL1"
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse("EVTSEKC[half cystine]LEMSC[half cystine]EFD") == [
        proforma.Proteoform(
            sequence="EVTSEKCLEMSCEFD",
            modifications=[
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="half cystine"
                        )
                    ],
                ),
                proforma.Modification(
                    position=11,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="half cystine"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "EVTSEKC[half cystine]LEMSC[half cystine]EFD", True
    ) == [
        proforma.Proteoform(
            sequence="EVTSEKCLEMSCEFD",
            modifications=[
                proforma.Modification(
                    mass=-1.007825,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00798",
                            name="half cystine",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=-1.007825,
                    position=11,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00798",
                            name="half cystine",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "EVTSEKC[MOD:00798]LEMSC[MOD:00798]EFDEVTSEK"
        "C[MOD:00798]LEMSC[MOD:00798]EFD"
    ) == [
        proforma.Proteoform(
            sequence="EVTSEKCLEMSCEFDEVTSEKCLEMSCEFD",
            modifications=[
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD", accession="MOD:00798"
                        )
                    ],
                ),
                proforma.Modification(
                    position=11,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD", accession="MOD:00798"
                        )
                    ],
                ),
                proforma.Modification(
                    position=21,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD", accession="MOD:00798"
                        )
                    ],
                ),
                proforma.Modification(
                    position=26,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD", accession="MOD:00798"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "EVTSEKC[MOD:00798]LEMSC[MOD:00798]EFDEVTSEK"
        "C[MOD:00798]LEMSC[MOD:00798]EFD",
        True,
    ) == [
        proforma.Proteoform(
            sequence="EVTSEKCLEMSCEFDEVTSEKCLEMSCEFD",
            modifications=[
                proforma.Modification(
                    mass=-1.007825,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00798",
                            name="half cystine",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=-1.007825,
                    position=11,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00798",
                            name="half cystine",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=-1.007825,
                    position=21,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00798",
                            name="half cystine",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=-1.007825,
                    position=26,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00798",
                            name="half cystine",
                        )
                    ],
                ),
            ],
        )
    ]


def test_proforma_delta_mass():
    # No prefix.
    assert proforma.parse("EM[+15.9949]EVEES[+79.9663]PEK") == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.Mass(mass=15.9949, controlled_vocabulary=None)
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.Mass(mass=79.9663, controlled_vocabulary=None)
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("EM[+15.9949]EVEES[+79.9663]PEK", True) == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    mass=15.9949,
                    position=1,
                    source=[
                        proforma.Mass(mass=15.9949, controlled_vocabulary=None)
                    ],
                ),
                proforma.Modification(
                    mass=79.9663,
                    position=6,
                    source=[
                        proforma.Mass(mass=79.9663, controlled_vocabulary=None)
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("EM[+15.995]EVEES[-18.01]PEK") == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.Mass(mass=15.995, controlled_vocabulary=None)
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.Mass(mass=-18.01, controlled_vocabulary=None)
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("EM[+15.995]EVEES[-18.01]PEK", True) == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    mass=15.995,
                    position=1,
                    source=[
                        proforma.Mass(mass=15.995, controlled_vocabulary=None)
                    ],
                ),
                proforma.Modification(
                    mass=-18.01,
                    position=6,
                    source=[
                        proforma.Mass(mass=-18.01, controlled_vocabulary=None)
                    ],
                ),
            ],
        )
    ]
    # CV prefix.
    assert proforma.parse("EM[U:+15.9949]EVEES[U:+79.9663]PEK") == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.Mass(
                            mass=15.9949, controlled_vocabulary="UNIMOD"
                        )
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.Mass(
                            mass=79.9663, controlled_vocabulary="UNIMOD"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("EM[U:+15.9949]EVEES[U:+79.9663]PEK", True) == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    mass=15.9949,
                    position=1,
                    source=[
                        proforma.Mass(
                            mass=15.9949, controlled_vocabulary="UNIMOD"
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.9663,
                    position=6,
                    source=[
                        proforma.Mass(
                            mass=79.9663, controlled_vocabulary="UNIMOD"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("EM[U:+15.995]EVEES[U:+79.966]PEK") == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.Mass(
                            mass=15.995, controlled_vocabulary="UNIMOD"
                        )
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.Mass(
                            mass=79.966, controlled_vocabulary="UNIMOD"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("EM[U:+15.995]EVEES[U:+79.966]PEK", True) == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    mass=15.995,
                    position=1,
                    source=[
                        proforma.Mass(
                            mass=15.995, controlled_vocabulary="UNIMOD"
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.966,
                    position=6,
                    source=[
                        proforma.Mass(
                            mass=79.966, controlled_vocabulary="UNIMOD"
                        )
                    ],
                ),
            ],
        )
    ]
    # Experimentally observed prefix.
    assert proforma.parse("EM[U:+15.995]EVEES[Obs:+79.978]PEK") == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.Mass(
                            mass=15.995, controlled_vocabulary="UNIMOD"
                        )
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.Mass(mass=79.978, controlled_vocabulary=None)
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("EM[U:+15.995]EVEES[Obs:+79.978]PEK", True) == [
        proforma.Proteoform(
            sequence="EMEVEESPEK",
            modifications=[
                proforma.Modification(
                    mass=15.995,
                    position=1,
                    source=[
                        proforma.Mass(
                            mass=15.995, controlled_vocabulary="UNIMOD"
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.978,
                    position=6,
                    source=[
                        proforma.Mass(mass=79.978, controlled_vocabulary=None)
                    ],
                ),
            ],
        )
    ]


def test_proforma_gap():
    assert proforma.parse("RTAAX[+367.0537]WT") == [
        proforma.Proteoform(
            sequence="RTAAXWT",
            modifications=[
                proforma.Modification(
                    position=4,
                    source=[
                        proforma.Mass(
                            mass=367.0537, controlled_vocabulary=None
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("RTAAX[+367.0537]WT", True) == [
        proforma.Proteoform(
            sequence="RTAAXWT",
            modifications=[
                proforma.Modification(
                    mass=367.0537,
                    position=4,
                    source=[
                        proforma.Mass(
                            mass=367.0537, controlled_vocabulary=None
                        )
                    ],
                ),
            ],
        )
    ]


def test_proforma_formula():
    assert proforma.parse("SEQUEN[Formula:C12H20O2]CE") == [
        proforma.Proteoform(
            sequence="SEQUENCE",
            modifications=[
                proforma.Modification(
                    position=5, source=[proforma.Formula(formula="C12H20O2")]
                ),
            ],
        )
    ]
    assert proforma.parse("SEQUEN[Formula:C12H20O2]CE", True) == [
        proforma.Proteoform(
            sequence="SEQUENCE",
            modifications=[
                proforma.Modification(
                    mass=196.14632988052,
                    position=5,
                    source=[proforma.Formula(formula="C12H20O2")],
                ),
            ],
        )
    ]
    assert proforma.parse("SEQUEN[Formula:C12 H20 O2]CE") == [
        proforma.Proteoform(
            sequence="SEQUENCE",
            modifications=[
                proforma.Modification(
                    position=5, source=[proforma.Formula(formula="C12 H20 O2")]
                ),
            ],
        )
    ]
    assert proforma.parse("SEQUEN[Formula:C12 H20 O2]CE", True) == [
        proforma.Proteoform(
            sequence="SEQUENCE",
            modifications=[
                proforma.Modification(
                    mass=196.14632988052,
                    position=5,
                    source=[proforma.Formula(formula="C12 H20 O2")],
                ),
            ],
        )
    ]
    # Mass calculation of isotopes is currently still unsupported.
    assert proforma.parse("SEQUEN[Formula:[13C2]CH6N]CE") == [
        proforma.Proteoform(
            sequence="SEQUENCE",
            modifications=[
                proforma.Modification(
                    position=5,
                    source=[
                        proforma.Formula(formula="CH6N", isotopes=["13C2"])
                    ],
                ),
            ],
        )
    ]
    # FIXME: Add isotope support to Pyteomics.
    with pytest.raises(NotImplementedError):
        proforma.parse("SEQUEN[Formula:[13C2]CH6N]CE", True)


def test_proforma_glycan():
    assert proforma.parse("SEQUEN[Glycan:HexNAc1Hex2]CE") == [
        proforma.Proteoform(
            sequence="SEQUENCE",
            modifications=[
                proforma.Modification(
                    position=5,
                    source=[
                        proforma.Glycan(
                            composition=[
                                proforma.Monosaccharide("HexNAc", 1),
                                proforma.Monosaccharide("Hex", 2),
                            ]
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("SEQUEN[Glycan:HexNAc1Hex2]CE", True) == [
        proforma.Proteoform(
            sequence="SEQUENCE",
            modifications=[
                proforma.Modification(
                    mass=527.185019356,
                    position=5,
                    source=[
                        proforma.Glycan(
                            composition=[
                                proforma.Monosaccharide("HexNAc", 1),
                                proforma.Monosaccharide("Hex", 2),
                            ]
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("SEQUEN[Glycan:HexNAc1 Hex2]CE") == [
        proforma.Proteoform(
            sequence="SEQUENCE",
            modifications=[
                proforma.Modification(
                    position=5,
                    source=[
                        proforma.Glycan(
                            composition=[
                                proforma.Monosaccharide("HexNAc", 1),
                                proforma.Monosaccharide("Hex", 2),
                            ]
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("SEQUEN[Glycan:HexNAc1 Hex2]CE", True) == [
        proforma.Proteoform(
            sequence="SEQUENCE",
            modifications=[
                proforma.Modification(
                    mass=527.185019356,
                    position=5,
                    source=[
                        proforma.Glycan(
                            composition=[
                                proforma.Monosaccharide("HexNAc", 1),
                                proforma.Monosaccharide("Hex", 2),
                            ]
                        )
                    ],
                ),
            ],
        )
    ]


def test_proforma_special():
    # N-terminal and C-terminal modifications.
    assert proforma.parse("[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK") == [
        proforma.Proteoform(
            sequence="EMEVNESPEK",
            modifications=[
                proforma.Modification(
                    position="N-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="iTRAQ4plex"
                        )
                    ],
                ),
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Oxidation"
                        )
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "[iTRAQ4plex]-EM[Oxidation]EVNE" "S[Phospho]PEK", True
    ) == [
        proforma.Proteoform(
            sequence="EMEVNESPEK",
            modifications=[
                proforma.Modification(
                    mass=144.102063,
                    position="N-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:214",
                            name="iTRAQ4plex",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.966331,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "[iTRAQ4plex]-EM[U:Oxidation]EVNES[Phospho]PE" "K[iTRAQ4plex]-[Methyl]"
    ) == [
        proforma.Proteoform(
            sequence="EMEVNESPEK",
            modifications=[
                proforma.Modification(
                    position="N-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="iTRAQ4plex"
                        )
                    ],
                ),
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD", name="Oxidation"
                        )
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        )
                    ],
                ),
                proforma.Modification(
                    position=9,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="iTRAQ4plex"
                        )
                    ],
                ),
                proforma.Modification(
                    position="C-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Methyl"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "[iTRAQ4plex]-EM[U:Oxidation]EVNES[Phospho]PE"
        "K[iTRAQ4plex]-[Methyl]",
        True,
    ) == [
        proforma.Proteoform(
            sequence="EMEVNESPEK",
            modifications=[
                proforma.Modification(
                    mass=144.102063,
                    position="N-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:214",
                            name="iTRAQ4plex",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.966331,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=144.102063,
                    position=9,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:214",
                            name="iTRAQ4plex",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=14.01565,
                    position="C-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:34",
                            name="Methyl",
                        )
                    ],
                ),
            ],
        )
    ]
    # Labile modifications.
    assert proforma.parse(
        "{Glycan:Hex}EM[U:Oxidation]EVNES[Phospho]PE" "K[iTRAQ4plex]"
    ) == [
        proforma.Proteoform(
            sequence="EMEVNESPEK",
            modifications=[
                proforma.Modification(
                    position="labile",
                    source=[
                        proforma.Glycan(
                            composition=[proforma.Monosaccharide("Hex", 1)]
                        )
                    ],
                ),
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD", name="Oxidation"
                        )
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        )
                    ],
                ),
                proforma.Modification(
                    position=9,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="iTRAQ4plex"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "{Glycan:Hex}EM[U:Oxidation]EVNES[Phospho]PE" "K[iTRAQ4plex]", True
    ) == [
        proforma.Proteoform(
            sequence="EMEVNESPEK",
            modifications=[
                proforma.Modification(
                    mass=162.052823418,
                    position="labile",
                    source=[
                        proforma.Glycan(
                            composition=[proforma.Monosaccharide("Hex", 1)]
                        )
                    ],
                ),
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.966331,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=144.102063,
                    position=9,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:214",
                            name="iTRAQ4plex",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "{Glycan:Hex}[iTRAQ4plex]-EM[Oxidation]EVNE"
        "S[Phospho]PEK[iTRAQ4plex]"
    ) == [
        proforma.Proteoform(
            sequence="EMEVNESPEK",
            modifications=[
                proforma.Modification(
                    position="labile",
                    source=[
                        proforma.Glycan(
                            composition=[proforma.Monosaccharide("Hex", 1)]
                        )
                    ],
                ),
                proforma.Modification(
                    position="N-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="iTRAQ4plex"
                        )
                    ],
                ),
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Oxidation"
                        )
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        )
                    ],
                ),
                proforma.Modification(
                    position=9,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="iTRAQ4plex"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "{Glycan:Hex}[iTRAQ4plex]-EM[Oxidation]EVNE"
        "S[Phospho]PEK[iTRAQ4plex]",
        True,
    ) == [
        proforma.Proteoform(
            sequence="EMEVNESPEK",
            modifications=[
                proforma.Modification(
                    mass=162.052823418,
                    position="labile",
                    source=[
                        proforma.Glycan(
                            composition=[proforma.Monosaccharide("Hex", 1)]
                        )
                    ],
                ),
                proforma.Modification(
                    mass=144.102063,
                    position="N-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:214",
                            name="iTRAQ4plex",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.966331,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=144.102063,
                    position=9,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:214",
                            name="iTRAQ4plex",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "{Glycan:Hex}[iTRAQ4plex]-EM[Oxidation]EVNE"
        "S[Phospho]PEK[iTRAQ4plex]-[Methyl]"
    ) == [
        proforma.Proteoform(
            sequence="EMEVNESPEK",
            modifications=[
                proforma.Modification(
                    position="labile",
                    source=[
                        proforma.Glycan(
                            composition=[proforma.Monosaccharide("Hex", 1)]
                        )
                    ],
                ),
                proforma.Modification(
                    position="N-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="iTRAQ4plex"
                        )
                    ],
                ),
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Oxidation"
                        )
                    ],
                ),
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        )
                    ],
                ),
                proforma.Modification(
                    position=9,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="iTRAQ4plex"
                        )
                    ],
                ),
                proforma.Modification(
                    position="C-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Methyl"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "{Glycan:Hex}[iTRAQ4plex]-EM[Oxidation]EVNE"
        "S[Phospho]PEK[iTRAQ4plex]-[Methyl]",
        True,
    ) == [
        proforma.Proteoform(
            sequence="EMEVNESPEK",
            modifications=[
                proforma.Modification(
                    mass=162.052823418,
                    position="labile",
                    source=[
                        proforma.Glycan(
                            composition=[proforma.Monosaccharide("Hex", 1)]
                        )
                    ],
                ),
                proforma.Modification(
                    mass=144.102063,
                    position="N-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:214",
                            name="iTRAQ4plex",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.966331,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=144.102063,
                    position=9,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:214",
                            name="iTRAQ4plex",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=14.01565,
                    position="C-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:34",
                            name="Methyl",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("{Glycan:Hex}{Glycan:NeuAc}EMEVNESPEK") == [
        proforma.Proteoform(
            sequence="EMEVNESPEK",
            modifications=[
                proforma.Modification(
                    position="labile",
                    source=[
                        proforma.Glycan(
                            composition=[proforma.Monosaccharide("Hex", 1)]
                        )
                    ],
                ),
                proforma.Modification(
                    position="labile",
                    source=[
                        proforma.Glycan(
                            composition=[proforma.Monosaccharide("NeuAc", 1)]
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("{Glycan:Hex}{Glycan:NeuAc}EMEVNESPEK", True) == [
        proforma.Proteoform(
            sequence="EMEVNESPEK",
            modifications=[
                proforma.Modification(
                    mass=162.052823418,
                    position="labile",
                    source=[
                        proforma.Glycan(
                            composition=[proforma.Monosaccharide("Hex", 1)]
                        )
                    ],
                ),
                proforma.Modification(
                    mass=291.095416506,
                    position="labile",
                    source=[
                        proforma.Glycan(
                            composition=[proforma.Monosaccharide("NeuAc", 1)]
                        )
                    ],
                ),
            ],
        )
    ]


def test_proforma_ambiguous_position():
    # Unknown modification position.
    assert proforma.parse("[Phospho]?EM[Oxidation]EVTSESPEK") == [
        proforma.Proteoform(
            sequence="EMEVTSESPEK",
            modifications=[
                proforma.Modification(
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        )
                    ],
                ),
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Oxidation"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("[Phospho]?EM[Oxidation]EVTSESPEK", True) == [
        proforma.Proteoform(
            sequence="EMEVTSESPEK",
            modifications=[
                proforma.Modification(
                    mass=79.966331,
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "[Phospho][Phospho]?[Acetyl]-EM[Oxidation]" "EVTSESPEK"
    ) == [
        proforma.Proteoform(
            sequence="EMEVTSESPEK",
            modifications=[
                proforma.Modification(
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        )
                    ],
                ),
                proforma.Modification(
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        )
                    ],
                ),
                proforma.Modification(
                    position="N-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Acetyl"
                        )
                    ],
                ),
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Oxidation"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "[Phospho][Phospho]?[Acetyl]-EM[Oxidation]" "EVTSESPEK", True
    ) == [
        proforma.Proteoform(
            sequence="EMEVTSESPEK",
            modifications=[
                proforma.Modification(
                    mass=79.966331,
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.966331,
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=42.010565,
                    position="N-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:1",
                            name="Acetyl",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "[Phospho]^2[Methyl]?[Acetyl]-EM[Oxidation]" "EVTSESPEK"
    ) == [
        proforma.Proteoform(
            sequence="EMEVTSESPEK",
            modifications=[
                proforma.Modification(
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        )
                    ],
                ),
                proforma.Modification(
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        )
                    ],
                ),
                proforma.Modification(
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Methyl"
                        )
                    ],
                ),
                proforma.Modification(
                    position="N-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Acetyl"
                        )
                    ],
                ),
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Oxidation"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "[Phospho]^2[Methyl]?[Acetyl]-EM[Oxidation]" "EVTSESPEK", True
    ) == [
        proforma.Proteoform(
            sequence="EMEVTSESPEK",
            modifications=[
                proforma.Modification(
                    mass=79.966331,
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.966331,
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=14.01565,
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:34",
                            name="Methyl",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=42.010565,
                    position="N-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:1",
                            name="Acetyl",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("[Phospho]^2?[Acetyl]-EM[Oxidation]EVTSESPEK") == [
        proforma.Proteoform(
            sequence="EMEVTSESPEK",
            modifications=[
                proforma.Modification(
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        )
                    ],
                ),
                proforma.Modification(
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        )
                    ],
                ),
                proforma.Modification(
                    position="N-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Acetyl"
                        )
                    ],
                ),
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Oxidation"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "[Phospho]^2?[Acetyl]-EM[Oxidation]" "EVTSESPEK", True
    ) == [
        proforma.Proteoform(
            sequence="EMEVTSESPEK",
            modifications=[
                proforma.Modification(
                    mass=79.966331,
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=79.966331,
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=42.010565,
                    position="N-term",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:1",
                            name="Acetyl",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
            ],
        )
    ]
    # Multiple possible modification positions.
    assert proforma.parse("EM[Oxidation]EVT[#g1]S[#g1]ES[Phospho#g1]PEK") == [
        proforma.Proteoform(
            sequence="EMEVTSESPEK",
            modifications=[
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Oxidation"
                        )
                    ],
                ),
                proforma.Modification(
                    position=4,
                    label=proforma.Label(proforma.LabelType.GENERAL, "g1"),
                ),
                proforma.Modification(
                    position=5,
                    label=proforma.Label(proforma.LabelType.GENERAL, "g1"),
                ),
                proforma.Modification(
                    position=7,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        )
                    ],
                    label=proforma.Label(proforma.LabelType.GENERAL, "g1"),
                ),
            ],
        )
    ]
    assert proforma.parse(
        "EM[Oxidation]EVT[#g1]S[#g1]ES[Phospho#g1]PEK", True
    ) == [
        proforma.Proteoform(
            sequence="EMEVTSESPEK",
            modifications=[
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
                proforma.Modification(
                    position=4,
                    label=proforma.Label(proforma.LabelType.GENERAL, "g1"),
                ),
                proforma.Modification(
                    position=5,
                    label=proforma.Label(proforma.LabelType.GENERAL, "g1"),
                ),
                proforma.Modification(
                    mass=79.966331,
                    position=7,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                    label=proforma.Label(proforma.LabelType.GENERAL, "g1"),
                ),
            ],
        )
    ]
    with pytest.raises(ValueError):
        proforma.parse("EM[Oxidation]EVT[#g1]S[Phospho#g1]ES[Phospho#g1]PEK")
    # Ranges of modification positions.
    assert proforma.parse("PROT(EOSFORMS)[+19.0523]ISK") == [
        proforma.Proteoform(
            sequence="PROTEOSFORMSISK",
            modifications=[
                proforma.Modification(
                    position=(4, 11), source=[proforma.Mass(mass=19.0523)]
                ),
            ],
        )
    ]
    assert proforma.parse("PROT(EOSFORMS)[+19.0523]ISK", True) == [
        proforma.Proteoform(
            sequence="PROTEOSFORMSISK",
            modifications=[
                proforma.Modification(
                    mass=19.0523,
                    position=(4, 11),
                    source=[proforma.Mass(mass=19.0523)],
                ),
            ],
        )
    ]
    assert proforma.parse("PROT(EOC[Carbamidomethyl]FORMS)[+19.0523]ISK") == [
        proforma.Proteoform(
            sequence="PROTEOCFORMSISK",
            modifications=[
                proforma.Modification(
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Carbamidomethyl"
                        )
                    ],
                ),
                proforma.Modification(
                    position=(4, 11), source=[proforma.Mass(mass=19.0523)]
                ),
            ],
        )
    ]
    assert proforma.parse(
        "PROT(EOC[Carbamidomethyl]FORMS)[+19.0523]ISK", True
    ) == [
        proforma.Proteoform(
            sequence="PROTEOCFORMSISK",
            modifications=[
                proforma.Modification(
                    mass=57.021464,
                    position=6,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:4",
                            name="Carbamidomethyl",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=19.0523,
                    position=(4, 11),
                    source=[proforma.Mass(mass=19.0523)],
                ),
            ],
        )
    ]
    assert proforma.parse("P(ROT(EOSFORMS)[+19.0523]IS)[+19.0523]K") == [
        proforma.Proteoform(
            sequence="PROTEOSFORMSISK",
            modifications=[
                proforma.Modification(
                    position=(4, 11), source=[proforma.Mass(mass=19.0523)]
                ),
                proforma.Modification(
                    position=(1, 13), source=[proforma.Mass(mass=19.0523)]
                ),
            ],
        )
    ]
    assert proforma.parse("P(ROT(EOSFORMS)[+19.0523]IS)[+19.0523]K", True) == [
        proforma.Proteoform(
            sequence="PROTEOSFORMSISK",
            modifications=[
                proforma.Modification(
                    mass=19.0523,
                    position=(4, 11),
                    source=[proforma.Mass(mass=19.0523)],
                ),
                proforma.Modification(
                    mass=19.0523,
                    position=(1, 13),
                    source=[proforma.Mass(mass=19.0523)],
                ),
            ],
        )
    ]
    # Modification localization score.
    assert proforma.parse(
        "EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]E" "S[Phospho#g1(0.90)]PEK"
    ) == [
        proforma.Proteoform(
            sequence="EMEVTSESPEK",
            modifications=[
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Oxidation"
                        )
                    ],
                ),
                proforma.Modification(
                    position=4,
                    label=proforma.Label(
                        proforma.LabelType.GENERAL, "g1", 0.01
                    ),
                ),
                proforma.Modification(
                    position=5,
                    label=proforma.Label(
                        proforma.LabelType.GENERAL, "g1", 0.09
                    ),
                ),
                proforma.Modification(
                    position=7,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        )
                    ],
                    label=proforma.Label(
                        proforma.LabelType.GENERAL, "g1", 0.90
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse(
        "EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]E" "S[Phospho#g1(0.90)]PEK",
        True,
    ) == [
        proforma.Proteoform(
            sequence="EMEVTSESPEK",
            modifications=[
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
                proforma.Modification(
                    position=4,
                    label=proforma.Label(
                        proforma.LabelType.GENERAL, "g1", 0.01
                    ),
                ),
                proforma.Modification(
                    position=5,
                    label=proforma.Label(
                        proforma.LabelType.GENERAL, "g1", 0.09
                    ),
                ),
                proforma.Modification(
                    mass=79.966331,
                    position=7,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                    label=proforma.Label(
                        proforma.LabelType.GENERAL, "g1", 0.90
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse(
        "[Phospho#s1]?EM[Oxidation]EVT[#s1(0.01)]"
        "S[#s1(0.09)]ES[#s1(0.90)]PEK"
    ) == [
        proforma.Proteoform(
            sequence="EMEVTSESPEK",
            modifications=[
                proforma.Modification(
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        )
                    ],
                    label=proforma.Label(proforma.LabelType.GENERAL, "s1"),
                ),
                proforma.Modification(
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Oxidation"
                        )
                    ],
                ),
                proforma.Modification(
                    position=4,
                    label=proforma.Label(
                        proforma.LabelType.GENERAL, "s1", 0.01
                    ),
                ),
                proforma.Modification(
                    position=5,
                    label=proforma.Label(
                        proforma.LabelType.GENERAL, "s1", 0.09
                    ),
                ),
                proforma.Modification(
                    position=7,
                    label=proforma.Label(
                        proforma.LabelType.GENERAL, "s1", 0.90
                    ),
                ),
            ],
        )
    ]
    assert proforma.parse(
        "[Phospho#s1]?EM[Oxidation]EVT[#s1(0.01)]"
        "S[#s1(0.09)]ES[#s1(0.90)]PEK",
        True,
    ) == [
        proforma.Proteoform(
            sequence="EMEVTSESPEK",
            modifications=[
                proforma.Modification(
                    mass=79.966331,
                    position="unknown",
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        )
                    ],
                    label=proforma.Label(proforma.LabelType.GENERAL, "s1"),
                ),
                proforma.Modification(
                    mass=15.994915,
                    position=1,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
                proforma.Modification(
                    position=4,
                    label=proforma.Label(
                        proforma.LabelType.GENERAL, "s1", 0.01
                    ),
                ),
                proforma.Modification(
                    position=5,
                    label=proforma.Label(
                        proforma.LabelType.GENERAL, "s1", 0.09
                    ),
                ),
                proforma.Modification(
                    position=7,
                    label=proforma.Label(
                        proforma.LabelType.GENERAL, "s1", 0.90
                    ),
                ),
            ],
        )
    ]


def test_proforma_global_modification():
    assert proforma.parse(
        "<[S-carboxamidomethyl-L-cysteine]@C>" "ATPEILTCNSIGCLK"
    ) == [
        proforma.Proteoform(
            sequence="ATPEILTCNSIGCLK",
            modifications=[
                proforma.Modification(
                    position=7,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None,
                            name="S-carboxamidomethyl-L-cysteine",
                        )
                    ],
                ),
                proforma.Modification(
                    position=12,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None,
                            name="S-carboxamidomethyl-L-cysteine",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "<[S-carboxamidomethyl-L-cysteine]@C>" "ATPEILTCNSIGCLK", True
    ) == [
        proforma.Proteoform(
            sequence="ATPEILTCNSIGCLK",
            modifications=[
                proforma.Modification(
                    mass=57.021464,
                    position=7,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:01060",
                            name="S-carboxamidomethyl-L-cysteine",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=57.021464,
                    position=12,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:01060",
                            name="S-carboxamidomethyl-L-cysteine",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("<[MOD:01090]@C>ATPEILTCNSIGCLK") == [
        proforma.Proteoform(
            sequence="ATPEILTCNSIGCLK",
            modifications=[
                proforma.Modification(
                    position=7,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD", accession="MOD:01090"
                        )
                    ],
                ),
                proforma.Modification(
                    position=12,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD", accession="MOD:01090"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("<[MOD:01090]@C>ATPEILTCNSIGCLK", True) == [
        proforma.Proteoform(
            sequence="ATPEILTCNSIGCLK",
            modifications=[
                proforma.Modification(
                    mass=57.021464,
                    position=7,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:01090",
                            name="iodoacetamide derivatized amino-terminal "
                            "residue",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=57.021464,
                    position=12,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:01090",
                            name="iodoacetamide derivatized amino-terminal "
                            "residue",
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("<[Oxidation]@C,M>MTPEILTCNSIGCLK") == [
        proforma.Proteoform(
            sequence="MTPEILTCNSIGCLK",
            modifications=[
                proforma.Modification(
                    position=0,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Oxidation"
                        )
                    ],
                ),
                proforma.Modification(
                    position=7,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Oxidation"
                        )
                    ],
                ),
                proforma.Modification(
                    position=12,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Oxidation"
                        )
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("<[Oxidation]@C,M>MTPEILTCNSIGCLK", True) == [
        proforma.Proteoform(
            sequence="MTPEILTCNSIGCLK",
            modifications=[
                proforma.Modification(
                    mass=15.994915,
                    position=0,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=15.994915,
                    position=7,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
                proforma.Modification(
                    mass=15.994915,
                    position=12,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:35",
                            name="Oxidation",
                        )
                    ],
                ),
            ],
        )
    ]


def test_proforma_info():
    assert proforma.parse("ELV[INFO:AnyString]IS") == [
        proforma.Proteoform(
            sequence="ELVIS",
            modifications=[
                proforma.Modification(
                    position=2, source=[proforma.Info("AnyString")]
                ),
            ],
        )
    ]
    assert proforma.parse("ELV[INFO:AnyString]IS", True) == [
        proforma.Proteoform(
            sequence="ELVIS",
            modifications=[
                proforma.Modification(
                    position=2, source=[proforma.Info("AnyString")]
                ),
            ],
        )
    ]
    assert proforma.parse("ELV[info:AnyString]IS") == [
        proforma.Proteoform(
            sequence="ELVIS",
            modifications=[
                proforma.Modification(
                    position=2, source=[proforma.Info("AnyString")]
                ),
            ],
        )
    ]
    assert proforma.parse("ELV[info:AnyString]IS", True) == [
        proforma.Proteoform(
            sequence="ELVIS",
            modifications=[
                proforma.Modification(
                    position=2, source=[proforma.Info("AnyString")]
                ),
            ],
        )
    ]
    assert proforma.parse("ELVIS[Phospho|INFO:newly discovered]K") == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        ),
                        proforma.Info("newly discovered"),
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("ELVIS[Phospho|INFO:newly discovered]K", True) == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    mass=79.966331,
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        ),
                        proforma.Info("newly discovered"),
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "ELVIS[Phospho|INFO:newly discovered|" "INFO:really awesome]K"
    ) == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        ),
                        proforma.Info("newly discovered"),
                        proforma.Info("really awesome"),
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "ELVIS[Phospho|INFO:newly discovered|" "INFO:really awesome]K", True
    ) == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    mass=79.966331,
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        ),
                        proforma.Info("newly discovered"),
                        proforma.Info("really awesome"),
                    ],
                ),
            ],
        )
    ]


def test_proforma_pipe():
    assert proforma.parse("ELVIS[U:Phospho|+79.966331]K") == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD", name="Phospho"
                        ),
                        proforma.Mass(mass=79.966331),
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("ELVIS[U:Phospho|+79.966331]K", True) == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    mass=79.966331,
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        ),
                        proforma.Mass(mass=79.966331),
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("ELVIS[U:Phospho|Obs:+79.978]K") == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD", name="Phospho"
                        ),
                        proforma.Mass(mass=79.978),
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("ELVIS[U:Phospho|Obs:+79.978]K", True) == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    mass=79.966331,
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        ),
                        proforma.Mass(mass=79.978),
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("ELVIS[Phospho|O-phospho-L-serine]K") == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        ),
                        proforma.CvEntry(
                            controlled_vocabulary=None,
                            name="O-phospho-L-serine",
                        ),
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("ELVIS[Phospho|O-phospho-L-serine]K", True) == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    mass=79.966331,
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        ),
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00046",
                            name="O-phospho-L-serine",
                        ),
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("ELVIS[UNIMOD:21|MOD:00046]K") == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                        ),
                        proforma.CvEntry(
                            controlled_vocabulary="MOD", accession="MOD:00046"
                        ),
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("ELVIS[UNIMOD:21|MOD:00046]K", True) == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    mass=79.966331,
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        ),
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00046",
                            name="O-phospho-L-serine",
                        ),
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("ELVIS[UNIMOD:21|Phospho]K") == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                        ),
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        ),
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("ELVIS[UNIMOD:21|Phospho]K", True) == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    mass=79.966331,
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        ),
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        ),
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "ELVIS[Phospho|O-phospho-L-serine|Obs:+79.966]K"
    ) == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        ),
                        proforma.CvEntry(
                            controlled_vocabulary=None,
                            name="O-phospho-L-serine",
                        ),
                        proforma.Mass(79.966),
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse(
        "ELVIS[Phospho|O-phospho-L-serine|Obs:+79.966]K", True
    ) == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    mass=79.966331,
                    position=4,
                    source=[
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        ),
                        proforma.CvEntry(
                            controlled_vocabulary="MOD",
                            accession="MOD:00046",
                            name="O-phospho-L-serine",
                        ),
                        proforma.Mass(79.966),
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("ELVIS[Obs:+79.966|Phospho|Sulfo]K") == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    position=4,
                    source=[
                        proforma.Mass(79.966),
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Phospho"
                        ),
                        proforma.CvEntry(
                            controlled_vocabulary=None, name="Sulfo"
                        ),
                    ],
                ),
            ],
        )
    ]
    assert proforma.parse("ELVIS[Obs:+79.966|Phospho|Sulfo]K", True) == [
        proforma.Proteoform(
            sequence="ELVISK",
            modifications=[
                proforma.Modification(
                    mass=79.966,
                    position=4,
                    source=[
                        proforma.Mass(79.966),
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:21",
                            name="Phospho",
                        ),
                        proforma.CvEntry(
                            controlled_vocabulary="UNIMOD",
                            accession="UNIMOD:40",
                            name="Sulfo",
                        ),
                    ],
                ),
            ],
        )
    ]
