import collections
import copy
import enum
import functools
import json
import os
import pickle
import urllib.request
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple, Union
from urllib.error import URLError

import fastobo
import lark
try:
    import pyteomics.cmass as pmass
except ImportError:
    import pyteomics.mass as pmass


# Set to None to disable caching.
cache_dir = os.path.join(os.path.expanduser('~'), '.cache', 'spectrum_utils')


LookupType = enum.Enum('LookupType', 'ACCESSION NAME MASS')
LabelType = enum.Enum('LabelType', 'XL BRANCH GENERAL')


@dataclass
class Ion:
    ion: str


@dataclass
class Charge:
    charge: int
    ions: Optional[List[Ion]] = None


@dataclass
class CvEntry:
    controlled_vocabulary: Optional[str]
    accession: Optional[str] = None
    name: Optional[str] = None


@dataclass
class Mass:
    mass: float
    controlled_vocabulary: Optional[str] = None


@dataclass
class Formula:
    formula: Optional[str] = None
    isotopes: Optional[List[str]] = None


@dataclass
class Monosaccharide:
    monosaccharide: str
    count: int


@dataclass
class Glycan:
    composition: List[Monosaccharide]


@dataclass
class Info:
    message: str


@dataclass
class Label:
    type: LabelType
    label: str
    score: Optional[float] = None


@dataclass
class Modification:
    mass: Optional[float] = None
    position: Optional[Union[int, Tuple[int, int], str]] = None
    source: List[Union[CvEntry, Mass, Formula, Glycan, Info]] = None
    label: Optional[Label] = None


@dataclass
class Proteoform:
    sequence: str
    modifications: List[Modification]
    charge: Optional[Charge] = None


# noinspection PyMethodMayBeStatic, PyPep8Naming
class ProFormaTransformer(lark.Transformer):
    sequence: List[str]
    modifications: List[Modification]
    global_modifications: Dict[str, List[Modification]]
    range_pos: List[int]

    def __init__(self):
        super().__init__()
        self.sequence, self.modifications = [], []
        self.global_modifications = collections.defaultdict(list)
        self.range_pos = []

    def proforma(self, tree) -> List[Proteoform]:
        return [proteoform for proteoform in tree
                if isinstance(proteoform, Proteoform)]

    def proteoform(self, tree) -> Proteoform:
        sequence = ''.join(self.sequence)
        # Apply global modifications to the relevant residues.
        for i, aa in enumerate(sequence):
            if aa in self.global_modifications:
                for mod in self.global_modifications[aa]:
                    mod = copy.copy(mod)
                    mod.position = i
                    self.modifications.append(mod)
        charge = tree[-1] if len(tree) > 1 else None
        proteoform = Proteoform(sequence=sequence,
                                modifications=self.modifications,
                                charge=charge)
        # Reset class variables.
        self.sequence, self.modifications = [], []
        self.global_modifications = collections.defaultdict(list)
        return proteoform

    def peptide(self, _) -> None:
        # Residues and modifications have been written to the class variables.
        pass

    def aa(self, tree) -> None:
        position = len(self.sequence)
        self.sequence.append(tree[0])
        # An amino acid token can be followed by (i) a modification on that
        # residue, or (ii) a label (linking it to another modified residue).
        if len(tree) == 2:
            if isinstance(tree[1], Label):
                self.modifications.append(
                    Modification(position=position, label=tree[1]))
            else:
                tree[1].position = position
                self.modifications.append(tree[1])

    def AA(self, token) -> str:
        return token.value

    def mod_global(self, mods) -> None:
        if len(mods) == 1:
            # Global isotope.
            self.modifications.append(Modification(
                position='global', source=[Formula(isotopes=mods)]))
        else:
            # Global modification on a specific residue.
            for aa in mods[1:]:
                self.global_modifications[aa].append(mods[0])

    def ISOTOPE(self, token) -> str:
        return token.value

    def mod_unknown_pos(self, mods) -> None:
        for mod in mods:
            if isinstance(mod, Modification):
                mod.position = 'unknown'
                self.modifications.append(mod)
            else:
                # Modification count.
                for _ in range(int(mod) - 1):
                    self.modifications.append(
                        copy.copy(self.modifications[-1]))

    def mod(self, mod_annotations) -> Modification:
        mod = Modification(source=[])
        for mod_annotation in mod_annotations:
            if isinstance(mod_annotation, Label):
                mod.label = mod_annotation
            else:
                mod.source.append(mod_annotation)
        return mod

    def mod_labile(self, mod_annotations) -> None:
        mod = self.mod(mod_annotations)
        mod.position = 'labile'
        self.modifications.append(mod)

    def MOD_COUNT(self, token) -> int:
        return int(token.value)

    def mod_n_term(self, mods) -> None:
        for mod in mods:
            if isinstance(mod, Label):
                mod = Modification(label=mod)
            mod.position = 'N-term'
            self.modifications.append(mod)

    def mod_c_term(self, mods) -> None:
        for mod in mods:
            if isinstance(mod, Label):
                mod = Modification(label=mod)
            mod.position = 'C-term'
            self.modifications.append(mod)

    def mod_range(self, tree) -> None:
        _, position, *mods = tree
        for mod in mods:
            mod.position = position
            self.modifications.append(mod)

    def mod_range_pos(self, _) -> Tuple[int, int]:
        return self.range_pos.pop(), len(self.sequence) - 1

    def MOD_RANGE_L(self, _) -> None:
        self.range_pos.append(len(self.sequence))

    def mod_name(self, tree) -> CvEntry:
        cv, name = tree if len(tree) == 2 else (None, tree[0])
        return CvEntry(controlled_vocabulary=cv, name=name)

    def CV_ABBREV_OPT(self, token) -> str:
        return self.CV_ABBREV(token)

    def CV_ABBREV(self, token) -> str:
        return {'U': 'UNIMOD', 'M': 'MOD', 'R': 'RESID', 'X': 'XLMOD',
                'G': 'GNO'}[token.value.upper()]

    def mod_accession(self, tree) -> CvEntry:
        return CvEntry(controlled_vocabulary=tree[0],
                       accession=f'{tree[0]}:{tree[1]}')

    def CV_NAME(self, token) -> str:
        return token.value

    def mod_mass(self, tree) -> Mass:
        if len(tree) == 1:
            return Mass(mass=tree[0])
        elif tree[0] == 'Obs':
            return Mass(mass=tree[1])
        else:
            return Mass(mass=tree[1], controlled_vocabulary=tree[0])

    def MOD_MASS_OBS(self, _) -> str:
        return 'Obs'

    def MOD_MASS(self, token) -> float:
        return float(token.value)

    def mod_formula(self, tree) -> Formula:
        *isotopes, formula = tree if len(tree) > 1 else (tree[0],)
        return Formula(
            formula=formula, isotopes=isotopes if isotopes else None)

    def FORMULA(self, token) -> str:
        return token.value

    def mod_glycan(self, tree) -> Glycan:
        # Ignore optional whitespace in monosaccharide composition.
        return Glycan(
            composition=[t for t in tree if isinstance(t, Monosaccharide)])

    def monosaccharide(self, tree) -> Monosaccharide:
        return Monosaccharide(
            tree[0].value, int(tree[1].value) if len(tree) > 1 else 1)

    def info(self, tree) -> Info:
        return Info(tree[0])

    def mod_label(self, tree) -> Label:
        return Label(type=tree[0][0], label=tree[0][1],
                     score=tree[1] if len(tree) > 1 else None)

    def MOD_LABEL_XL(self, token) -> Tuple[LabelType, str]:
        return LabelType.XL, token.value

    def MOD_LABEL_BRANCH(self, token) -> Tuple[LabelType, str]:
        return LabelType.BRANCH, token.value

    def MOD_LABEL(self, token) -> Tuple[LabelType, str]:
        return LabelType.GENERAL, token.value

    def MOD_SCORE(self, token) -> float:
        return float(token)

    def charge(self, tree) -> Charge:
        return Charge(charge=int(tree[0].value),
                      ions=tree[1] if len(tree) > 1 else None)

    def ion(self, tree) -> List[Ion]:
        return [Ion(ion) for ion in tree]

    def TEXT(self, token) -> str:
        return token.value


def parse(proforma: str, resolve_mods: bool = False) -> List[Proteoform]:
    """
    Parse a ProForma-encoded string.

    The string is parsed by building an abstract syntax tree to interpret the
    ProForma language. Next, modifications are localized to residue positions
    and translated into Python objects.

    Parameters
    ----------
    proforma : str
        The ProForma string.
    resolve_mods : bool
        Resolve modifications by retrieving the modification masses from the
        corresponding controlled vocabularies or computing the masses from
        the specified molecular or glycan compositions.

    Returns
    -------
    List[Proteoform]
        A list of proteoforms with their modifications (and additional)
        information specified by the ProForma string.

    Raises
    ------
    URLError
        If a controlled vocabulary could not be retrieved from its online
        resource.
    KeyError
        If a term was not found in its controlled vocabulary.
    NotImplementedError
        Mass calculation using a molecular formula that contains isotopes (not
        supported by Pyteomics mass calculation).
    ValueError
        If no mass was specified for a GNO term or its parent terms.
    """
    with open('spectrum_utils/proforma.ebnf') as f_in:
        parser = lark.Lark(f_in.read(), start='proforma', parser='earley',
                           lexer='dynamic_complete')
        proteoforms = ProFormaTransformer().transform(parser.parse(proforma))
        if resolve_mods:
            for proteoform in proteoforms:
                for mod in proteoform.modifications:
                    if mod.source is None:
                        # Only a label without a modification mass/source.
                        continue
                    for source in reversed(mod.source):
                        if isinstance(source, CvEntry):
                            if source.controlled_vocabulary is None:
                                for cv in ('UNIMOD', 'MOD'):
                                    try:
                                        source.controlled_vocabulary = cv
                                        mod.mass = _resolve_cv(source)
                                        break
                                    except KeyError:
                                        pass
                                else:
                                    raise KeyError(
                                        f'Term "{source.name}" not found in '
                                        f'UNIMOD or PSI-MOD')
                            else:
                                mod.mass = _resolve_cv(source)
                        elif isinstance(source, Mass):
                            mod.mass = _resolve_mass(source)
                        elif isinstance(source, Formula):
                            mod.mass = _resolve_formula(source)
                        elif isinstance(source, Glycan):
                            mod.mass = _resolve_glycan(source)
        return proteoforms


def _resolve_cv(cv_entry: CvEntry) -> float:
    """
    Retrieve the modification mass from a controlled vocabulary entry.

    Parameters
    ----------
    cv_entry : CvEntry
        The CV entry whose modification mass will be retrieved.

    Returns
    -------
    float
        The modification mass.

    Raises
    ------
    KeyError
        If the term was not found in its controlled vocabulary.
    URLError
        If the controlled vocabulary could not be retrieved from its online
        resource.
    ValueError
        - If an unknown controlled vocabulary identifier is specified.
        - If no mass was specified for a GNO term or its parent terms.
    """
    cv_by_accession, cv_by_name = _import_cv(
        cv_entry.controlled_vocabulary, cache_dir)
    key = lookup_type = None
    try:
        if cv_entry.accession is not None:
            key, lookup_type = cv_entry.accession, 'accession'
            cv_entry.name = cv_by_accession[key][1]
            return cv_by_accession[key][0]
        elif cv_entry.name is not None:
            key, lookup_type = cv_entry.name, 'name'
            cv_entry.accession = cv_by_name[key][1]
            return cv_by_name[key][0]
    except KeyError:
        raise KeyError(
            f'Term {key} not found in the {cv_entry.controlled_vocabulary} '
            f'controlled vocabulary by {lookup_type}')


def _resolve_mass(mass: Mass) -> float:
    """
    Return the mass of a mass-based modification.

    Parameters
    ----------
    mass : Mass
        The mass-based modification.

    Returns
    -------
    float
        The modification mass.
    """
    return mass.mass


def _resolve_formula(formula: Formula) -> float:
    """
    Calculate the mass of a molecular formula.

    Parameters
    ----------
    formula : Formula
        The molecular formula whose mass will be computed.

    Returns
    -------
    float
        The modification mass.

    Raises
    ------
    NotImplementedError
        Mass calculation using a molecular formula that contains isotopes (not
        supported by Pyteomics mass calculation).
    """
    if formula.isotopes is not None:
        # FIXME: Add isotope support to Pyteomics.
        raise NotImplementedError('Mass calculation of molecular formulas with'
                                  ' isotopes is currently not supported')
    # Make sure there are no spaces in the molecular formula (Pyteomics mass
    # calculation can't handle those).
    return pmass.calculate_mass(formula=formula.formula.replace(' ', ''))


def _resolve_glycan(glycan: Glycan) -> float:
    """
    Calculate the mass of a glycan based on its monosaccharide composition.

    Parameters
    ----------
    glycan : Glycan
        The glycan composition whose mass will be computed.

    Returns
    -------
    float
        The modification mass.

    Raises
    ------
    URLError
        If the monosaccharide definitions could not be retrieved from its
        online resource.
    """
    mono = _load_from_cache(cache_dir, 'mono.pkl')
    if mono is None:
        with urllib.request.urlopen(
                'https://raw.githubusercontent.com/HUPO-PSI/ProForma/master/'
                'monosaccharides/mono.obo.json') as response:
            if 400 <= response.getcode() < 600:
                raise URLError('Failed to retrieve the monosaccharide '
                               'definitions from its online resource')
            mono = json.loads(response.read().decode(
                response.info().get_param('charset', 'utf-8')))
        mono = {term['name']: float(term['has_monoisotopic_mass'])
                for term in mono['terms'].values()}
        _store_in_cache(cache_dir, 'mono.pkl', mono)
    return sum([mono[m.monosaccharide] * m.count for m in glycan.composition])


@functools.lru_cache
def _import_cv(cv_id: str, cache: Optional[str]) \
        -> Tuple[Dict[str, Tuple[float, str]], Dict[str, Tuple[float, str]]]:
    """
    Import a ProForma controlled vocabulary from its online resource.

    Parameters
    ----------
    cv_id : str
        The controlled vocabulary identifier.
    cache : Optional[str]
        Directory used to cache downloaded controlled vocabularies, or None to
        disable caching.

    Returns
    -------
    Tuple[Dict[str, Tuple[float, str]], Dict[str, Tuple[float, str]]]
        A tuple with mappings (i) from term accession to modification mass and
        term name, and (ii) from term name to modification mass and term
        accession.

    Raises
    ------
    URLError
        If the controlled vocabulary could not be retrieved from its online
        resource.
    ValueError
        - If an unknown controlled vocabulary identifier is specified.
        - If no mass was specified for a GNO term or its parent terms.
    """
    if cv_id == 'UNIMOD':
        url = 'http://www.unimod.org/obo/unimod.obo'
    elif cv_id in ('MOD', 'RESID'):
        # RESID is not available as a separate entity anymore but is part of
        # PSI-MOD.
        url = ('https://raw.githubusercontent.com/HUPO-PSI/psi-mod-CV/master/'
               'PSI-MOD.obo')
    elif cv_id == 'XLMOD':
        url = ('https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/'
               'cv/XLMOD.obo')
    elif cv_id == 'GNO':
        url = ('https://github.com/glygen-glycan-data/GNOme/releases/latest/'
               'download/GNOme.obo')
    else:
        raise ValueError(f'Unknown controlled vocabulary: {cv_id}')
    # Try to retrieve from the cache.
    cv_cached = _load_from_cache(cache, f'{cv_id}.pkl')
    if cv_cached is not None:
        return cv_cached
    # Read from the online resource if not found in the cache.
    cv_by_accession, cv_by_name, gno_graph = {}, {}, {}
    with urllib.request.urlopen(url) as response:
        if 400 <= response.getcode() < 600:
            raise URLError(f'Failed to retrieve the {cv_id} controlled '
                           f'vocabulary from its online resource')
        for frame in fastobo.load(response):
            term_accession, term_name, term_mass = str(frame.id), None, None
            if isinstance(frame, fastobo.term.TermFrame):
                for clause in frame:
                    if isinstance(clause, fastobo.term.NameClause):
                        term_name = clause.name.strip()
                        if cv_id == 'GNO' and 'molecular weight' in term_name:
                            term_mass = float(
                                term_name[term_name.rindex('weight ') + 7:
                                          term_name.rindex(' Da')])
                    elif (cv_id == 'RESID' and
                          isinstance(clause, fastobo.term.DefClause)):
                        for xref in clause.xrefs:
                            if xref.id.prefix == 'RESID':
                                term_accession = str(xref.id)
                                break
                    elif isinstance(clause, fastobo.term.XrefClause):
                        term_xref = clause.raw_value()
                        if ((cv_id == 'UNIMOD' and
                             'delta_mono_mass' in term_xref) or
                                (cv_id in ('MOD', 'RESID') and
                                 'DiffMono' in term_xref and
                                 not term_xref.endswith('"none"'))):
                            term_mass = float(
                                term_xref[term_xref.index('"') + 1:
                                          term_xref.rindex('"')])
                    elif (cv_id == 'XLMOD' and
                          isinstance(clause, fastobo.term.PropertyValueClause)
                          and (clause.property_value.relation.prefix ==
                               'monoIsotopicMass')):
                        term_mass = float(clause.property_value.value)
                    elif (cv_id == 'GNO' and
                          isinstance(clause, fastobo.term.IsAClause)):
                        gno_graph[term_accession] = term_name, str(clause.term)
                if (term_accession.startswith(f'{cv_id}:') and
                        term_mass is not None):
                    cv_by_accession[term_accession] = term_mass, term_name
                    if term_name is not None:
                        cv_by_name[term_name] = term_mass, term_accession
    if cv_id == 'GNO':
        for term_accession, (term_name, parent_accession) in gno_graph.items():
            if term_accession in cv_by_accession:
                continue
            terms = [(term_accession, term_name)]
            while (parent_accession not in cv_by_accession and
                   parent_accession in gno_graph):
                parent_name, grandparent_accession = \
                    gno_graph[parent_accession]
                terms.append((parent_accession, parent_name))
                parent_accession = grandparent_accession
            if parent_accession not in cv_by_accession:
                raise ValueError(f'No mass found for term {term_accession}'
                                 f' in the GNO controlled vocabulary')
            term_mass = cv_by_accession[parent_accession][0]
            for add_accession, add_name in terms:
                cv_by_accession[add_accession] = term_mass, add_name
                cv_by_name[add_name] = term_mass, add_accession
    # Save to the cache if enabled.
    _store_in_cache(cache, f'{cv_id}.pkl', (cv_by_accession, cv_by_name))
    return cv_by_accession, cv_by_name


def _store_in_cache(cache: Optional[str], filename: str, obj: Any) -> None:
    """
    Store an object in the cache.

    Parameters
    ----------
    cache : Optional[str]
        Directory where the cached objects are stored, or None to disable
        caching.
    filename : str
        Filename of the stored cache object.
    obj : Any
        Object to store in the cache.
    """
    if cache is not None:
        os.makedirs(cache, exist_ok=True)
        with open(os.path.join(cache, filename), 'wb') as f_out:
            pickle.dump(obj, f_out, pickle.HIGHEST_PROTOCOL)


def _load_from_cache(cache: Optional[str], filename: str) -> Optional[Any]:
    """
    Load an object from the cache.

    Parameters
    ----------
    cache : Optional[str]
        Directory where the cached objects are stored, or None to disable
        caching.
    filename : str
        Filename of the stored cache object.

    Returns
    -------
    Any
        The object retrieved from the cache, or None if not found.
    """
    if cache is not None:
        cache_filename = os.path.join(cache, filename)
        if os.path.isfile(cache_filename):
            with open(cache_filename, 'rb') as f_in:
                return pickle.load(f_in)
    return None


def clear_cache(res_ids: Optional[Tuple[str]] = None) -> None:
    """
    Clear the cache of downloaded resources.

    Parameters
    ----------
    res_ids : Optional[Tuple[str]]
        Identifiers of the resources to remove from the cache. If None, all
        known files will be removed.
    """
    # Clear the in-memory cache.
    _import_cv.cache_clear()
    # Clear the on-disk cache.
    if cache_dir is not None:
        if res_ids is None:
            res_ids = ('UNIMOD', 'MOD', 'RESID', 'XLMOD', 'GNO', 'mono')
        for cv_id in res_ids:
            filename = os.path.join(cache_dir, f'{cv_id}.pkl')
            if os.path.isfile(filename):
                os.remove(filename)
