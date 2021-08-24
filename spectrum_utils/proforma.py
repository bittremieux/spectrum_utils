import collections
import copy
import enum
import functools
import json
import os
import pickle
import re
import urllib.request
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple, Union
from urllib.error import URLError

import fastobo
import lark
from pyteomics.auxiliary.structures import PyteomicsError
try:
    from pyteomics import cmass as mass
except ImportError:
    from pyteomics import mass

from spectrum_utils.spectrum import aa_mass


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
    controlled_vocabulary: str
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

    def __init__(self):
        super().__init__()
        self.sequence, self.modifications = [], []
        self.global_modifications = collections.defaultdict(list)

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
        position, *mods = tree
        for mod in mods:
            mod.position = position
            self.modifications.append(mod)

    def mod_range_pos(self, tree) -> Tuple[int, int]:
        return len(self.sequence) - len(tree), len(self.sequence) - 1

    def mod_name(self, tree) -> CvEntry:
        cv, name = tree if len(tree) == 2 else (None, tree[0])
        return CvEntry(controlled_vocabulary=cv, name=name)

    def CV_ABBREV_OPT(self, token) -> str:
        return self.CV_ABBREV(token)

    def CV_ABBREV(self, token) -> str:
        return {'U': 'UNIMOD', 'M': 'MOD', 'R': 'RESID', 'X': 'XLMOD',
                'G': 'GNO'}[token.value.upper()]

    def mod_accession(self, tree) -> CvEntry:
        return CvEntry(controlled_vocabulary=tree[0], accession=tree[1])

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
        return Glycan(composition=tree)

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


def parse(proforma: str) -> List[Proteoform]:
    """
    Parse a ProForma-encoded string.

    The string is parsed by building an abstract syntax tree to interpret the
    ProForma language. Next, modifications are localized to residue positions
    and translated into Python objects.

    Parameters
    ----------
    proforma : str
        The ProForma string.

    Returns
    -------
    List[Proteoform]
        A list of proteoforms with their modifications (and additional)
        information specified by the ProForma string.
    """
    with open('spectrum_utils/proforma.ebnf') as f_in:
        parser = lark.Lark(f_in.read(), start='proforma', parser='earley',
                           lexer='dynamic_complete')
        # TODO: Interpret modifications encoded by a CV.
        return ProFormaTransformer().transform(parser.parse(proforma))


def parse_regex(proforma: str) -> Tuple[str, Dict[Union[int, str], float]]:
    """
    Parse a ProForma-encoded string.

    Parameters
    ----------
    proforma : str
        The ProForma string.

    Returns
    -------
    TODO
    Tuple[str, Dict[Union[int, str], float]]
        (1) The base peptide string without modifications, and (2) a mapping of
        modification positions and mass differences. Modification positions are
        specified as their amino acid index in the peptide (0-based), 'N-term',
        or 'C-term'.

    Raises
    ------
    URLError
        If a controlled vocabulary could not be retrieved from its online
        resource.
    KeyError
        If a specified modification could not be found in the corresponding
        controlled vocabulary.
    NotImplementedError
        - Inter-chain crosslinks are not supported (specified with "\\").
        - Isotopes in a molecular formula are not supported.
        - Modification positions ranges are not supported.
    ValueError
        - If the same label is used for multiple preferred positions.
        - If an unknown character (amino acid) is encountered.
        - If a modification annotation couldn't be successfully parsed.
        - If no mass was specified for a GNO term or its parent terms.
    """
    peptide_base, modifications = [], {}
    preferred_positions, global_modifications = set(), {}
    if '\\' in peptide:
        raise NotImplementedError('Inter-chain crosslinks are not supported')
    # Extract global modifications.
    if peptide.startswith('<'):
        end_index = peptide.index('>')
        modification, residues = peptide[1:end_index].split('@')
        mod_mass = _parse_modification(modification[1:-1])
        if mod_mass is not None:
            global_modifications = {residue.upper(): mod_mass
                                    for residue in residues.split(',')}
        peptide = peptide[end_index + 1:]
    # Extract modifications from the peptide string.
    i = 0
    while i < len(peptide):
        # Parse modification on residue.
        if peptide[i] == '[':
            pos, j, left_bracket_count, right_bracket_count = None, i + 1, 1, 0
            while right_bracket_count < left_bracket_count:
                j += 1
                left_bracket_count += peptide[j] == '['
                right_bracket_count += peptide[j] == ']'
            if len(peptide_base) == 0:
                # Modification as first element followed by a hyphen indicates
                # an N-terminal modification.
                if peptide[j + 1] == '-':
                    pos = 'N-term'
                # Modification as first element followed by "?" (with optional
                # cardinality) indicates unknown modification positions.
                # These will be ignored because they can't be assigned to
                # specific positions in the peptide.
                elif peptide[j + 1] in '?^[':
                    j += 1
                    while peptide[j] != '?':
                        j += 1
            # Modification as final element preceded by a hyphen indicates a
            # C-terminal modification.
            elif i > 0 and peptide[i - 1] == '-':
                pos = 'C-term'
            # Modification on a specific residue.
            else:
                pos = len(peptide_base) - 1
            # Process modifications with a valid positions that don't signify
            # secondary locations ("#" pointing to another preferred location).
            if pos is not None and peptide[i + 1] != '#':
                modification = peptide[i + 1:j]
                # Extract the preferred position tag (if applicable).
                pref_pos = re.match(
                    r'^(.+)(#[a-zA-z0-9]+)(\(([0-9]*[.])?[0-9]+\))?$',
                    modification)
                if pref_pos is not None:
                    modification = pref_pos[1]
                    if pref_pos[2] in preferred_positions:
                        raise ValueError(
                            f'Multiple preferred positions specified with '
                            f'label "{pref_pos[2]}"')
                    else:
                        preferred_positions.add(pref_pos[2])
                # If multiple modification options are specified, just take the
                # first one.
                for mod in modification.split('|'):
                    mod_mass = _parse_modification(mod)
                    if mod_mass is not None:
                        modifications[pos] = mod_mass
                        break
            i = j + 1
        # Skip N-term / C-term connectors.
        elif peptide[i] == '-':
            i += 1
        # Ignore labile modifications.
        elif peptide[i] == '{':
            j, left_bracket_count, right_bracket_count = i + 1, 1, 0
            while right_bracket_count < left_bracket_count:
                j += 1
                left_bracket_count += peptide[j] == '{'
                right_bracket_count += peptide[j] == '}'
            i = j + 1
        # Modification position ranges are currently not supported because they
        # can't be assigned to a specific residue.
        elif peptide[i] == '(':
            raise NotImplementedError('Modification position ranges are not '
                                      'supported')
        # Standard amino acid.
        elif peptide[i] in aa_mass.keys():
            peptide_base.append(peptide[i])
            i += 1
        # Give an error for amino acids whose mass can't be computed.
        else:
            raise ValueError(f'Unknown/unsupported amino acid "{peptide[i]}"')
    peptide_base = ''.join(peptide_base).upper()
    # Apply any global modifications.
    for i, residue in enumerate(peptide_base):
        if residue in global_modifications:
            modifications[i] = global_modifications[residue]

    return peptide_base, modifications


def _parse_modification(modification: str) -> Optional[float]:
    """
    Convert a ProForma modification string to the corresponding mass
    difference.

    Parameters
    ----------
    modification : str
        The ProForma modification string.

    Returns
    -------
    Optional[float]
        The modification mass difference, or `None` for info tags.

    Raises
    ------
    URLError
        If the controlled vocabulary could not be retrieved from its online
        resource.
    KeyError
        If the modification could not be found in the corresponding controlled
        vocabulary.
    NotImplementedError
        Isotopes in a molecular formula (not supported by Pyteomics mass
        calculation).
    ValueError
        - If the modification couldn't be successfully parsed.
        - If no mass was specified for a GNO term or its parent terms.
    """
    # Info tag.
    if modification.upper().startswith('INFO:'):
        return None
    # Modification specified by name or mass only (no / invalid prefix).
    if (':' not in modification or
            modification[:modification.index(':')].upper() not in (
                    'U', 'UNIMOD', 'M', 'MOD', 'R', 'RESID', 'X', 'XLMOD',
                    'G', 'GNO', 'OBS', 'FORMULA', 'GLYCAN')):
        # Numerical mass difference.
        if modification[0] in '+-':
            return float(modification)
        # Retrieve from Unimod or PSI-MOD by name.
        else:
            try:
                return _cv_lookup(f'U:{modification}')
            except KeyError:
                try:
                    return _cv_lookup(f'M:{modification}')
                except KeyError:
                    raise KeyError(f'Term {modification} not found in the '
                                   f'UNIMOD or PSI-MOD reference controlled '
                                   f'vocabularies')
    else:
        # Numerical mass difference (controlled vocabulary or "Obs" prefix).
        if (modification.title().startswith(('U:', 'M:', 'R:', 'X:', 'G:',
                                             'Obs:'))
                and modification[modification.index(':') + 1] in '+-'):
            return float(modification[modification.index(':') + 1:])
        # Retrieve from controlled vocabulary by name or accession.
        elif modification.upper().startswith(
                ('U:', 'UNIMOD:', 'M:', 'MOD:', 'R:', 'RESID:',
                 'X:', 'XLMOD:', 'G:', 'GNO:')):
            return _cv_lookup(modification)
        # Calculate mass from molecular formula.
        elif modification.title().startswith('Formula:'):
            try:
                # Make sure there are no spaces in the molecular formula
                # (Pyteomics mass calculation can't handle those).
                return mass.calculate_mass(
                    formula=(modification[modification.index(':') + 1:]
                             .replace(' ', '')))
            except PyteomicsError:
                # TODO: Add isotope support to Pyteomics?
                raise NotImplementedError('Isotopes in molecular formulas are '
                                          'not supported')
        # Calculate mass from glycan composition.
        elif modification.title().startswith('Glycan:'):
            mono = _load_from_cache(cache_dir, 'mono.pkl')
            if mono is None:
                with urllib.request.urlopen(
                        'https://raw.githubusercontent.com/HUPO-PSI/ProForma/'
                        'master/monosaccharides/mono.obo.json') as response:
                    if 400 <= response.getcode() < 600:
                        raise URLError('Failed to retrieve the monosaccharide '
                                       'definitions from its online resource')
                    mono = json.loads(response.read().decode(
                        response.info().get_param('charset') or 'utf-8'))
                mono = {term['name']: float(term['has_monoisotopic_mass'])
                        for term in mono['terms'].values()}
                _store_in_cache(cache_dir, 'mono.pkl', mono)
            return sum([mono[m[1]] * int(m[2]) for m in re.finditer(
                fr'({"|".join([re.escape(m) for m in mono.keys()])})(\d+)',
                modification[modification.index(':') + 1:])])


def _cv_lookup(lookup: str) -> Optional[float]:
    """
    Retrieve a term from a controlled vocabulary by identifier or name.

    Parameters
    ----------
    lookup : str
        The ProForma-formatted term to look up in the controlled vocabulary.

    Returns
    -------
    Optional[float]
        The mass difference of the modification as specified in the controlled
        vocabulary, or None if not found.

    Raises
    ------
    URLError
        If the controlled vocabulary could not be retrieved from its online
        resource.
    KeyError
        If the term was not found in its controlled vocabulary.
    ValueError
        - If an unknown controlled vocabulary identifier is given.
        - If no mass was specified for a GNO term or its parent terms.
    """
    cv_id, term_lookup = lookup.split(':', 1)
    # Name lookups are specified using single-letter abbreviations, ID lookups
    # use the full controlled vocabulary name.
    lookup_by_id = len(cv_id) > 1
    cv_id = {'U': 'UNIMOD', 'M': 'MOD', 'R': 'RESID', 'X': 'XLMOD',
             'G': 'GNO'}.get(cv_id, cv_id)
    cv_from_id, cv_from_name = _import_cv(cv_id, cache_dir)
    try:
        return (cv_from_id[f'{cv_id}:{term_lookup}'] if lookup_by_id else
                cv_from_name[term_lookup])
    except KeyError:
        raise KeyError(f'Term {lookup} not found in the {cv_id} '
                       f'controlled vocabulary '
                       f'{"by ID" if lookup_by_id else "by name"}')


@functools.lru_cache
def _import_cv(cv_id: str, cache: Optional[str]) \
        -> Tuple[Dict[str, float], Dict[str, float]]:
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
    Tuple[Dict[str, float], Dict[str, float]]
        A tuple with mappings (i) from term ID to modification mass, and
        (ii) from term name to modification mass.

    Raises
    ------
    URLError
        If the controlled vocabulary could not be retrieved from its online
        resource.
    ValueError
        - If an unknown controlled vocabulary identifier is given.
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
    cv_by_id, cv_by_name, gno_graph = {}, {}, {}
    with urllib.request.urlopen(url) as response:
        if 400 <= response.getcode() < 600:
            raise URLError(f'Failed to retrieve the {cv_id} controlled '
                           f'vocabulary from its online resource')
        for frame in fastobo.load(response):
            term_id, term_name, term_mass = str(frame.id), None, None
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
                                term_id = str(xref.id)
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
                        gno_graph[term_id] = str(clause.term), term_name
                if term_id.startswith(f'{cv_id}:') and term_mass is not None:
                    cv_by_id[term_id] = term_mass
                    if term_name is not None:
                        cv_by_name[term_name] = term_mass
    if cv_id == 'GNO':
        for term_id, (parent_id, term_name) in gno_graph.items():
            if term_id not in cv_by_id:
                terms = [(term_id, term_name)]
                while parent_id not in cv_by_id and parent_id in gno_graph:
                    parent_id_new, parent_name = gno_graph[parent_id]
                    terms.append((parent_id, parent_name))
                    parent_id = parent_id_new
                if parent_id in cv_by_id:
                    term_mass = cv_by_id[parent_id]
                    for add_id, add_name in terms:
                        cv_by_id[add_id] = cv_by_name[add_name] = term_mass
                else:
                    raise ValueError(f'No mass found for term {term_id} in '
                                     f'the GNO controlled vocabulary')
    # Save to the cache if enabled.
    _store_in_cache(cache, f'{cv_id}.pkl', (cv_by_id, cv_by_name))
    return cv_by_id, cv_by_name


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
