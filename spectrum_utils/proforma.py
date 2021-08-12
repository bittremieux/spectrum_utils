import functools
import json
import os
import pickle
import re
import urllib.request
from typing import Dict, Optional, Tuple, Union
from urllib.error import URLError

import fastobo
from pyteomics.auxiliary.structures import PyteomicsError
try:
    from pyteomics import cmass as mass
except ImportError:
    from pyteomics import mass

from spectrum_utils.spectrum import aa_mass


# Set to None to disable caching.
cache_dir = os.path.join(os.path.expanduser('~'), '.cache', 'spectrum_utils')


def parse(peptide: str) -> Tuple[str, Dict[Union[int, str], float]]:
    """
    Parse a ProForma-encoded peptide string.

    Parameters
    ----------
    peptide : str
        The ProForma peptide string.

    Returns
    -------
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
            j, left_bracket_count, right_bracket_count = i + 1, 1, 0
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
                    pos = None
                else:
                    pos = None      # This shouldn't occur, so just ignore.
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
                # TODO
                raise NotImplementedError('Isotopes in molecular formulas are '
                                          'not supported')
        # Calculate mass from glycan composition.
        elif modification.title().startswith('Glycan:'):
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
                return sum([mono[m[1]] * int(m[2]) for m in re.finditer(
                    fr'({"|".join([re.escape(m) for m in mono.keys()])})(\d+)',
                    modification[modification.index(':') + 1:])])
    raise ValueError(f'Unknown ProForma modification: {modification}')


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
    if cache is not None:
        cache_filename = os.path.join(cache, f'{cv_id}.pkl')
        if os.path.isfile(cache_filename):
            with open(cache_filename, 'rb') as f_in:
                return pickle.load(f_in)
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
    if cache is not None:
        os.makedirs(cache, exist_ok=True)
        with open(os.path.join(cache, f'{cv_id}.pkl'), 'wb') as f_out:
            pickle.dump((cv_by_id, cv_by_name), f_out, pickle.HIGHEST_PROTOCOL)
    return cv_by_id, cv_by_name


def clear_cache(cv_ids: Optional[Tuple[str]] = None) -> None:
    """
    Clear the cache of downloaded controlled vocabularies.

    Parameters
    ----------
    cv_ids : Optional[Tuple[str]]
        Identifiers of the controlled vocabularies to remove from the cache.
        If None, all known files will be removed.
    """
    # Clear the in-memory cache.
    _import_cv.cache_clear()
    # Clear the on-disk cache.
    if cache_dir is not None:
        if cv_ids is None:
            cv_ids = ('UNIMOD', 'MOD', 'RESID', 'XLMOD', 'GNO')
        for cv_id in cv_ids:
            filename = os.path.join(cache_dir, f'{cv_id}.pkl')
            if os.path.isfile(filename):
                os.remove(filename)
