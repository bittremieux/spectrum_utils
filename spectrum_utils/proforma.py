import functools
import re
from typing import Dict, Optional, Tuple, Union

import pronto
import requests
from pyteomics.auxiliary.structures import PyteomicsError
try:
    from pyteomics import cmass as mass
except ImportError:
    from pyteomics import mass

import spectrum_utils.spectrum as sus


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
        elif peptide[i] in sus.aa_mass.keys():
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
    # Modification specified by name or mass only (no / no valid prefix).
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
                return _cv_lookup(f'M:{modification}')
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
                raise NotImplementedError('Isotopes in molecular formulas are '
                                          'not supported')
        # Calculate mass from glycan composition.
        elif modification.title().startswith('Glycan:'):
            # FIXME: Cache web resource.
            monosaccharides = {
                term['name']: float(term['has_monoisotopic_mass'])
                for term in (requests.get(
                                'https://raw.githubusercontent.com/HUPO-PSI/'
                                'ProForma/master/monosaccharides/'
                                'mono.obo.json')
                             .json()['terms'].values())}
            return sum([monosaccharides[gly[1]] * int(gly[2])
                        for gly in re.finditer(
                            r'([a-zA-Z]+)(\d+)',
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
    KeyError
        If the term was not found in its controlled vocabulary.
    ValueError
        - If an unknown controlled vocabulary identifier is given.
        - If no mass was specified for a GNO term or its parent terms.
    """
    cv_id, term_lookup = lookup.split(':', 1)
    # Name lookups are specified using single-letter abbrevations, ID lookups
    # use the full controlled vocabulary name.
    lookup_by_id = len(cv_id) > 1
    cv_id = {'U': 'UNIMOD', 'M': 'MOD', 'R': 'RESID', 'X': 'XLMOD',
             'G': 'GNO'}.get(cv_id, cv_id)
    cv, cv_name_to_id = _import_cv(cv_id)
    try:
        # Look up the term by its ID.
        if lookup_by_id:
            term_id = f'{cv_id}:{term_lookup}'
            # Translate RESID ids (RESID is not available as a separate entity
            # anymore but is part of PSI-MOD).
            if cv_id == 'RESID' and term_id in cv_name_to_id:
                term_id = cv_name_to_id[term_id]
        # Translate the term name to its ID.
        elif term_lookup in cv_name_to_id:
            term_id = cv_name_to_id[term_lookup]
        # Look up by name, but term name not found in the CV.
        else:
            raise KeyError(f'Term {lookup} not found in the {cv_id} '
                           f'controlled vocabulary by name')
        # Retrieve the term from the CV.
        term = cv.get_term(term_id)
    except KeyError:
        raise KeyError(f'Term {lookup} not found in the {cv_id} '
                       f'controlled vocabulary '
                       f'{"by ID" if lookup_by_id else "by name"}')
    # Extract the mass difference from the term based on the originating CV.
    if cv_id == 'UNIMOD':
        for xref in term.xrefs:
            if xref.id == 'delta_mono_mass':
                return float(xref.description)
    elif cv_id in ('MOD', 'RESID'):
        for xref in term.xrefs:
            if xref.id == 'DiffMono:':
                return float(xref.description)
    elif cv_id == 'XLMOD':
        for annotation in term.annotations:
            if (isinstance(annotation, pronto.pv.LiteralPropertyValue) and
                    annotation.property == 'monoIsotopicMass:'):
                return float(annotation.literal)
    elif cv_id == 'GNO':
        if 'molecular weight' not in term.name:
            for superclass in term.superclasses():
                if 'molecular weight' in superclass.name:
                    term = superclass
                    break
            else:
                raise ValueError(f'No mass found for term {lookup} in the '
                                 f'GNO controlled vocabulary')
        return float(term.name[term.name.rindex('weight ') + 7:
                               term.name.rindex(' Da')])
    raise KeyError(f'Term {lookup} not found in the {cv_id} controlled '
                   f'vocabulary {"by ID" if lookup_by_id else "by name"}')


@functools.lru_cache
def _import_cv(cv_id: str) -> Tuple[pronto.Ontology, Dict]:
    """
    Import a ProForma controlled vocabulary from its online resource.

    Parameters
    ----------
    cv_id : str
        The controlled vocabulary identifier.

    Returns
    -------
    Tuple[pronto.Ontology, Dict]
        A tuple of (i) the controlled vocabulary, and (ii) a mapping from term
        names to term IDs.

    Raises
    ------
    ValueError
        If an unknown controlled vocabulary identifier is given.
    """
    # FIXME: Cache web resource.
    if cv_id == 'UNIMOD':
        cv = pronto.Ontology('http://www.unimod.org/obo/unimod.obo')
    elif cv_id in ('MOD', 'XLMOD', 'GNO'):
        cv = pronto.Ontology.from_obo_library(f'{cv_id.lower()}.obo')
    elif cv_id == 'RESID':
        # RESID is not available as a separate entity anymore but is part of
        # PSI-MOD.
        cv = pronto.Ontology.from_obo_library('mod.obo')
    else:
        raise ValueError(f'Unknown controlled vocabulary: {cv_id}')
    name_to_id = {term.name.strip(): term.id for term in cv.terms()}
    # Translate from RESID identifiers to PSI-MOD identifiers.
    if cv_id == 'RESID':
        for term in cv.terms():
            for xref in term.definition.xrefs:
                if xref.id.startswith('RESID'):
                    name_to_id[xref.id] = term.id
    return cv, name_to_id
