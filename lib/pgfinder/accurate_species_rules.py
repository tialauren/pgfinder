"""
Biochemically accurate species-specific donor/acceptor rules for peptidoglycan.

Rules are evaluated against actual amino acid sequences (e.g., "gm-AQKAA"),
parsed from PGLang structure notation, rather than generic structural names.
"""

import logging
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from pgfinder.logs.logs import LOGGER_NAME

LOGGER = logging.getLogger(LOGGER_NAME)


@dataclass
class PeptideInfo:
    """
    Parsed information about a peptide structure.

    Attributes
    ----------
    sequence : List[str]
        Amino acid sequence (e.g., ['A', 'Q', 'K', 'A', 'A'])
    length : int
        Number of amino acids in stem
    has_mdap_pos3 : bool
        mDAP (diaminopimelic acid) at position 3
    has_lys_pos3 : bool
        Lysine at position 3
    glycine_bridge_length : int
        Number of glycines in bridge (0 if no bridge)
    has_dasp_bridge : bool
        Has D-aspartate bridge
    original_structure : str
        Original structure string from CSV
    """

    sequence: List[str]
    length: int
    has_mdap_pos3: bool
    has_lys_pos3: bool
    glycine_bridge_length: int
    has_dasp_bridge: bool
    original_structure: str


_MODIFICATION_TAG_PATTERN = re.compile(r"^(.*?)( \([^)]*\))$")


def split_modification_tag(structure_base: str) -> Tuple[str, str]:
    """
    Split a structure (with the ``|N`` oligomerization suffix already removed)
    into its (stem_and_glycan_part, modification_tag).

    ``pgfinder.matching.modification_generator`` appends modifications like
    Anhydro-MurNAc/Amidation/Deacetylation as a trailing " (ABBR)" tag (e.g.
    "gm-AEJAA (Anh)"). These act on the glycan or side chains, not on stem
    length/position-3 identity, so donor/acceptor classification and
    structure-trimming should operate on the stem alone - the tag is split
    off here and re-attached by the caller once that's done.

    Examples
    --------
    >>> split_modification_tag("gm-AEJAA (Anh)")
    ('gm-AEJAA', ' (Anh)')
    >>> split_modification_tag("gm-AEJAA")
    ('gm-AEJAA', '')

    Returns
    -------
    Tuple[str, str]
        (structure with any modification tag removed, the tag itself
        including its leading space, or "" if there was none)
    """
    match = _MODIFICATION_TAG_PATTERN.match(structure_base)
    if match:
        return match.group(1), match.group(2)
    return structure_base, ""


def parse_peptide_structure(structure: str) -> Optional[PeptideInfo]:
    """
    Parse peptidoglycan structure to extract peptide sequence and bridge info.

    Examples
    --------
    >>> parse_peptide_structure("gm-AQKAA|1")
    PeptideInfo(sequence=['A','Q','K','A','A'], length=5, has_mdap_pos3=False, has_lys_pos3=True, ...)

    >>> parse_peptide_structure("gm-AQK[GGGGG]AA|1")
    PeptideInfo(sequence=['A','Q','K','A','A'], glycine_bridge_length=5, ...)

    Parameters
    ----------
    structure : str
        Structure string from CSV (e.g., "gm-AQKAA|1" or "gm-AQK[GGGGG]AA|1")

    Returns
    -------
    Optional[PeptideInfo]
        Parsed peptide information, or None if parsing fails
    """
    if not structure or not isinstance(structure, str):
        return None

    # Remove oligomerization state (|1, |2, etc.)
    base_structure = structure.split("|")[0] if "|" in structure else structure

    # Strip any appended modification tag (e.g. "gm-AEJAA (Anh)" -> "gm-AEJAA") -
    # see `split_modification_tag` for why this is safe to discard here.
    base_structure, _ = split_modification_tag(base_structure)

    # Remove disaccharide prefix (gm-, m-, g-, g(-Ac)m-, etc.)
    # Pattern: prefix-SEQUENCE or prefix-SEQUENCE[BRIDGE]
    # Prefix tokens are g/m, each optionally carrying a modification in parentheses
    # (e.g. the de-N-acetylated "g(-Ac)m-" used in the C. diff database).
    match = re.match(r"^(?:[gm](?:\([^)]*\))?)+-([A-Z]+(?:\[[A-Z]+\])?[A-Z]*)$", base_structure, re.IGNORECASE)
    if not match:
        LOGGER.warning(f"Could not parse structure: {structure}")
        return None

    peptide_part = match.group(1)

    # Extract bridge if present (e.g., [GGGGG] or [D])
    bridge_match = re.search(r"\[([A-Z]+)\]", peptide_part)
    bridge_sequence = ""
    if bridge_match:
        bridge_sequence = bridge_match.group(1)
        # Remove bridge from peptide sequence
        peptide_part = re.sub(r"\[[A-Z]+\]", "", peptide_part)

    # Convert to list of amino acids
    sequence = list(peptide_part.upper())

    if not sequence:
        LOGGER.warning(f"Empty sequence after parsing: {structure}")
        return None

    # Analyze the peptide
    length = len(sequence)

    # Check position 3 (index 2) for mDAP or Lys.
    # pgfinder's own mass databases use 'J' for meso-DAP and 'K' for Lysine
    # at this position (e.g. "gm-AEJAA" for E. coli) - they are distinct
    # residues with different masses, not interchangeable.
    has_mdap_pos3 = False
    has_lys_pos3 = False

    if length >= 3:
        pos3 = sequence[2]  # Position 3 (0-indexed as 2)
        if pos3 == "J":
            has_mdap_pos3 = True  # meso-DAP (Gram-negative / DAP-type)
        elif pos3 == "K":
            has_lys_pos3 = True  # Lysine (Gram-positive Lys-type)

    # Analyze bridge
    glycine_bridge_length = 0
    has_dasp_bridge = False

    if bridge_sequence:
        if all(aa == "G" for aa in bridge_sequence):
            glycine_bridge_length = len(bridge_sequence)
        elif bridge_sequence == "D":
            has_dasp_bridge = True

    return PeptideInfo(
        sequence=sequence,
        length=length,
        has_mdap_pos3=has_mdap_pos3,
        has_lys_pos3=has_lys_pos3,
        glycine_bridge_length=glycine_bridge_length,
        has_dasp_bridge=has_dasp_bridge,
        original_structure=structure,
    )


@dataclass
class CustomCrosslinkRule:
    """
    User-defined donor/acceptor rule for species not covered by ``SPECIES_INFO``.

    Attributes
    ----------
    label : str
        Used for the "[label]" tag in output dimer naming.
    donor_pattern : str
        Regex matched against the bridge-stripped stem sequence (e.g. "AEJAA").
    acceptor_pattern : str
        Regex matched against the bridge-stripped stem sequence.
    loses_terminal_ala : bool
        Whether the donor loses a terminal D-Ala during crosslink formation
        (matches the only two lost-residue cases seen across the built-in
        species: D-Ala, or nothing).
    acceptor_bridge_type : str
        One of 'none', 'glycine', 'dasp'. 'none' means bridge composition is
        not checked (matches whatever the acceptor_pattern alone says) - this
        is the right choice for direct (bridge-less) mDAP-type crosslinks.
        'glycine' requires a polyglycine bridge whose length falls within
        [acceptor_min_glycine_bridge, acceptor_max_glycine_bridge] (the S.
        aureus mechanism). 'dasp' requires a single D-Asp/D-Asn bridge residue
        (the E. faecalis mechanism).
    acceptor_min_glycine_bridge : int
        Minimum glycine bridge length required when acceptor_bridge_type is
        'glycine'. Ignored otherwise.
    acceptor_max_glycine_bridge : int
        Maximum glycine bridge length required when acceptor_bridge_type is
        'glycine'. Ignored otherwise.
    """

    label: str
    donor_pattern: str
    acceptor_pattern: str
    loses_terminal_ala: bool
    acceptor_bridge_type: str = "none"
    acceptor_min_glycine_bridge: int = 0
    acceptor_max_glycine_bridge: int = 0


def match_custom_donor(structure: str, rule: CustomCrosslinkRule) -> bool:
    """Check if a structure matches a user-supplied custom donor pattern."""
    info = parse_peptide_structure(structure)
    if info is None:
        return False
    try:
        return re.match(rule.donor_pattern, "".join(info.sequence)) is not None
    except re.error:
        return False


def match_custom_acceptor(structure: str, rule: CustomCrosslinkRule) -> bool:
    """Check if a structure matches a user-supplied custom acceptor pattern,
    including any bridge requirement."""
    info = parse_peptide_structure(structure)
    if info is None:
        return False
    try:
        if re.match(rule.acceptor_pattern, "".join(info.sequence)) is None:
            return False
    except re.error:
        return False

    if rule.acceptor_bridge_type == "glycine":
        return rule.acceptor_min_glycine_bridge <= info.glycine_bridge_length <= rule.acceptor_max_glycine_bridge
    elif rule.acceptor_bridge_type == "dasp":
        return info.has_dasp_bridge
    return True


class AccurateSpeciesRules:
    """
    Biochemically accurate species-specific crosslinking rules.

    Donor/acceptor validity is determined from the parsed amino acid
    sequence and position-specific checks, documented per-species in
    ``SPECIES_INFO`` and the methods below.
    """

    # Species database with descriptive information
    SPECIES_INFO = {
        "ecoli": {
            "name": "Escherichia coli",
            "description": "Gram-negative DAP-type with D,D and L,D-transpeptidases",
            "crosslink_types": ["4-3", "3-3"],
        },
        "fusobacterium": {
            "name": "Fusobacterium",
            "description": "Gram-negative DAP-type",
            "crosslink_types": ["4-3", "3-3"],
        },
        "cdiff": {
            "name": "Clostridioides difficile",
            "description": "Gram-positive DAP-type",
            "crosslink_types": ["4-3", "3-3"],
        },
        "bsubtilis": {
            "name": "Bacillus subtilis",
            "description": "Strain-dependent (mDAP-type common, Lys-type rare)",
            "crosslink_types": ["4-3", "3-3"],
        },
        "saureus": {
            "name": "Staphylococcus aureus",
            "description": "Lys-type with pentaglycine bridge",
            "crosslink_types": ["4-3-bridge"],
        },
        "efaecalis": {
            "name": "Enterococcus faecalis",
            "description": "Lys-type with D-Asp bridge and 3-3 crosslinks",
            "crosslink_types": ["4-3-bridge", "3-3"],
        },
    }

    @classmethod
    def get_species_info(cls, species_code: str) -> Dict:
        """Get descriptive information about a species."""
        species_code = species_code.lower()
        if species_code not in cls.SPECIES_INFO:
            available = ", ".join(cls.SPECIES_INFO.keys())
            raise ValueError(f"Unknown species: '{species_code}'. Available: {available}")
        return cls.SPECIES_INFO[species_code]

    @classmethod
    def list_available_species(cls) -> List[Tuple[str, str]]:
        """Return list of available species with their full names."""
        return [(code, data["name"]) for code, data in cls.SPECIES_INFO.items()]

    @classmethod
    def is_valid_donor(cls, structure: str, species_code: str, crosslink_type: str) -> bool:
        """
        Check if a structure can act as a donor for specified crosslink type.

        Parameters
        ----------
        structure : str
            Structure string (e.g., "gm-AQKAA|1")
        species_code : str
            Species identifier
        crosslink_type : str
            Type of crosslink ('4-3', '3-3', '4-3-bridge')

        Returns
        -------
        bool
            True if structure can be a donor
        """
        info = parse_peptide_structure(structure)
        if info is None:
            return False

        species_code = species_code.lower()

        # E. coli, Fusobacterium, C. diff (Gram-negative/DAP-type)
        if species_code in ["ecoli", "fusobacterium", "cdiff"]:
            if crosslink_type == "4-3":
                return info.length == 5 and info.has_mdap_pos3 and info.sequence[3] == "A" and info.sequence[4] == "A"
            elif crosslink_type == "3-3":
                # L,D-transpeptidase cleaves the mDAP3–D-Ala4 bond specifically;
                # position 4 must be D-Ala (A) for the enzyme to act on it.
                return info.length == 4 and info.has_mdap_pos3 and info.sequence[3] == "A"

        # B. subtilis (handles both mDAP and Lys types)
        elif species_code == "bsubtilis":
            if crosslink_type == "4-3":
                is_dap = info.length == 5 and info.has_mdap_pos3 and info.sequence[3] == "A" and info.sequence[4] == "A"
                is_lys = info.length == 5 and info.has_lys_pos3 and info.sequence[3] == "A" and info.sequence[4] == "A"
                return is_dap or is_lys
            elif crosslink_type == "3-3":
                # L,D-transpeptidase cleaves the pos3–D-Ala4 bond; position 4 must be D-Ala.
                is_dap = info.length == 4 and info.has_mdap_pos3 and info.sequence[3] == "A"
                is_lys = info.length == 4 and info.has_lys_pos3 and info.sequence[3] == "A"
                return is_dap or is_lys

        # S. aureus (pentaglycine bridge)
        elif species_code == "saureus":
            if crosslink_type == "4-3-bridge":
                return info.length == 5 and info.has_lys_pos3 and info.sequence[3] == "A" and info.sequence[4] == "A"

        # E. faecalis (D-Asp bridge)
        elif species_code == "efaecalis":
            if crosslink_type == "4-3-bridge":
                return (
                    info.length == 5
                    and info.has_lys_pos3
                    and info.sequence[3] == "A"
                    and info.sequence[4] == "A"
                    and not info.has_dasp_bridge
                )
            elif crosslink_type == "3-3":
                # L,D-transpeptidase cleaves the Lys3–D-Ala4 bond; position 4 must be D-Ala.
                return info.length == 4 and info.has_lys_pos3 and info.sequence[3] == "A" and not info.has_dasp_bridge

        return False

    @classmethod
    def is_valid_acceptor(
        cls, structure: str, species_code: str, crosslink_type: str, strict_mode: bool = True
    ) -> bool:
        """
        Check if a structure can act as an acceptor for specified crosslink type.

        Parameters
        ----------
        structure : str
            Structure string (e.g., "gm-AQK[GGGGG]AA|1")
        species_code : str
            Species identifier
        crosslink_type : str
            Type of crosslink ('4-3', '3-3', '4-3-bridge')
        strict_mode : bool
            If True: Enforce exact literature rules (default)
            If False: Allow variants for novel discovery

        Returns
        -------
        bool
            True if structure can be an acceptor
        """
        info = parse_peptide_structure(structure)
        if info is None:
            return False

        species_code = species_code.lower()

        # E. coli, Fusobacterium, C. diff
        if species_code in ["ecoli", "fusobacterium", "cdiff"]:
            if crosslink_type == "4-3":
                return info.has_mdap_pos3 and info.length >= 3
            elif crosslink_type == "3-3":
                return info.has_mdap_pos3 and info.length >= 3

        # B. subtilis
        elif species_code == "bsubtilis":
            if crosslink_type == "4-3":
                return (info.has_mdap_pos3 or info.has_lys_pos3) and info.length >= 3
            elif crosslink_type == "3-3":
                return (info.has_mdap_pos3 or info.has_lys_pos3) and info.length >= 3

        # S. aureus (CRITICAL: Bridge length rules)
        elif species_code == "saureus":
            if crosslink_type == "4-3-bridge":
                if strict_mode:
                    # STRICT: Exactly 5 glycines (literature standard)
                    return info.has_lys_pos3 and info.length >= 3 and info.glycine_bridge_length == 5
                else:
                    # PERMISSIVE: 1-5 glycines (for novel discovery)
                    return info.has_lys_pos3 and info.length >= 3 and 1 <= info.glycine_bridge_length <= 5

        # E. faecalis
        elif species_code == "efaecalis":
            if crosslink_type == "4-3-bridge":
                return info.has_lys_pos3 and info.length >= 3 and info.has_dasp_bridge
            elif crosslink_type == "3-3":
                return info.length == 4 and info.has_lys_pos3 and not info.has_dasp_bridge

        return False

    @classmethod
    def get_crosslink_types(cls, species_code: str) -> List[str]:
        """
        Get all crosslink types supported by a species.

        Returns
        -------
        List[str]
            List of crosslink types (e.g., ['4-3', '3-3'])
        """
        info = cls.get_species_info(species_code)
        return info["crosslink_types"]

    @classmethod
    def get_lost_residues(cls, crosslink_type: str) -> List[str]:
        """
        Get residues lost during crosslink formation.

        Parameters
        ----------
        crosslink_type : str
            Type of crosslink

        Returns
        -------
        List[str]
            List of amino acids lost (e.g., ['A'] for terminal D-Ala)
        """
        # 4-3/4-3-bridge: D,D-transpeptidase cleaves the D-Ala4-D-Ala5 bond of a
        # pentapeptide donor, releasing free D-Ala5.
        # 3-3: L,D-transpeptidase cleaves the mDAP3-D-Ala4 bond of a tetrapeptide
        # donor, releasing free D-Ala4 - same lost residue, one position earlier.
        if crosslink_type in ["4-3", "4-3-bridge", "3-3"]:
            return ["A"]
        return []
