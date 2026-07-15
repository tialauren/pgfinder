"""
HCD fragment ion prediction for muropeptide monomers.

Residue-level mass arithmetic — no RDKit required.
Fragmentation rules mirror hcd_rules.kdl in TheLostLambda/smithereens.

Bond cleavage model
-------------------
For a fragment containing k residues (R1..Rk), the **neutral mass** is:

    neutral = sum(free_residue_mass(Ri) for i in 1..k)
              - (k - 1) * H2O          # condensation bonds within the fragment
              + terminus_correction     # depends on which end was cut

Terminus corrections (from hcd_rules.kdl):
    b / B terminus  lost="OHH"          → correction = -H2O
    y / Y / C       lost="H" gained="H" → correction = 0   (net zero)
    Z               lost="OHH"          → correction = -H2O

Then add PROTON per unit charge for [M+nH]^n+ ions.

Reduction (AUTO_MODS "Red")
---------------------------
pgfinder mass databases assume NaBH4 reduction of the reducing-end saccharide
(+H2, +2.01565 Da).  Fragments that CONTAIN the reducing end of MurNAc carry
this offset; fragments that don't (B-type glycan ions, b-lac ions, all peptide
ions) do not.
"""

import re
from typing import List, Optional

import pandas as pd

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

PROTON = 1.00727646677
H2O = 18.01056469
H2_REDUCTION = 2.01565007  # NaBH4 reduction: gained "H2"
CO = 27.99491462  # for immonium ions (b+y single-residue - CO)

# ---------------------------------------------------------------------------
# Monoisotopic atomic masses
# ---------------------------------------------------------------------------

_ATOMS: dict = {
    "C": 12.00000000000,
    "H": 1.00782503207,
    "N": 14.00307400480,
    "O": 15.99491461956,
    "S": 31.97207100000,
    "P": 30.97376163000,
}

_FORMULA_RE = re.compile(r"([A-Z][a-z]?)(\d*)")


def _mass(formula: str) -> float:
    """Monoisotopic mass from molecular formula string, e.g. 'C8H15NO6'."""
    total = 0.0
    for element, count in _FORMULA_RE.findall(formula):
        if element in _ATOMS:
            total += _ATOMS[element] * (int(count) if count else 1)
    return total


# ---------------------------------------------------------------------------
# Residue free-molecule masses
# Compositions from polymer_database.kdl (TheLostLambda/smithereens).
# These are the FULL monoisotopic masses before condensation (each bond
# removes one H2O via the `lost "H2O"` bond definition in the KDL).
# ---------------------------------------------------------------------------

_RESIDUE_MASS: dict = {
    # Monosaccharides
    "g": _mass("C8H15NO6"),  # GlcNAc  221.08994
    "m": _mass("C11H19NO8"),  # MurNAc  293.11107
    # Amino acids
    "A": _mass("C3H7NO2"),  # Ala      89.04768
    "B": _mass("C4H10N2O2"),  # DAB     118.07423
    "C": _mass("C3H7NO2S"),  # Cys     121.01975
    "D": _mass("C4H7NO4"),  # Asp     133.03751
    "E": _mass("C5H9NO4"),  # Glu     147.05316
    "F": _mass("C9H11NO2"),  # Phe     165.07898
    "G": _mass("C2H5NO2"),  # Gly      75.03203
    "H": _mass("C6H9N3O2"),  # His     155.06948
    "I": _mass("C6H13NO2"),  # Ile     131.09463
    "J": _mass("C7H14N2O4"),  # mDAP   190.09535
    "K": _mass("C6H14N2O2"),  # Lys    146.10553
    "L": _mass("C6H13NO2"),  # Leu    131.09463
    "M": _mass("C5H11NO2S"),  # Met    149.05105
    "N": _mass("C4H8N2O3"),  # Asn    132.05349
    "O": _mass("C5H12N2O2"),  # Orn    132.08988
    "P": _mass("C5H9NO2"),  # Pro    115.06333
    "Q": _mass("C5H10N2O3"),  # Gln    146.06914
    "R": _mass("C6H14N4O2"),  # Arg    174.11168
    "S": _mass("C3H7NO3"),  # Ser    105.04259
    "T": _mass("C4H9NO3"),  # Thr    119.05824
    "U": _mass("C4H9NO3"),  # HomoSer119.05824
    "V": _mass("C5H11NO2"),  # Val    117.07898
    "W": _mass("C11H12N2O2"),  # Trp    204.08988
    "Y": _mass("C9H11NO3"),  # Tyr    181.07389
    "Z": _mass("C5H9NO5"),  # HydroxyGlu 163.04808
    "mL": _mass("C6H12N2O4S"),  # Lanthionine208.05175
}

# Lactoyl group mass: the C3H4O2 moiety in the MurNAc lactyl ether.
# b-lac terminus lost="C3H4O2OHH": loses C3H4O2 AND H2O.
# y-lac terminus lost="H" gained="C3H4O2H": gains C3H4O2 (net).
_LACTOYL = _mass("C3H4O2")  # 72.02113

# Secondary losses from GlcNAc (g) freed fragments — [M+H]+ values given in
# comments. Each is the neutral mass of the freed fragment.
# From hcd_rules.kdl residue "g" block.
_GLC_FREED: list = [
    _mass("C8H11NO4"),  # [M+H]+ 186.076
    _mass("C8H9NO3"),  # [M+H]+ 168.066
    _mass("C6H9NO3"),  # [M+H]+ 144.066
    _mass("C7H7NO2"),  # [M+H]+ 138.055
    _mass("C6H7NO2"),  # [M+H]+ 126.055
]
_GLC_FREED_LABELS = ["Glc-1", "Glc-2", "Glc-3a", "Glc-3b", "Glc-4"]

# Secondary losses from MurNAc (m) freed fragments.
_MUR_FREED: list = [_mass("C7H7NO2")]  # [M+H]+ 138.055
_MUR_FREED_LABELS = ["Mur-1"]

# ---------------------------------------------------------------------------
# Structure string parsing
# ---------------------------------------------------------------------------

# Matches standard monomer: e.g. "gm-AEJA|1", "gm-AEJ[GGGGG]AA|1",
# "gm-AEJAA (Anh)|1".  Does not handle dimers (those contain '=').
_MONOMER_RE = re.compile(
    r"^(?P<glycan>[gm]+)"  # glycan prefix: g, m, gm, gmgm ...
    r"-(?P<stem>[A-Z]+(?:\[[A-Z]+\][A-Z]*)?)"  # stem peptide, with optional bridge
    r"(?P<mod> \([^)]+\))?"  # optional modification tag
    r"\|(?P<oligo>\d+)$"  # oligomerisation state
)


def _parse_glycan(glycan_str: str) -> List[str]:
    """Split glycan prefix string into residue list, e.g. 'gm' -> ['g', 'm']."""
    return list(glycan_str)


def _parse_stem(stem_str: str) -> List[str]:
    """
    Split stem peptide string into residue list.

    Handles bridge notation: 'AEJ[GGGGG]AA' -> ['A','E','J','G','G','G','G','G','A','A']
    The bridge is inlined at the crosslink position for fragmentation purposes.
    """
    # Replace bridge notation [XYZ] with just XYZ inline
    flattened = re.sub(r"\[([A-Z]+)\]", r"\1", stem_str)
    residues = []
    i = 0
    while i < len(flattened):
        # Try two-letter codes first
        two = flattened[i : i + 2]
        if two in _RESIDUE_MASS:
            residues.append(two)
            i += 2
        else:
            residues.append(flattened[i])
            i += 1
    return residues


def _parse_monomer(structure: str):
    """
    Parse a monomer structure string into components.

    Returns (glycan_residues, stem_residues, modification_tag, oligo_state)
    or None if the string is not a parseable monomer.
    """
    m = _MONOMER_RE.match(structure.strip())
    if m is None:
        return None
    glycan = _parse_glycan(m.group("glycan"))
    stem = _parse_stem(m.group("stem"))
    mod = m.group("mod") or ""
    oligo = int(m.group("oligo"))
    return glycan, stem, mod, oligo


# ---------------------------------------------------------------------------
# Fragment mass helpers
# ---------------------------------------------------------------------------


def _neutral_fragment(residues: List[str], has_reducing_end: bool = False) -> float:
    """
    Neutral mass of a sub-sequence of residues after all condensation bonds
    within that fragment, before terminus corrections.

    has_reducing_end: True if this fragment contains the reducing end of MurNAc
    (applies the +H2 AUTO_MODS reduction offset).
    """
    if not residues:
        return 0.0
    total = sum(_RESIDUE_MASS[r] for r in residues)
    bonds = (len(residues) - 1) * H2O
    reduction = H2_REDUCTION if has_reducing_end else 0.0
    return total - bonds + reduction


def _mz(neutral: float, charge: int) -> float:
    return (neutral + charge * PROTON) / charge


# ---------------------------------------------------------------------------
# Ion generators
# ---------------------------------------------------------------------------


def _glycan_ions(glycan: List[str], stem: List[str], reduced: bool = True) -> List[dict]:
    """
    B, C, Y and Z ions from glycosidic bond cleavages in the glycan chain.
    Domon-Costello nomenclature (Glycoconj. J. 1988).

    B ions: non-reducing-end fragment, oxonium form (lost OH at cleavage).
    C ions: non-reducing-end fragment, with OH intact (C = B + H2O).
    Y ions: reducing-end fragment + stem, with H restored at cleavage.
    Z ions: reducing-end fragment, without H restored (Z = Y − H2O).

    reduced : bool
        If True (default), adds the NaBH4 reduction offset (+H2) to fragments
        containing the reducing end of MurNAc, matching the pgfinder mass
        convention (AUTO_MODS "Red"). Set to False for non-reduced spectra
        (e.g. Kwan 2024 experimental data).
    """
    ions = []
    n = len(glycan)

    # Full mass of stem piece (attaches to reducing-end m via Stem bond)
    stem_mass_with_stem_bond = (
        (
            sum(_RESIDUE_MASS[r] for r in stem)
            - (len(stem) - 1) * H2O  # peptide bonds within stem
            - H2O  # stem (lactoyl ether) bond connecting m to stem
        )
        if stem
        else 0.0
    )

    for i in range(1, n):  # cut between glycan[i-1] and glycan[i]
        # B ion: glycan[0..i-1], non-reducing end, lost OH (oxonium ion)
        b_residues = glycan[:i]
        b_neutral = _neutral_fragment(b_residues, has_reducing_end=False) - H2O
        label_b = f"B{i}"

        # C ion: same fragment, OH intact → C = B + H2O
        c_neutral = b_neutral + H2O
        label_c = f"C{i}"

        # Y ion: glycan[i..n-1] + stem, reducing end (optionally +H2 reduced)
        y_glycan = glycan[i:]
        y_glycan_neutral = _neutral_fragment(y_glycan, has_reducing_end=reduced)
        # y-terminus: lost H gained H → net 0
        y_neutral = y_glycan_neutral + stem_mass_with_stem_bond
        label_y = f"Y{n - i}"

        # Z ion: same fragment, without the restored H → Z = Y − H2O
        z_neutral = y_neutral - H2O
        label_z = f"Z{n - i}"

        ions.append({"label": label_b, "neutral": b_neutral, "ion_series": "B"})
        ions.append({"label": label_c, "neutral": c_neutral, "ion_series": "C"})
        ions.append({"label": label_y, "neutral": y_neutral, "ion_series": "Y"})
        ions.append({"label": label_z, "neutral": z_neutral, "ion_series": "Z"})

    return ions


def _stem_lac_ions(glycan: List[str], stem: List[str], reduced: bool = True) -> List[dict]:
    """
    b-lac and y-lac ions from cleavage of the lactoyl (Stem) bond between
    MurNAc and the N-terminus of the stem peptide.

    b-lac: full glycan fragment (contains reducing end when reduced=True).
    y-lac: stem peptide + lactoyl group (never contains reducing end; unaffected
    by the reduced flag).
    """
    if not glycan or "m" not in glycan or not stem:
        return []

    # b-lac: full glycan chain, minus the lactyl group and H2O at the cut site
    # b-lac terminal: lost="C3H4O2OHH" = lost C3H4O2 + H2O
    b_lac_neutral = _neutral_fragment(glycan, has_reducing_end=reduced) - _LACTOYL - H2O
    label_b_lac = "b-lac"

    # y-lac: stem peptide with lactyl group at N-terminus
    # y-lac terminal: lost="H" gained="C3H4O2H" → net gain of C3H4O2
    stem_neutral = _neutral_fragment(stem, has_reducing_end=False)
    y_lac_neutral = stem_neutral + _LACTOYL
    label_y_lac = "y-lac"

    return [
        {"label": label_b_lac, "neutral": b_lac_neutral, "ion_series": "b-lac"},
        {"label": label_y_lac, "neutral": y_lac_neutral, "ion_series": "y-lac"},
    ]


def _peptide_ions(stem: List[str]) -> List[dict]:
    """
    b and y ions from peptide bond cleavages within the stem peptide.

    b ions: N-terminal fragments.  b-terminus correction: -H2O.
    y ions: C-terminal fragments.  y-terminus correction: 0.

    b1 is the smallest meaningful b ion (one residue).
    y1 is the smallest meaningful y ion (one residue).
    A b_max or y_max ion equal to the full stem is omitted (that would be
    the intact peptide, not a fragment).
    """
    ions = []
    n = len(stem)

    for i in range(1, n):  # cut after position i (b_i / y_(n-i))
        # b ion: stem[0..i-1]
        b_neutral = _neutral_fragment(stem[:i], has_reducing_end=False) - H2O  # b-terminus: -H2O
        ions.append({"label": f"b{i}", "neutral": b_neutral, "ion_series": "b"})

        # y ion: stem[i..n-1]
        y_neutral = _neutral_fragment(stem[i:], has_reducing_end=False)  # y-terminus: 0 correction
        ions.append({"label": f"y{n - i}", "neutral": y_neutral, "ion_series": "y"})

    return ions


def _immonium_ions(stem: List[str]) -> List[dict]:
    """
    Immonium ions from single-residue internal fragments (b+y type, then -CO).

    These are the smallest diagnostic ions — one residue, with both a b and y
    terminus applied, then CO lost.  formula: residue_free_mass - H2O(b) + 0(y) - CO - H2O(bond).

    Simplified: single residue with b-terminus (-H2O) minus CO = residue - H2O - CO.
    Since there's only one residue, there are no condensation bonds.
    """
    seen = set()
    ions = []
    for r in stem:
        if r in seen:
            continue
        seen.add(r)
        # Single-residue fragment: no condensation bonds
        # b-terminus (-H2O) applied, then -CO for immonium
        neutral = _RESIDUE_MASS[r] - H2O - CO
        ions.append({"label": f"Im({r})", "neutral": neutral, "ion_series": "immonium"})
    return ions


def _secondary_losses(glycan: List[str]) -> List[dict]:
    """
    Generate standalone freed-residue diagnostic ions from GlcNAc and MurNAc.

    These are small, fixed-mass ions that appear in any spectrum containing
    that monosaccharide, independent of the cleavage pattern.
    """
    extra = []
    has_glc = "g" in glycan
    has_mur = "m" in glycan

    if has_glc:
        for label, neutral in zip(_GLC_FREED_LABELS, _GLC_FREED):
            extra.append({"label": label, "neutral": neutral, "ion_series": "glycan-freed"})
    if has_mur:
        for label, neutral in zip(_MUR_FREED_LABELS, _MUR_FREED):
            extra.append({"label": label, "neutral": neutral, "ion_series": "glycan-freed"})

    return extra


# Losses for E/Q residues in y-type fragments (from hcd_rules.kdl).
# residue "E": lost "H2O" (e1), lost "H2OCONH2" (e2)
# residue "Q": lost "NH3" (q1), lost "NH3CONH2" (q2)
_NH3 = _mass("NH3")  # 17.02655
_CONH2 = _mass("CH2NO")  # 44.01380 — CONH2 written as formula CH2NO (C1H2N1O1)
_E_LOSSES = [H2O, H2O + _CONH2]
_E_LABELS = ["e1", "e2"]
_Q_LOSSES = [_NH3, _NH3 + _CONH2]
_Q_LABELS = ["q1", "q2"]


def _eq_losses(stem: List[str], peptide_ions: List[dict]) -> List[dict]:
    """
    Apply E and Q secondary losses to y-type stem peptide ions that contain
    those residues (hcd_rules.kdl residue "E"/"Q" blocks).

    For each y_k ion whose C-terminal k residues include an E or Q, generate:
        E → -H2O (e1) and -(H2O + CONH2) (e2)
        Q → -NH3 (q1) and -(NH3 + CONH2) (q2)
    """
    extra = []
    n = len(stem)
    for ion in peptide_ions:
        label = ion["label"]
        if not label.startswith("y"):
            continue
        try:
            k = int(label[1:])
        except ValueError:
            continue
        # y_k covers the last k residues: stem[n-k..n-1]
        y_residues = stem[n - k :]
        for residue, losses, ion_labels in (
            ("E", _E_LOSSES, _E_LABELS),
            ("Q", _Q_LOSSES, _Q_LABELS),
        ):
            if residue in y_residues:
                for loss, loss_label in zip(losses, ion_labels):
                    neutral = ion["neutral"] - loss
                    extra.append(
                        {
                            "label": f"{label}-{loss_label}",
                            "neutral": neutral,
                            "ion_series": loss_label,
                        }
                    )
    return extra


def _terminal_losses(precursor_neutral: float) -> List[dict]:
    """
    Terminal group losses from the precursor (hcd_rules.kdl group blocks).

    group "Carboxyl" at "C-terminal" lost "H2O"  → precursor - H2O
    group "Amino"    at "N-terminal" lost "NH3"   → precursor - NH3
    """
    return [
        {"label": "precursor-H2O", "neutral": precursor_neutral - H2O, "ion_series": "terminal-loss"},
        {"label": "precursor-NH3", "neutral": precursor_neutral - _NH3, "ion_series": "terminal-loss"},
    ]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def predict_hcd_fragments(
    structure: str,
    charge_max: int = 3,
    reduced: bool = True,
) -> Optional[pd.DataFrame]:
    """
    Predict HCD fragment ions for a muropeptide monomer structure string.

    Parameters
    ----------
    structure : str
        pgfinder structure string, e.g. "gm-AEJA|1" or "gm-AEJAA (Anh)|1".
        Dimers and non-monomer structures are not yet supported and return None.
    charge_max : int
        Maximum charge state to generate.  Default 3 covers most HCD spectra.
        Charges 1..charge_max are generated for each fragment.
    reduced : bool
        If True (default), the reducing-end saccharide carries the NaBH4
        reduction offset (+H2, +2.01565 Da), matching the pgfinder mass
        convention (AUTO_MODS "Red").  Set to False for non-reduced spectra
        such as the Kwan 2024 experimental dataset.

    Returns
    -------
    pd.DataFrame or None
        Columns: ion_type, ion_series, mz, charge, neutral_mass.
        Sorted by mz ascending.
        Returns None if the structure cannot be parsed.

    Notes
    -----
    - Intact precursor mass is included as ion_type "precursor".
    - Adducts beyond protons (Na+, K+) are not generated.
    """
    if "=" in structure:
        return None  # dimers not yet supported

    parsed = _parse_monomer(structure)
    if parsed is None:
        return None

    glycan, stem, *_ = parsed

    # Validate all residue codes are known
    all_residues = glycan + stem
    unknown = [r for r in all_residues if r not in _RESIDUE_MASS]
    if unknown:
        return None

    # Collect all fragment ions
    fragments: List[dict] = []
    fragments.extend(_glycan_ions(glycan, stem, reduced=reduced))
    fragments.extend(_stem_lac_ions(glycan, stem, reduced=reduced))
    peptide_ions = _peptide_ions(stem)
    fragments.extend(peptide_ions)
    fragments.extend(_eq_losses(stem, peptide_ions))
    fragments.extend(_immonium_ions(stem))
    fragments.extend(_secondary_losses(glycan))

    # Intact precursor
    precursor_neutral = (
        _neutral_fragment(glycan, has_reducing_end=reduced)
        + sum(_RESIDUE_MASS[r] for r in stem)
        - (len(stem) - 1) * H2O  # peptide bonds within stem
        - H2O  # stem (lactoyl) bond
    )
    fragments.append({"label": "precursor", "neutral": precursor_neutral, "ion_series": "precursor"})
    fragments.extend(_terminal_losses(precursor_neutral))

    # Expand to charge states
    rows = []
    for frag in fragments:
        neutral = frag["neutral"]
        if neutral <= 0:
            continue
        for z in range(1, charge_max + 1):
            mz_val = _mz(neutral, z)
            if mz_val < 50:  # below instrument detection range
                continue
            rows.append(
                {
                    "ion_type": frag["label"],
                    "ion_series": frag["ion_series"],
                    "mz": round(mz_val, 6),
                    "charge": z,
                    "neutral_mass": round(neutral, 6),
                }
            )

    if not rows:
        return pd.DataFrame(columns=["ion_type", "ion_series", "mz", "charge", "neutral_mass"])

    df = pd.DataFrame(rows)
    df.sort_values("mz", inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df
