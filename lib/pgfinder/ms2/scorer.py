"""
Cosine similarity scorer for MS2 spectra.

Matches an experimental peak list against the theoretical fragment ions
predicted by fragmentor.predict_hcd_fragments and returns a score plus
an annotated peak table.

Scoring model
-------------
A modified cosine similarity is used, standard in metabolomics/proteomics
(e.g. used by GNPS, MS-DIAL, Sirius):

    score = Σ(sqrt(I_exp_i) × sqrt(I_theo_i))
            / (sqrt(Σ I_exp_i) × sqrt(Σ I_theo_i))

where the sum runs over MATCHED peaks only, and I_theo is uniform (1.0 per
theoretical fragment) since we have no predicted relative intensities yet.

Using sqrt(intensity) down-weights dominant peaks and prevents a single high
peak from dominating the score — this is the standard "square-root cosine"
used in GNPS/MassBank matching.

Tolerance
---------
Default tolerance is 20 ppm (appropriate for Orbitrap HCD data).  A minimum
absolute tolerance of 0.02 Da is also applied so that very low m/z ions are
not penalised by an unrealistically tight window.

Matching is greedy: each experimental peak is matched to at most one
theoretical peak (the closest one within tolerance).  Each theoretical peak
can also be matched at most once.
"""

from __future__ import annotations

import math
from typing import Optional

import pandas as pd

from pgfinder.ms2.fragmentor import predict_hcd_fragments


def _ppm_window(mz: float, ppm: float, min_da: float) -> float:
    """Return the larger of the ppm-derived tolerance and the minimum Da window."""
    return max(mz * ppm / 1e6, min_da)


def score_spectrum(
    structure: str,
    experimental_peaks: pd.DataFrame,
    ppm_tolerance: float = 20.0,
    min_da_tolerance: float = 0.02,
    charge_max: int = 3,
    reduced: bool = True,
) -> Optional[dict]:
    """
    Match a structure's predicted HCD fragments against an experimental MS2
    spectrum and return a cosine similarity score.

    Parameters
    ----------
    structure : str
        pgfinder structure string, e.g. "gm-AEJA|1".
    experimental_peaks : pd.DataFrame
        Must have columns "mz" and "intensity".  Intensities are normalised
        internally so their absolute scale does not matter.
    ppm_tolerance : float
        Mass tolerance in ppm for peak matching.  Default 20 ppm (Orbitrap).
    min_da_tolerance : float
        Minimum absolute tolerance in Da.  Prevents over-strict windows at
        very low m/z.  Default 0.02 Da.
    charge_max : int
        Maximum charge state to predict (passed to predict_hcd_fragments).

    Returns
    -------
    dict or None
        Keys:
            score           float 0–1  cosine similarity (sqrt-weighted)
            matched_peaks   int        number of experimental peaks matched
            total_theo      int        total theoretical fragment peaks
            coverage        float 0–1  fraction of theoretical peaks matched
            annotated       pd.DataFrame  experimental peaks with matched
                            ion_type, ion_series, theo_mz columns added;
                            unmatched rows have NaN in those columns.
        Returns None if the structure cannot be fragmented (parse failure).
    """
    theo = predict_hcd_fragments(structure, charge_max=charge_max, reduced=reduced)
    if theo is None or theo.empty:
        return None

    exp = experimental_peaks[["mz", "intensity"]].copy()
    if exp.empty:
        return None

    # Normalise experimental intensities to max=1 so scale doesn't matter
    exp = exp.copy()
    exp["intensity"] = exp["intensity"] / exp["intensity"].max()

    # Sort both by m/z for greedy matching
    theo_sorted = theo.sort_values("mz").reset_index(drop=True)
    exp_sorted = exp.sort_values("mz").reset_index(drop=True)

    # --- Greedy nearest-neighbour matching, theoretical → experimental ---
    # We iterate over THEORETICAL peaks and find the best experimental match
    # for each. This prevents a noise peak at e.g. 204.067 from consuming the
    # theoretical B1 at 204.087 before the real high-intensity peak at 204.086
    # gets a chance to claim it.
    theo_mz = theo_sorted["mz"].to_numpy()
    exp_mz = exp_sorted["mz"].to_numpy()

    matched_exp_idx = set()
    matched_theo_idx = set()

    matches = []  # (exp_idx, theo_idx, delta_ppm)

    e = 0  # pointer into exp_mz
    for t_idx, t_mz in enumerate(theo_mz):
        tol = _ppm_window(t_mz, ppm_tolerance, min_da_tolerance)

        # Advance exp pointer to first peak in range
        while e < len(exp_mz) and exp_mz[e] < t_mz - tol:
            e += 1

        # Find the closest unmatched experimental peak within tolerance
        best_dist = float("inf")
        best_e = -1
        j = e
        while j < len(exp_mz) and exp_mz[j] <= t_mz + tol:
            if j not in matched_exp_idx:
                dist = abs(exp_mz[j] - t_mz)
                if dist < best_dist:
                    best_dist = dist
                    best_e = j
            j += 1

        if best_e >= 0:
            matched_exp_idx.add(best_e)
            matched_theo_idx.add(t_idx)
            delta_ppm = (exp_mz[best_e] - t_mz) / t_mz * 1e6
            matches.append((best_e, t_idx, delta_ppm))

    # --- Score ---
    # sqrt-weighted cosine: Σ(√I_exp × √I_theo) / (√Σ(I_exp) × √N_theo_total)
    # I_theo = 1.0 for every theoretical peak
    numerator = sum(math.sqrt(exp_sorted.iloc[e]["intensity"]) for e, _, _ in matches)
    denom_exp = math.sqrt(sum(exp_sorted["intensity"]))
    denom_theo = math.sqrt(len(theo_sorted))
    score = numerator / (denom_exp * denom_theo) if (denom_exp * denom_theo) > 0 else 0.0

    coverage = len(matched_theo_idx) / len(theo_sorted) if len(theo_sorted) > 0 else 0.0

    # --- Annotated output ---
    # Add ion annotations to the experimental peak table
    annotated = exp_sorted.copy()
    annotated["ion_type"] = None
    annotated["ion_series"] = None
    annotated["theo_mz"] = float("nan")
    annotated["delta_ppm"] = float("nan")

    for e_idx, t_idx, dppm in matches:
        annotated.at[e_idx, "ion_type"] = theo_sorted.iloc[t_idx]["ion_type"]
        annotated.at[e_idx, "ion_series"] = theo_sorted.iloc[t_idx]["ion_series"]
        annotated.at[e_idx, "theo_mz"] = theo_sorted.iloc[t_idx]["mz"]
        annotated.at[e_idx, "delta_ppm"] = round(dppm, 2)

    return {
        "score": round(score, 6),
        "matched_peaks": len(matches),
        "total_theo": len(theo_sorted),
        "coverage": round(coverage, 4),
        "annotated": annotated.reset_index(drop=True),
    }


def score_structures(
    candidates: list[str],
    experimental_peaks: pd.DataFrame,
    ppm_tolerance: float = 20.0,
    min_da_tolerance: float = 0.02,
    charge_max: int = 3,
    reduced: bool = True,
) -> pd.DataFrame:
    """
    Score multiple candidate structures against one experimental MS2 spectrum.

    Parameters
    ----------
    candidates : list[str]
        List of pgfinder structure strings to score.
    experimental_peaks : pd.DataFrame
        Must have columns "mz" and "intensity".
    ppm_tolerance, min_da_tolerance, charge_max
        Passed through to score_spectrum.

    Returns
    -------
    pd.DataFrame
        Columns: structure, score, matched_peaks, total_theo, coverage.
        Sorted by score descending.
    """
    rows = []
    for structure in candidates:
        result = score_spectrum(
            structure,
            experimental_peaks,
            ppm_tolerance=ppm_tolerance,
            min_da_tolerance=min_da_tolerance,
            charge_max=charge_max,
            reduced=reduced,
        )
        if result is None:
            continue
        rows.append(
            {
                "structure": structure,
                "score": result["score"],
                "matched_peaks": result["matched_peaks"],
                "total_theo": result["total_theo"],
                "coverage": result["coverage"],
            }
        )

    if not rows:
        return pd.DataFrame(columns=["structure", "score", "matched_peaks", "total_theo", "coverage"])

    df = pd.DataFrame(rows).sort_values("score", ascending=False).reset_index(drop=True)
    return df
