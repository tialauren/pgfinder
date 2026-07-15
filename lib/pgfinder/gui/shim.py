"""Module for gluing the WebUI together."""

import json
import re

from pgfinder import COLUMNS, matching, pgio, validation
from pgfinder.accurate_species_rules import (
    AccurateSpeciesRules,
    CustomCrosslinkRule,
    match_custom_acceptor,
    match_custom_donor,
)
from pgfinder.gui.internal import (
    MASS_LIB_DIR,
    ms_upload_reader,
    theo_masses_upload_reader,
)

CUSTOM_SPECIES_VALUE = "__custom__"


def mass_library_index():
    return open(MASS_LIB_DIR / "index.json").read()


def species_index():
    return json.dumps(dict(AccurateSpeciesRules.list_available_species()))


def allowed_modifications():
    return validation.allowed_modifications()


def preview_custom_rule():
    """Show example structures matched by a user-supplied donor/acceptor regex pair.

    Lets the WebUI surface a bad pattern (matches nothing, matches the wrong
    things) before the user commits to a full analysis run.
    """
    from pyio import (
        customAcceptorBridgeType,
        customAcceptorPattern,
        customDonorPattern,
        customMaxGlycineBridge,
        customMinGlycineBridge,
        massLibrary,
    )

    try:
        re.compile(customDonorPattern)
        re.compile(customAcceptorPattern)
    except re.error as e:
        return json.dumps({"error": str(e)})

    rule = CustomCrosslinkRule(
        label="preview",
        donor_pattern=customDonorPattern,
        acceptor_pattern=customAcceptorPattern,
        loses_terminal_ala=False,
        acceptor_bridge_type=customAcceptorBridgeType,
        acceptor_min_glycine_bridge=customMinGlycineBridge,
        acceptor_max_glycine_bridge=customMaxGlycineBridge,
    )

    theo_masses = theo_masses_upload_reader(upload=massLibrary.to_py())
    structure_column = COLUMNS["pgfinder"]["inferred"]["structure"]

    donor_matches = []
    acceptor_matches = []
    for structure in theo_masses[structure_column]:
        if match_custom_donor(structure, rule):
            donor_matches.append(structure)
        if match_custom_acceptor(structure, rule):
            acceptor_matches.append(structure)

    return json.dumps(
        {
            "donorMatches": donor_matches[:10],
            "donorCount": len(donor_matches),
            "acceptorMatches": acceptor_matches[:10],
            "acceptorCount": len(acceptor_matches),
            "error": None,
        }
    )


def _custom_rule_from_fields(
    donor_pattern, acceptor_pattern, loses_terminal_ala, acceptor_bridge_type, min_glycine_bridge, max_glycine_bridge
):
    return CustomCrosslinkRule(
        label="Custom",
        donor_pattern=donor_pattern,
        acceptor_pattern=acceptor_pattern,
        loses_terminal_ala=loses_terminal_ala,
        acceptor_bridge_type=acceptor_bridge_type,
        acceptor_min_glycine_bridge=min_glycine_bridge,
        acceptor_max_glycine_bridge=max_glycine_bridge,
    )


def generate_theoretical_dimers():
    """Generate just the theoretical dimer list - the same donor/acceptor
    combinations a full dimer-matching run would search for - without matching
    them back against the raw MS1 data.

    This lets a user preview the dimer search space (and the masses PGFinder
    would actually use) without waiting for a full analysis to complete.
    """
    from pyio import (
        customAcceptorBridgeType,
        customAcceptorPattern,
        customDonorPattern,
        customLosesTerminalAla,
        customMaxGlycineBridge,
        customMinGlycineBridge,
        donorAbundanceThreshold,
        enabledModifications,
        massLibrary,
        msData,
        permissiveMode,
        ppmTolerance,
        species,
    )

    theo_masses = theo_masses_upload_reader(upload=massLibrary.to_py())
    donor_abundance_threshold = donorAbundanceThreshold / 100

    def generate(virt_file):
        ms_data = ms_upload_reader(virt_file)
        if species == CUSTOM_SPECIES_VALUE:
            custom_rule = _custom_rule_from_fields(
                customDonorPattern,
                customAcceptorPattern,
                customLosesTerminalAla,
                customAcceptorBridgeType,
                customMinGlycineBridge,
                customMaxGlycineBridge,
            )
            theoretical_dimers, donors_used = matching.generate_theoretical_dimers(
                ms_data,
                theo_masses,
                enabledModifications,
                ppmTolerance,
                custom_rule=custom_rule,
                donor_abundance_threshold=donor_abundance_threshold,
            )
        else:
            theoretical_dimers, donors_used = matching.generate_theoretical_dimers(
                ms_data,
                theo_masses,
                enabledModifications,
                ppmTolerance,
                species_code=species,
                strict_mode=not permissiveMode,
                donor_abundance_threshold=donor_abundance_threshold,
            )
        return {"dimers": theoretical_dimers.to_csv(index=False), "donors": donors_used.to_csv(index=False)}

    return {f["name"]: generate(f) for f in msData.to_py()}


def run_analysis():
    from pyio import (
        cleanupWindow,
        consolidationPpm,
        customAcceptorBridgeType,
        customAcceptorPattern,
        customDonorPattern,
        customLosesTerminalAla,
        customMaxGlycineBridge,
        customMinGlycineBridge,
        donorAbundanceThreshold,
        enableDimers,
        enabledModifications,
        massLibrary,
        msData,
        permissiveMode,
        ppmTolerance,
        species,
    )

    theo_masses = theo_masses_upload_reader(upload=massLibrary.to_py())
    donor_abundance_threshold = donorAbundanceThreshold / 100

    def analyze(virt_file):
        ms_data = ms_upload_reader(virt_file)
        if enableDimers and species == CUSTOM_SPECIES_VALUE:
            custom_rule = _custom_rule_from_fields(
                customDonorPattern,
                customAcceptorPattern,
                customLosesTerminalAla,
                customAcceptorBridgeType,
                customMinGlycineBridge,
                customMaxGlycineBridge,
            )
            matched, donors_used = matching.data_analysis_with_dimers(
                ms_data,
                theo_masses,
                cleanupWindow,
                enabledModifications,
                ppmTolerance,
                consolidationPpm,
                enable_dimers=True,
                custom_rule=custom_rule,
                donor_abundance_threshold=donor_abundance_threshold,
            )
        elif enableDimers:
            matched, donors_used = matching.data_analysis_with_dimers(
                ms_data,
                theo_masses,
                cleanupWindow,
                enabledModifications,
                ppmTolerance,
                consolidationPpm,
                enable_dimers=True,
                species_code=species,
                strict_mode=not permissiveMode,
                donor_abundance_threshold=donor_abundance_threshold,
            )
        else:
            matched = matching.data_analysis(
                ms_data, theo_masses, cleanupWindow, enabledModifications, ppmTolerance, consolidationPpm
            )
            donors_used = None

        result = {virt_file["name"]: pgio.dataframe_to_csv_metadata(matched)}
        if donors_used is not None:
            basename = virt_file["name"].rsplit(".", 1)[0]
            result[f"{basename}_donors_used.csv"] = donors_used.to_csv(index=False)
        return result

    combined = {}
    for f in msData.to_py():
        combined.update(analyze(f))
    return combined
