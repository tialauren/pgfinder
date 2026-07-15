#!/usr/bin/env python3
import argparse as arg
import ast
import importlib.resources as pkg_resources
import logging
import warnings
from pathlib import Path

import yaml

from pgfinder.errors import UserError
from pgfinder.logs.logs import LOGGER_NAME, setup_logger
from pgfinder.matching import data_analysis, data_analysis_with_dimers, generate_theoretical_dimers
from pgfinder.pgio import (
    dataframe_to_csv_metadata,
    ms_file_reader,
    read_yaml,
    theo_masses_reader,
)
from pgfinder.utils import update_config

LOGGER = setup_logger()
LOGGER = logging.getLogger(LOGGER_NAME)


def create_parser() -> arg.ArgumentParser:
    """Create a parser for reading options."""
    parser = arg.ArgumentParser(
        description="Process sample. Additional arguments over-ride those in the configuration file."
    )
    parser.add_argument(
        "-c", "--config_file", dest="config_file", required=False, help="Path to a YAML configuration file."
    )
    parser.add_argument("--input_file", dest="input_file", required=False, help="Input File")
    parser.add_argument("--ppm_tolerance", dest="ppm_tolerance", type=float, required=False, help="PPM Toleraance.")
    parser.add_argument(
        "--consolidation_ppm",
        dest="consolidation_ppm",
        type=float,
        required=False,
        help="Maximum absolute ppm distance between consolidated structures.",
    )
    parser.add_argument("--masses_file", dest="masses_file", type=str, required=False, help="Theoretical masses file.")
    parser.add_argument("--time_delta", dest="time_delta", type=int, required=False, help="Time delta.")
    parser.add_argument(
        "--mod_list", dest="mod_list", type=ast.literal_eval, required=False, help="Modifications to include."
    )
    parser.add_argument("--output_dir", dest="output_dir", type=str, required=False, help="Output directory.")
    parser.add_argument("--warnings", dest="warnings", type=str, required=False, help="Whether to ignore warnings.")
    parser.add_argument("--quiet", dest="quiet", type=bool, required=False, help="Supress output.")
    parser.add_argument(
        "--float_format", dest="float_format", type=int, required=False, help="Decimal places in output."
    )

    #  Dimer matching options
    parser.add_argument(
        "--enable_dimers", dest="enable_dimers", action="store_true", help="Enable dimer matching (requires --species)"
    )
    parser.add_argument(
        "--species",
        dest="species",
        type=str,
        required=False,
        choices=["ecoli", "saureus", "efaecalis", "bsubtilis", "cdiff", "fusobacterium"],
        help="Bacterial species for dimer matching (e.g., saureus, ecoli)",
    )
    parser.add_argument(
        "--permissive_mode",
        dest="permissive_mode",
        action="store_true",
        help="Use permissive mode for dimer matching (allows bridge variants)",
    )
    parser.add_argument(
        "--donor_abundance_threshold",
        dest="donor_abundance_threshold",
        type=float,
        required=False,
        help="Percentage (0-100) of the eligible donor pool's cumulative intensity to cover when selecting donors",
    )

    return parser


def process_file(
    input_file,
    masses_file,
    ppm_tolerance: float,
    consolidation_ppm: float,
    time_delta: int,
    mod_list: list,
    output_dir,
    float_format: int = 2,
    to_csv: dict = None,
    enable_dimers: bool = False,
    species: str = None,
    permissive_mode: bool = False,
    donor_abundance_threshold: float = 90,
):
    """Process files

    Parameters
    ----------
    input_file : str | Path
        Mass Spectrometry input file to process.
    masses_file : str | Path
        Input file of known masses.
    mod_list : list
        Modifications to include.
    ppm_tolerance : float
        Parts Per Million tolerance for matching.
    time_delta : int
        Time difference.
    output_dir : str | Path
        Output directory where results are written to.
    float_format : int
       Decimal places to use in CSV files.
    to_csv: dict
       Dictionary of options to pass to pd.to_csv(), primarly used to overwrite existing files.
    enable_dimers : bool
       Enable dimer matching (requires species)
    species : str
       Bacterial species for dimer matching
    permissive_mode : bool
       Use permissive mode for dimer matching (allows bridge variants)
    donor_abundance_threshold : float
       Percentage (0-100) of the eligible donor pool's cumulative intensity to cover when selecting donors
    """
    input_file = Path(input_file)
    masses_file = Path(masses_file)
    output_dir = Path(output_dir)

    df = ms_file_reader(input_file)
    masses = theo_masses_reader(masses_file)

    LOGGER.info(f"PPM Tolerance                      : {ppm_tolerance}")
    LOGGER.info(f"Time Delta                         : {time_delta}")

    if enable_dimers:
        if species is None:
            raise UserError("--species is required when using --enable_dimers")

        LOGGER.info("Dimer matching                     : ENABLED")
        LOGGER.info(f"Species                            : {species}")
        LOGGER.info(f"Mode                               : {'PERMISSIVE' if permissive_mode else 'STRICT'}")
        LOGGER.info(f"Donor abundance threshold          : {donor_abundance_threshold}%")

        # Generate theoretical dimers (detects monomers internally to pick donors)
        LOGGER.info(f"Generating theoretical dimers for {species}...")
        theoretical_dimers_df, donors_used_df = generate_theoretical_dimers(
            raw_data_df=df,
            theo_masses_df=masses,
            enabled_mod_list=mod_list,
            ppm_tolerance=ppm_tolerance,
            species_code=species,
            strict_mode=not permissive_mode,
            donor_abundance_threshold=donor_abundance_threshold / 100,
        )

        # Output theoretical dimers database
        input_basename = input_file.stem
        theoretical_dimer_filename = f"{input_basename}_theoretical_dimers.csv"
        theoretical_dimers_df.to_csv(
            output_dir / theoretical_dimer_filename, index=False, float_format=f"%.{float_format}f"
        )
        LOGGER.info(f"Theoretical dimers database saved to: {output_dir}/{theoretical_dimer_filename}")
        LOGGER.info(f"Generated {len(theoretical_dimers_df)} theoretical dimers")

        # Output the donors actually used to build them
        donors_used_filename = f"{input_basename}_donors_used.csv"
        donors_used_df.to_csv(output_dir / donors_used_filename, index=False, float_format=f"%.{float_format}f")
        LOGGER.info(f"Donors used database saved to: {output_dir}/{donors_used_filename}")

        #  full dimer analysis
        results, _ = data_analysis_with_dimers(
            raw_data_df=df,
            theo_masses_df=masses,
            rt_window=time_delta,
            enabled_mod_list=mod_list,
            ppm_tolerance=ppm_tolerance,
            consolidation_ppm=consolidation_ppm,
            enable_dimers=True,
            species_code=species,
            strict_mode=not permissive_mode,
            donor_abundance_threshold=donor_abundance_threshold / 100,
        )

        monomers = results[results["Inferred structure"].str.contains(r"\|1$", na=False, regex=True)]
        dimers = results[results["Inferred structure"].str.contains(r"\|2$", na=False, regex=True)]
        LOGGER.info(f"Monomers found                     : {len(monomers)}")
        LOGGER.info(f"Dimers found                       : {len(dimers)}")
    else:
        LOGGER.info("Dimer matching                     : DISABLED")
        results = data_analysis(
            raw_data_df=df,
            theo_masses_df=masses,
            rt_window=time_delta,
            enabled_mod_list=mod_list,
            ppm_tolerance=ppm_tolerance,
            consolidation_ppm=consolidation_ppm,
        )

    LOGGER.info("Processing complete!")

    input_basename = input_file.stem

    if enable_dimers:
        # I have put in 2 options here either the user can have everything together
        filename = f"{input_basename}_with_dimers.csv"
        dataframe_to_csv_metadata(
            save_filepath=output_dir,
            output_dataframe=results,
            filename=filename,
            float_format=f"%.{float_format}f",
        )
        LOGGER.info(f"Combined results saved to: {output_dir}/{filename}")

        # Save dimers separately
        dimers_only = results[results["Inferred structure"].str.contains(r"\|2$", na=False, regex=True)]
        if not dimers_only.empty:
            dimers_filename = f"{input_basename}_dimers_only.csv"
            dataframe_to_csv_metadata(
                save_filepath=output_dir,
                output_dataframe=dimers_only,
                filename=dimers_filename,
                float_format=f"%.{float_format}f",
            )
            LOGGER.info(f"Dimers-only file saved to: {output_dir}/{dimers_filename}")

        # Count results
        monomers = results[results["Inferred structure"].str.contains(r"\|1$", na=False, regex=True)]
        dimers = results[results["Inferred structure"].str.contains(r"\|2$", na=False, regex=True)]
        LOGGER.info(f"Monomers found: {len(monomers)}")
        LOGGER.info(f"Dimers found: {len(dimers)}")
    else:
        # Monomero nly mode
        filename = f"{input_basename}.csv"
        dataframe_to_csv_metadata(
            save_filepath=output_dir,
            output_dataframe=results,
            filename=filename,
            float_format=f"%.{float_format}f",
        )
        LOGGER.info(f"Results saved to: {output_dir}/{filename}")


def main():
    """Run processing."""
    try:

        parser = create_parser()
        args = parser.parse_args()

        if args.config_file is not None:
            config = read_yaml(args.config_file)
            LOGGER.info(f"Configuration file loaded from     : {args.config_file}")
        else:
            default_config = pkg_resources.open_text(__package__, "default_config.yaml")
            config = yaml.safe_load(default_config.read())
            LOGGER.info("Default configuration file loaded.")

        config = update_config(config, args)

        if config["warnings"] == "ignore":
            warnings.filterwarnings("ignore")
            LOGGER.info("NB : All warnings have been turned off for this run.")
        elif config["warnings"] == "deprecated":

            def fxn():
                warnings.warn("deprecated", DeprecationWarning, stacklevel=2)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fxn()

        if config["quiet"]:
            LOGGER.setLevel("ERROR")

        process_file(
            input_file=config["input_file"],
            masses_file=config["masses_file"],
            ppm_tolerance=config["ppm_tolerance"],
            consolidation_ppm=config["consolidation_ppm"],
            time_delta=config["time_delta"],
            mod_list=config["mod_list"],
            output_dir=config["output_dir"],
            float_format=config["float_format"],
            enable_dimers=config.get("enable_dimers", False),
            species=config.get("species", None),
            permissive_mode=config.get("permissive_mode", False),
            donor_abundance_threshold=config.get("donor_abundance_threshold", 90),
        )

    except UserError as e:

        LOGGER.error(e)


if __name__ == "__main__":
    main()
