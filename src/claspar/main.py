#!/usr/bin/env python3

"""
Example main script that pulls in functions from module(s) and runs the code. Suggested layout
and example helper functions provided below. These can be amended as required.
"""

# Imports - ordered (can use ruff to do this automatically)
import argparse
import logging
import sys
from importlib import resources
from pathlib import Path

from claspar import bacteria, setup, virus


# Arg parse setup
def get_args():
    parser = argparse.ArgumentParser(
        prog="claspar",
        description="""ClasPar: the friendly classifier parser that parses, filters and writes classifier results to 
        analysis tables.""",
    )
    parser.add_argument("--sample_id", "-i", type=str, required=True, help="Sample ID")
    parser.add_argument("--output_dir", "-o", type=str, required=True, help="directory to save results to.")
    parser.add_argument(
        "--config",
        "-c",
        type=str,
        required=False,
        help="Path to yaml file with filtering thresholds",
    )
    parser.add_argument(
        "--server",
        "-s",
        type=str,
        required=True,
        choices=["mscape", "synthscape"],
        help="Specify server code is being run on - helpful if developing on synthscape and running on mscape",
    )
    parser.add_argument(
        "--dry-run",
        "-d",
        required=False,
        action="store_true",
        help="Use this option if code includes a step that writes to onyx so that it can be tested",
    )

    return parser.parse_args()


# Logger set up
def set_up_logger(stdout_file):
    """Example logger set up which can be amended as required. In this example,
    all logging messages go to a stdout log file, and error messages also go to
    stderr log. If the component runs correctly, stderr is empty. The logger is
    set to append mode so logs from older runs are not overwritten.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s")

    out_handler = logging.FileHandler(stdout_file, mode="a")
    out_handler.setFormatter(formatter)
    logger.addHandler(out_handler)

    return logger


# Main function
def main():
    """ """
    # Main function description here

    #########
    # Setup #
    #########
    main_exitcode = 0
    # Retrieve command line arguments:
    args = get_args()  # noqa: F841

    # Set up log file:
    log_file = Path(args.output) / f"{args.input}_claspar_log.txt"
    set_up_logger(log_file)

    # Set up thresholds:
    threshold_dict = {}
    # Use default filtering thresholds yaml file if custom file is not supplied
    if not args.config:
        with resources.as_file(resources.files("claspar.lib").joinpath("filter_thresholds.yaml")) as config_file:
            threshold_dict, exit_codes = setup.read_config_file(config_file)
            if any(exit_codes):
                # logging happens in the function
                main_exitcode = 1
                return main_exitcode

        logging.info("No custom filtering thresholds yaml file specified, using default parameters from included file.")
    else:
        logging.info("Reading filtering thresholds from custom yaml file provided: %s", args.config)

        # Read in filtering thresholds from yaml file
        try:
            threshold_dict, config_exit_codes = setup.read_config_file(args.config)
            # If any of the exit_codes are 1, then some filters were missing - need to exit.
            # These are logged in the check_filters function:
            # --> Exit if any issues with the config:
            if any(config_exit_codes):
                # logging happens in the function
                main_exitcode = 1
                return main_exitcode
        # --> Exit if config file not found:
        except FileNotFoundError:
            logging.error("Specified filtering thresholds yaml file %s not found, exiting program.", args.config)
            main_exitcode = 1
            return main_exitcode

    # Set up data needed (query Onyx once here)
    viral_aligner_input_df, sylph_input_df, classifier_calls_df = setup.get_input_data(args.sample_id, args.server)

    ####################
    # The Actual Thing #
    ####################

    ############
    # Bacteria #
    ############

    ##########
    # Kraken #
    ##########

    # Add in rest of code including logging messages:
    logging.info("Parsing the Kraken Bacteria classifications...")

    # Instantiate the kraken bacterial parser class
    kraken_bacteria_parser = bacteria.KrakenBacteria(
        sample_id=args.sample_id,
        original_classifier_df=classifier_calls_df,
        kraken_bacteria_thresholds_dict=threshold_dict["kraken_bacterial_filters"],
        server=args.server,
    )

    # Get the analysis table:
    kraken_bacterial_analysis_table = kraken_bacteria_parser.get_kraken_bacteria_analysis_table()

    # Check it's all going smoothly in there:
    if kraken_bacteria_parser.exitcode == 1:
        main_exitcode = 1
        return main_exitcode

    # All good so far, let's write analysis table to json:
    kraken_bacteria_json_path = Path(args.output) / f"{args.sample_id}_kraken_bacteria_analysis_fields.json"
    kraken_bacterial_analysis_table.write_analysis_to_json(result_file=kraken_bacteria_json_path)  # type: ignore

    logging.info(
        "Kraken Bacterial Classifications for Onyx analysis fields written to file %s", kraken_bacteria_json_path
    )

    # Write files to csv:
    kraken_bacteria_parser.save_outputs_to_csv(args.output_dir)
    logging.info("All Kraken bacterial species and genera data written to csv in %s", args.output_dir)
    logging.info("Finished parsing Kraken Bacterial results.")

    #########
    # Sylph #
    #########

    # Add in rest of code including logging messages:
    logging.info("Parsing the Sylph classifications...")

    # Instantiate the kraken bacterial parser class
    sylph_parser = bacteria.SylphBacteria(
        sample_id=args.sample_id,
        original_sylph_df=sylph_input_df,
        sylph_bacteria_thresholds_dict=threshold_dict["sylph_filters"],
        server=args.server,
    )

    # Get the analysis table:
    sylph_analysis_table = sylph_parser.get_sylph_analysis_table()

    # Check it's all going smoothly in there:
    if sylph_parser.exitcode == 1:
        main_exitcode = 1
        return main_exitcode

    # All good so far, let's write analysis table to json:
    sylph_json_path = Path(args.output) / f"{args.sample_id}_sylph_analysis_fields.json"
    sylph_analysis_table.write_analysis_to_json(result_file=sylph_json_path)  # type: ignore

    logging.info("Viral Aligner Onyx analysis fields written to file %s", kraken_bacteria_json_path)

    # Write files to csv:
    kraken_bacteria_parser.save_outputs_to_csv(args.output_dir)
    logging.info("All Kraken bacterial species and genera data written to csv in %s", args.output_dir)
    logging.info("Finished parsing Kraken Bacterial results.")

    ###########
    # Viruses #
    ###########

    # Add in rest of code including logging messages:
    logging.info("Parsing the viral aligner classifications.")

    viral_aligner = virus.VirusClasPar(
        sample_id=args.sample_id,
        original_viral_aligner_df=viral_aligner_input_df,
        virus_thresholds_dict=threshold_dict["viral_aligner_filters"],
        server=args.server,
    )

    viral_aligner_analysis_table = viral_aligner.get_virus_analysis_table()

    # Check it's all going smoothly in there:
    if viral_aligner.exitcode == 1:
        main_exitcode = 1
        return main_exitcode

    # All good so far, let's write analysis table to json:
    viral_aligner_json_path = Path(args.output) / f"{args.sample_id}_viral_aligner_analysis_fields.json"
    viral_aligner_analysis_table.write_analysis_to_json(result_file=viral_aligner_json_path)  # type: ignore

    logging.info("Viral Aligner Onyx analysis fields written to file %s", viral_aligner_json_path)

    # Save csv files to go to s3:
    viral_aligner.save_outputs_to_csv(results_dir=args.output_dir)
    logging.info("Viral Aligner filtered data written to csv in %s", args.output_dir)
    logging.info("Finished parsing Viral Aligner results.")
    ########
    # End! #
    ########

    # Finally, we are finished, just return the exitcode.
    return main_exitcode


# Run
if __name__ == "__main__":
    sys.exit(main())
