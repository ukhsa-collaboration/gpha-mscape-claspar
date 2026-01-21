#!/usr/bin/env python3

"""
Example main script that pulls in functions from module(s) and runs the code. Suggested layout
and example helper functions provided below. These can be amended as required.
"""

# Imports - ordered (can use ruff to do this automatically)
import argparse
import logging
import sys
from datetime import datetime
from importlib import resources
from pathlib import Path

from claspar import bacteria, setup, virus

today = datetime.today().strftime("%Y-%m-%d-%H%M")


# Arg parse setup
def get_args():
    parser = argparse.ArgumentParser(
        prog="claspar",
        description="""
        ClasPar: the friendly classifier parser that parses, filters and writes classifier results to analysis tables.
        """,
    )
    parser.add_argument("--sample_id", "-i", dest="sample_id", type=str, required=True, help="Climb-ID for sample.")
    parser.add_argument(
        "--output_dir",
        "-o",
        dest="output_dir",
        type=str,
        required=True,
        help="Path to directory where results will be saved to. Directory will be created if it does not exist.",
    )
    parser.add_argument(
        "--config",
        "-c",
        dest="config",
        type=str,
        required=False,
        help="Optional - Path to yaml file with filtering thresholds",
    )
    parser.add_argument(
        "--server",
        "-s",
        dest="server",
        type=str,
        required=True,
        choices=["mscape", "synthscape"],
        help="Specify server code is being run on - helpful if developing on synthscape and running on mscape",
    )
    parser.add_argument(
        "--samplesheet",
        "-t",
        dest="samplesheet_path",
        type=str,
        required=False,
        help="Optional - Path to samplesheet. Must be tsv, should have header 'full_Onyx_json' or 2 columns, 1 row.",
    )
    parser.add_argument(
        "--log-file",
        "-l",
        dest="log_file",
        type=str,
        required=False,
        help=(
            """
            Optional - Path to log file. Default will be a file called '/sample-id/_claspar_/date-time/.log' in the 
            output directory (where sample_id is the climb-id for the sample and /date-time/ is the date and time of 
            running).
            """
        ),
    )

    return parser.parse_args()


# Logger set up
def set_up_logger(log_filepath):
    """Example logger set up which can be amended as required. In this example,
    all logging messages go to a stdout log file, and error messages also go to
    stderr log. If the component runs correctly, stderr is empty. The logger is
    set to append mode so logs from older runs are not overwritten.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s")

    out_handler = logging.FileHandler(log_filepath, mode="a")
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

    # Set up output dir:
    setup.setup_outdir(args.output_dir)

    # Set up log file:
    log_file = (
        Path(args.log_file) if args.log_file else Path(args.output_dir) / f"{args.sample_id}_{today}_claspar_log.txt"
    )
    set_up_logger(log_file)

    # Set up thresholds:
    threshold_dict = {}
    # Use default filtering thresholds yaml file if custom file is not supplied
    if not args.config:
        with resources.as_file(resources.files("claspar.data").joinpath("filter_thresholds.yaml")) as config_file:
            threshold_dict, exit_codes = setup.read_config_file(config_file)
            if any(exit_codes):
                # logging happens in the function
                # --> Exit if any issues with the config:
                main_exitcode = 1
                return main_exitcode

        logging.info("No custom filtering thresholds yaml file specified, using default parameters from included file.")

    else:
        logging.info("Reading the filtering thresholds from custom yaml file provided: %s" % (args.config))  # noqa

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
            logging.error("Specified filtering thresholds yaml file %s not found, exiting program." % (args.config))  # noqa
            main_exitcode = 1
            return main_exitcode

    # Set up data here - samplesheet or onyx:
    if args.samplesheet_path:
        samplesheet_exitcode, dataframes = setup.read_samplesheet(args.samplesheet_path)
        if samplesheet_exitcode == 1:
            # --> Exit if issues with samplesheet parsing:
            logging.error("Cannot parse samplesheet provided, exiting.")
            main_exitcode = 1
            return main_exitcode

    else:
        # Set up data needed (query Onyx once here)
        onyx_exitcode, dataframes = setup.get_input_data(args.sample_id, args.server)
        if onyx_exitcode == 1:
            # --> Exit if issues with Onyx:
            logging.error("Exiting due to issues with Onyx.")
            main_exitcode = 1
            return main_exitcode

    # Unpack dataframes into variables:
    try:
        viral_aligner_input_df, sylph_input_df, classifier_calls_df = dataframes
    except ValueError as e:
        # --> Exit if cannot unpack dataframes:
        logging.error("Cannot unpack the dataframes into the variables: %s" % (e))  # noqa
        main_exitcode = 1
        return main_exitcode

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
    kraken_bacteria_json_path = Path(args.output_dir) / f"{args.sample_id}_kraken_bacteria_analysis_fields.json"
    kraken_bacterial_analysis_table.write_analysis_to_json(result_file=kraken_bacteria_json_path)  # type: ignore

    logging.info(
        "Kraken Bacterial Classifications for Onyx analysis fields written to file %s" % (kraken_bacteria_json_path)  # noqa
    )

    # Write files to csv:
    kraken_bacteria_parser.save_outputs_to_csv(args.output_dir)
    logging.info("All Processed Kraken bacterial species and genera data written to csv in %s" % (args.output_dir))  # noqa
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
    sylph_json_path = Path(args.output_dir) / f"{args.sample_id}_sylph_analysis_fields.json"
    sylph_analysis_table.write_analysis_to_json(result_file=sylph_json_path)  # type: ignore

    logging.info("Sylph Bacterial Classifications for Onyx analysis fields written to file %s" % (sylph_json_path))  # noqa

    # Write files to csv:
    sylph_parser.save_outputs_to_csv(args.output_dir)
    logging.info("All processed Sylph data written to csv in %s" % (args.output_dir))  # noqa

    logging.info("Finished parsing Sylph results.")

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
    viral_aligner_json_path = Path(args.output_dir) / f"{args.sample_id}_viral_aligner_analysis_fields.json"
    viral_aligner_analysis_table.write_analysis_to_json(result_file=viral_aligner_json_path)  # type: ignore

    logging.info("Viral Aligner Onyx analysis fields written to file %s", viral_aligner_json_path)

    # Save csv files to go to s3:
    viral_aligner.save_outputs_to_csv(results_dir=args.output_dir)
    logging.info("Viral Aligner filtered data written to csv in %s" % (args.output_dir))  # noqa
    logging.info("Finished parsing Viral Aligner results.")

    ########
    # End! #
    ########

    # Finally, we are finished, just return the exitcode.
    return main_exitcode


# Run
if __name__ == "__main__":
    sys.exit(main())
