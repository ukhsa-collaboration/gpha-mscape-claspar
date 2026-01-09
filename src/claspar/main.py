#!/usr/bin/env python3

"""
Example main script that pulls in functions from module(s) and runs the code. Suggested layout
and example helper functions provided below. These can be amended as required.
"""

# Imports - ordered (can use ruff to do this automatically)
import argparse
import logging
import sys
import yaml
import os
from pathlib import Path
from importlib import resources
from claspar import setup, bacteria

import mscape_template.mscape_functions as mf  # noqa: F401


# Arg parse setup
def get_args():
    parser = argparse.ArgumentParser(
        prog="claspar",
        description="""ClasPar: the friendly classifier parser that parses, filters and writes classifier results to 
        analysis tables.""",
    )
    parser.add_argument("--input", "-i", type=str, required=True, help="Sample ID")
    parser.add_argument(
        "--output", "-o", type=str, required=True, help="Folder to save results to"
    )
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

    # Retrieve command line arguments:
    args = get_args()  # noqa: F841

    # Set up log file:
    log_file = Path(args.output) / f"{args.input}_claspar_log.txt"
    set_up_logger(log_file)

    # Set up thresholds:
    threshold_dict = {}
    # Use default filtering thresholds yaml file if custom file is not supplied
    if not args.config:
        config_path = resources.files("claspar.lib").joinpath("filter_thresholds.yaml")
        logging.info(
            f"No custom filtering thresholds yaml file specified, using default parameters from file: {config_path}"
        )
    else:
        config_path = args.config
        logging.info(
            f"Reading filtering thresholds from custom yaml file provided: {args.config}"
        )

    # Read in filtering thresholds from yaml file
    try:
        threshold_dict, exit_codes = setup.read_config_file(config_path)
        # If any of the exit_codes are 1, then some filters were missing - need to exit.
        # These are logged in the check_filters function:
        if any(exit_codes):
            exitcode = 1
            return exitcode
    except FileNotFoundError:
        logging.error(
            f"Specified filtering thresholds yaml file {config_path} not found, exiting program."
        )
        exitcode = 1
        return exitcode

    ####################
    # The Actual Thing #
    ####################

    # Add in rest of code including logging messages:
    logging.info(
        "mscape template code beginning"
    )  # Example only - add more informative logging messages

    # Additional code in here

    # Write to logs if component finished successfully (or not):
    logging.info("mscape template code successfully completed")


# Run
if __name__ == "__main__":
    sys.exit(main())
