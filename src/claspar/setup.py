import logging
import yaml
import os
from pathlib import Path


def check_filters(filters: list[str], threshold_dict: dict) -> int:
    """
    Checks that the filter threshold dictionary contains all the keys as needed. If there are any extra keys, these are
    logged but exitcode is 0. If any keys are missing, these are logged as an error and exitcode 1 is returned.

    :param filters: list of strings: the filters expected in the filter threshold dictionary.
    :param threshold_dict: the filter threshold dictionary.
    :return: exitcode (int)
    """
    exitcode = 0
    # First check for filters in the threshold dict compared to expected (extras)
    for f1 in threshold_dict:
        if f1 not in filters:
            logging.info(f"Filter '{f1}' is not set to be used.")

    # Second check for filters not in the threshold dict that are expected (missing)
    for f2 in filters:
        if f2 not in threshold_dict:
            logging.error(
                f"Filter {f2} has not been provided and is required. Exiting."
            )
            exitcode = 1
    return exitcode


def read_config_file(config_file: str | os.PathLike) -> tuple[dict, list]:
    """
    Read config file to get thresholds to filter each of the classifier results. Check that all filters are accounted
    for and log any that are not used or missing.
    :param config_file: path to yaml file containing filter thresholds.
    :returns: dict, nested dictionary of thresholds and list of exit codes (list of ints)
    """
    exit_codes = []

    with Path(config_file).open("r") as file:
        thresholds = yaml.safe_load(file)

    expected_kraken_filters = [
        "READ_THRESHOLD",
        "GENUS_RANK_THRESHOLD",
        "GENUS_READ_PCT_THRESHOLD",
    ]
    exit_codes.append(
        check_filters(expected_kraken_filters, thresholds["kraken_bacterial_filters"])
    )

    expected_sylph_filters = [
        "CONTAINMENT_INDEX_THRESHOLD",
        "EFFECTIVE_COVERAGE_THRESHOLD",
    ]
    exit_codes.append(
        check_filters(expected_sylph_filters, thresholds["sylph_filters"])
    )

    expected_viral_aligner_filters = [
        "EVENNESS_VALUE",
        "COVERAGE_1X",
        "UNIQUELY_MAPPED_READS",
        "MEAN_READ_IDENTITY",
        "MEAN_ALIGNMENT_LENGTH",
    ]
    exit_codes.append(
        check_filters(
            expected_viral_aligner_filters, thresholds["viral_aligner_filters"]
        )
    )
    return thresholds, exit_codes
