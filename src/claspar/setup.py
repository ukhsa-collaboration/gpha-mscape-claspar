import json
import logging
import os
from pathlib import Path

import pandas as pd
import yaml
from onyx import OnyxClient, OnyxConfig, OnyxEnv
from onyx_analysis_helper import onyx_analysis_helper_functions as oa

# Set up onyx config
CONFIG = OnyxConfig(
    domain=os.environ[OnyxEnv.DOMAIN],
    token=os.environ[OnyxEnv.TOKEN],
)


@oa.call_to_onyx
def get_input_data(sample_id: str, server: str) -> tuple[int, list[pd.DataFrame]]:
    """
    Get the input data from Onyx. Decorated to handle errors suitably.
    :param sample_id: ID of the sample (climb-id).
    :param server: the server to query.
    :return: tuple of exitcode (int) and dataframes (list). List of dataframes consists of the alignment results,
    the sylph results and the classifier calls. Note that any of these could be empty dataframes!
    """
    with OnyxClient(CONFIG) as client:
        record = client.get(
            project="mscape",
            climb_id="C-0125714289",
            include=["climb_id", "classifier_calls", "alignment_results", "sylph_results"],
        )

    try:
        alignment_results_df = pd.DataFrame(record["alignment_results"])
        sylph_results_df = pd.DataFrame(record["sylph_results"])
        classifier_calls_df = pd.DataFrame(record["classifier_calls"])
        exitcode = 0
    except KeyError as e:
        logging.error(f"Could not find key {e} in Onyx Record. Exiting cleanly.")
        exitcode = 1

    return exitcode, [alignment_results_df, sylph_results_df, classifier_calls_df]


# def read_samplesheet(path_to_samplesheet: os.PathLike | str) -> list[pd.DataFrame]:
#     """
#     NOT FUNCTIONAL
#     """
#     samplesheet_df = pd.read_csv(path_to_samplesheet, sep="\t")
#     record = samplesheet_df.join(samplesheet_df["record_json"].apply(json.loads).apply(pd.Series))

#     alignment_results_df = pd.DataFrame(record["alignment_results"])
#     sylph_results_df = pd.DataFrame(record["sylph_results"])
#     classifier_calls_df = pd.DataFrame(record["classifier_calls"])

#     return [alignment_results_df, sylph_results_df, classifier_calls_df]


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
            logging.info("Filter '%s' is not set to be used.", f1)

    # Second check for filters not in the threshold dict that are expected (missing)
    for f2 in filters:
        if f2 not in threshold_dict:
            logging.error("Filter %s has not been provided and is required. Exiting.", f2)
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
    exit_codes.append(check_filters(expected_kraken_filters, thresholds["kraken_bacterial_filters"]))

    expected_sylph_filters = [
        "CONTAINMENT_INDEX_THRESHOLD",
        "EFFECTIVE_COVERAGE_THRESHOLD",
    ]
    exit_codes.append(check_filters(expected_sylph_filters, thresholds["sylph_filters"]))

    expected_viral_aligner_filters = [
        "EVENNESS_VALUE",
        "COVERAGE_1X",
        "UNIQUELY_MAPPED_READS",
        "MEAN_READ_IDENTITY",
        "MEAN_ALIGNMENT_LENGTH",
    ]
    exit_codes.append(check_filters(expected_viral_aligner_filters, thresholds["viral_aligner_filters"]))
    return thresholds, exit_codes
