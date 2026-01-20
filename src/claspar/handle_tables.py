import os
from pathlib import Path

import pandas as pd
from onyx_analysis_helper import onyx_analysis_helper_functions as oa

###################
# Analysis tables #


def create_analysis_fields(
    *,
    domain: str,
    classifier: str,
    record_id: str,
    thresholds: dict[str, int],
    headline_result: str,
    results: dict,
    server: str,
) -> tuple[oa.OnyxAnalysis, int]:
    """
    Set up fields dictionary used to populate analysis table containing
    ClasPar outputs.
    Arguments:
        domain -- str, one of 'bacteria', 'virus', 'fungi' etc
        classifier -- the type of classifier being reported in the table (kraken or sylph)
        record_id -- Climb ID for sample
        thresholds -- Dictionary containing criteria used to filter
        headline_result -- Short description of main result
        results -- Dictionary containing results
        server -- Server code is running on, one of "mscape" or "synthscape"
    Returns:
        onyx_analysis -- Class containing required fields for input to onyx
                         analysis table
        exitcode -- Exit code for checks - will be 0 if all checks passed, 1 if any checks failed
    """
    onyx_analysis = oa.OnyxAnalysis()  # set up class
    # Add analysis details
    onyx_analysis.add_analysis_details(
        analysis_name=f"{domain}-classifier-parser",
        analysis_description=f"This is an analysis to parse and filter the {domain} classifications from {classifier}",
    )
    # Add metadata about the pipeline/package
    onyx_analysis.add_package_metadata(package_name="claspar")
    # Check that the methods were parsed by the class
    methods_fail = onyx_analysis.add_methods(methods_dict=thresholds)
    # Check that the results were parsed by the class
    results_fail = onyx_analysis.add_results(top_result=headline_result, results_dict=results)
    # Add info about sample and server (mscape/synthscape)
    onyx_analysis.add_server_records(sample_id=record_id, server_name=server)
    # Check the final object using the helper method
    required_field_fail, attribute_fail = onyx_analysis.check_analysis_object(publish_analysis=False)
    # If any fail, raise exit code.
    if any(  # noqa: SIM108
        [methods_fail, results_fail, required_field_fail, attribute_fail]
    ):  # noqa SIM108
        exitcode = 1
    else:
        exitcode = 0

    return onyx_analysis, exitcode


def write_df_to_csv(*, df: pd.DataFrame, filename: str, results_dir: str | os.PathLike) -> os.PathLike:
    """
    Write results dataframe to csv.
    :param df: dataframe of results to save.
    :param filename: str, unique name of file WITHOUT extension. (csv gets added).
    :param results_dir: Directory to save results to.
    :returns: os.path of saved csv file.
    """

    result_file_path = Path(results_dir) / f"{filename}.csv"

    df.to_csv(result_file_path)

    return result_file_path
