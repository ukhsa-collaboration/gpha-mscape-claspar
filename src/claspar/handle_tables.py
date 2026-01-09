from onyx_analysis_helper import onyx_analysis_helper_functions as oa

###################
# Analysis tables #


def create_bacterial_analysis_fields(
    *,
    domain: str,
    classifier: str,
    record_id: str,
    thresholds: dict,
    headline_result: str,
    results: dict,
    server: str,
) -> dict:
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
    onyx_analysis = oa.OnyxAnalysis()
    onyx_analysis.add_analysis_details(
        analysis_name=f"{domain}-classifier-parser",
        analysis_description=f"This is an analysis to parse and filter the {domain} classifications from {classifier}",
    )
    onyx_analysis.add_package_metadata(package_name="claspar")
    methods_fail = onyx_analysis.add_methods(methods_dict=thresholds)
    results_fail = onyx_analysis.add_results(top_result=headline_result, results_dict=results)
    onyx_analysis.add_server_records(sample_id=record_id, server_name=server)
    required_field_fail, attribute_fail = onyx_analysis.check_analysis_object(publish_analysis=False)

    if any(  # noqa: SIM108
        [methods_fail, results_fail, required_field_fail, attribute_fail]
    ):  # noqa SIM108
        exitcode = 1
    else:
        exitcode = 0

    return onyx_analysis, exitcode


# def write_qc_results_to_json(
#     qc_dict: dict, sample_id: str, results_dir: os.path
# ) -> os.path:
#     """Write qc results dictionary to json output file.
#     Arguments:
#         qc_dict -- Dictionary containing qc results
#         sample_id -- Sample ID to use in file name
#         results_dir -- Directory to save results to
#     Returns:
#         os.path of saved json file
#     """
#     result_file = Path(results_dir) / f"{sample_id}_qc_results.json"
#
#     with Path(result_file).open("w") as file:
#         json.dump(qc_dict, file)
#
#     return result_file
