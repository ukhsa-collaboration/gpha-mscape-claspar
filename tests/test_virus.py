import json

import pandas as pd
import pytest

from claspar import virus


@pytest.fixture
def metadata_json():
    with open("tests/test_metadata.json") as json_file:  # noqa: PTH123
        metadata = json.loads(json_file.read())
    assert metadata["alignment_results"]
    return metadata


@pytest.fixture
def viral_aligner_test_data(metadata_json):
    return pd.DataFrame(metadata_json["alignment_results"])


@pytest.fixture
def viral_aligner_test_thresholds():
    thresholds = {
        "EVENNESS_VALUE": 25,
        "COVERAGE_1X": 25,
        "UNIQUELY_MAPPED_READS": 0,
        "MEAN_READ_IDENTITY": 90,
        "MEAN_ALIGNMENT_LENGTH": 500,
    }
    return thresholds


def test_filter_viral_aligner(viral_aligner_test_data, viral_aligner_test_thresholds):
    virus_parser = virus.VirusClasPar(
        sample_id="ID-12345678",
        original_viral_aligner_df=viral_aligner_test_data,
        virus_thresholds_dict=viral_aligner_test_thresholds,
        server="testing_server",
    )
    ec, filtered_df = virus_parser._filter_viral_aligner()
    assert ec == 0, f"Expected exitcode to be 0, got {ec}"
    assert (s := filtered_df.shape[0]) == 2, f"Expected 2 results to pass the filter, got {s}"
    print(f"\nFunction filter_viral_aligner with test data correctly gives \n{filtered_df}.")


def test_get_viral_aligner_results(viral_aligner_test_data, viral_aligner_test_thresholds):
    virus_parser = virus.VirusClasPar(
        sample_id="ID-12345678",
        original_viral_aligner_df=viral_aligner_test_data,
        virus_thresholds_dict=viral_aligner_test_thresholds,
        server="testing_server",
    )
    headline, results = virus_parser._get_viral_aligner_results()
    assert (s := "Sample ID-12345678 has 2 viral taxa classified") in headline, (
        f"Expected {s} in headline, got {headline}"
    )
    assert (x := len(results.keys())) == 2, f"Expected 2 keys in the results dict, got {x}"
    assert results[0]["human_readable"] == "rhinovirus A1", (
        "Something has changed about how the results dictionary is made."
    )
    print(f"\nHeadline is: '{headline}'\n and results are: {results}")


def test_get_results_with_no_results(viral_aligner_test_thresholds):
    # If the alignment has been run but there are no results, or it has not been run, both return an empty dataframe
    empty_alignment_df = pd.DataFrame()
    virus_parser = virus.VirusClasPar(
        sample_id="ID-87654321",
        original_viral_aligner_df=empty_alignment_df,
        virus_thresholds_dict=viral_aligner_test_thresholds,
        server="testing_server",
    )
    headline, results = virus_parser._get_viral_aligner_results()
    assert (
        "ID-87654321 has 0 viral taxa classified by the Viral Aligner that passed the filters (out of a total of 0)."
        in headline
    )
    assert results == {}, f"Expected empty dict, got {results}"
    print(f"\nHeadline is '{headline}'\n and results are: {results}")


def test_filter_viral_aligner_no_results(viral_aligner_test_thresholds):
    empty_alignment_df = pd.DataFrame()
    virus_parser = virus.VirusClasPar(
        sample_id="ID-87654321",
        original_viral_aligner_df=empty_alignment_df,
        virus_thresholds_dict=viral_aligner_test_thresholds,
        server="testing_server",
    )
    exitcode, filtered_df = virus_parser._filter_viral_aligner()
    assert exitcode == 0, f"Expected exitcode 0, got {exitcode}"
    assert filtered_df.empty, f"Expected empty dataframe, got {filtered_df}"
    print(f"\nFunction filter_viral_aligner with empty dataframe correctly gives: {filtered_df}.")


def test_get_virus_analysis_table(viral_aligner_test_data, viral_aligner_test_thresholds):
    virus_parser = virus.VirusClasPar(
        sample_id="ID-12345678",
        original_viral_aligner_df=viral_aligner_test_data,
        virus_thresholds_dict=viral_aligner_test_thresholds,
        server="mscape",
    )
    analysis_table = virus_parser.get_virus_analysis_table()
    assert (e := virus_parser.exitcode) == 0, f"Expected exitcode of 0, got {e}"
    assert (p := analysis_table.pipeline_name) == "ClasPar", f'Expected pipeline name "ClasPar", got "{p}"'
    assert (n := analysis_table.name) == "virus-classifier-parser", f'Expected name "virus-classifier-parser", got {n}'


def test_nothing_breaks_if_missing_column(viral_aligner_test_data, viral_aligner_test_thresholds):
    broken_test_data = viral_aligner_test_data.drop("mean_alignment_length", axis=1)
    virus_parser = virus.VirusClasPar(
        sample_id="ID-9999999",
        original_viral_aligner_df=broken_test_data,
        virus_thresholds_dict=viral_aligner_test_thresholds,
        server="mscape",
    )
    assert (e := virus_parser.exitcode) == 1, f"Expected exitcode to be 1 for broken input data, got {e}"
    assert (r := virus_parser.results) == {}, f"Expected empty results dict, got {r}"
    virus_parser.get_virus_analysis_table()
    ec = virus_parser.exitcode
    assert ec == 1, f"Expected exitcode of 1, but got {ec}"
    print(
        "\n",
        "When testing with input data that is missing a column, the exitcode should be 1 before and after getting "
        "the analysis table - success.",
    )
