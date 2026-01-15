from claspar import virus
import json
import pytest
import pandas as pd


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


def test_process_viral_aligner(viral_aligner_test_data, viral_aligner_test_thresholds):
    filtered_df = virus.process_viral_aligner(
        viral_aligner_test_data, viral_aligner_test_thresholds
    )
    assert (
        filtered_df.shape[0] == 2
    )  # should have only two results that pass the filter


def test_get_viral_aligner_results(
    viral_aligner_test_data, viral_aligner_test_thresholds
):
    headline, results, dfs = virus.get_viral_aligner_results(
        sample_id="ID-12345678",
        original_viral_aligner_results=viral_aligner_test_data,
        viral_aligner_thresholds_dict=viral_aligner_test_thresholds,
    )
    assert (
        "Sample ID-12345678 has 2 viral taxa classified by the Viral Aligner that passed the filters."
        in headline
    )
    print(headline)
