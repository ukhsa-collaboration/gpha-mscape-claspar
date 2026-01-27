import json

import pandas as pd
import pytest
from claspar import virus


@pytest.fixture(scope="module")
def metadata_json():
    with open("tests/test_metadata.json") as json_file:  # noqa: PTH123
        metadata = json.loads(json_file.read())
    assert metadata["alignment_results"]
    return metadata


class TestViralParser:
    @pytest.fixture(autouse=True)
    def viral_aligner_test_data(self, metadata_json):
        self.test_input_data = pd.DataFrame(metadata_json["alignment_results"])

    @pytest.fixture(autouse=True)
    def viral_aligner_test_thresholds(self):
        thresholds = {
            "EVENNESS_VALUE": 25,
            "COVERAGE_1X": 25,
            "UNIQUELY_MAPPED_READS": 0,
            "MEAN_READ_IDENTITY": 90,
            "MEAN_ALIGNMENT_LENGTH": 500,
        }
        self.thresholds = thresholds

    @pytest.fixture(autouse=True)
    def test_instance_1(self, viral_aligner_test_thresholds, viral_aligner_test_data):
        self.test_class_instance = virus.VirusClasPar(
            sample_id="ID-12345678",
            original_viral_aligner_df=self.test_input_data,
            virus_thresholds_dict=self.thresholds,  # type: ignore
            server="mscape",
        )

    def test_test_class(self):
        assert (
            s := self.test_class_instance.sample_id
        ) == "ID-12345678", f"Sanity checking the test class setup failed - expected ID-12345678, got {s}"

    def test_filter_viral_aligner(self):
        ec, filtered_df = self.test_class_instance._filter_viral_aligner()
        assert ec == 0, f"Expected exitcode to be 0, got {ec}"
        assert (s := filtered_df.shape[0]) == 2, f"Expected 2 results to pass the filter, got {s}"
        print(f"\nFunction filter_viral_aligner with test data correctly gives \n{filtered_df}.")

    def test_get_viral_aligner_results(self):
        headline, results = self.test_class_instance._get_viral_aligner_results()
        assert (
            s := "Sample ID-12345678 has 2 viral taxa classified"
        ) in headline, f"Expected {s} in headline, got {headline}"
        assert (x := len(results.keys())) == 2, f"Expected 2 keys in the results dict, got {x}"
        assert (
            results[0]["human_readable"] == "rhinovirus A1"
        ), "Something has changed about how the results dictionary is made."
        print(f"\nHeadline is: '{headline}'\n and results are: {results}")

    def test_get_virus_analysis_table(self):
        analysis_table = self.test_class_instance.get_virus_analysis_table()
        assert (e := self.test_class_instance.exitcode) == 0, f"Expected exitcode of 0, got {e}"
        assert (p := analysis_table.pipeline_name) == "ClasPar", f'Expected pipeline name "ClasPar", got "{p}"'
        assert (
            n := analysis_table.name
        ) == "virus-classifier-parser", f'Expected name "virus-classifier-parser", got {n}'
        print(f"Resulting analysis table is formatted correctly: {analysis_table}")

    def test_save_outputs_to_csv(self, tmp_path):
        self.test_class_instance.save_outputs_to_csv(tmp_path)
        print(f"Saving to {tmp_path}")

    def test_write_to_json(self, tmp_path):
        self.test_class_instance.get_virus_analysis_table()
        filename = tmp_path / f"{self.test_class_instance.sample_id}_viral_aligner_analysis_fields.json"
        self.test_class_instance.analysis_table.write_analysis_to_json(filename)  # type: ignore
        print(f"Saving json to {filename}")


class TestNoVirus:
    @pytest.fixture(autouse=True)
    def viral_aligner_test_thresholds(self):
        thresholds = {
            "EVENNESS_VALUE": 25,
            "COVERAGE_1X": 25,
            "UNIQUELY_MAPPED_READS": 0,
            "MEAN_READ_IDENTITY": 90,
            "MEAN_ALIGNMENT_LENGTH": 500,
        }
        self.thresholds = thresholds

    @pytest.fixture(autouse=True)
    def test_instance_1(self, viral_aligner_test_thresholds):
        self.test_no_data_instance = virus.VirusClasPar(
            sample_id="ID-87654321",
            original_viral_aligner_df=pd.DataFrame(),
            virus_thresholds_dict=self.thresholds,  # type: ignore
            server="mscape",
        )

    def test_get_results_with_no_results(self):
        # If the alignment has been run but there are no results, or it has not been run, both return an empty dataframe
        headline, results = self.test_no_data_instance._get_viral_aligner_results()
        assert (
            "ID-87654321 has 0 viral taxa classified by the Viral Aligner that passed the filters (out of a total of 0)."
            in headline
        )
        assert results == {}, f"Expected empty dict, got {results}"
        print(f"\nHeadline is '{headline}'\n and results are: {results}")

    def test_filter_viral_aligner_no_results(self):
        exitcode, filtered_df = self.test_no_data_instance._filter_viral_aligner()
        assert exitcode == 0, f"Expected exitcode 0, got {exitcode}"
        assert filtered_df.empty, f"Expected empty dataframe, got {filtered_df}"
        print(f"\nFunction filter_viral_aligner with empty dataframe correctly gives: {filtered_df}.")

    def test_save_outputs_to_csv(self, tmp_path):
        self.test_no_data_instance.save_outputs_to_csv(tmp_path)
        print(f"Saving broken data csvs {tmp_path}")

    def test_write_to_json(self, tmp_path):
        self.test_no_data_instance.get_virus_analysis_table()

        filename = tmp_path / f"{self.test_no_data_instance.sample_id}_no_viral_aligner_analysis_fields.json"
        self.test_no_data_instance.analysis_table.write_analysis_to_json(filename)  # type: ignore
        print(f"Saving json to {filename}")


class TestBrokenInput:
    @pytest.fixture(autouse=True)
    def viral_aligner_test_thresholds(self):
        thresholds = {
            "EVENNESS_VALUE": 25,
            "COVERAGE_1X": 25,
            "UNIQUELY_MAPPED_READS": 0,
            "MEAN_READ_IDENTITY": 90,
            "MEAN_ALIGNMENT_LENGTH": 500,
        }
        self.thresholds = thresholds

    @pytest.fixture(autouse=True)
    def viral_aligner_broken_test_data(self, metadata_json):
        viral_aligner_test_data = pd.DataFrame(metadata_json["alignment_results"])
        broken_test_data = viral_aligner_test_data.drop("mean_alignment_length", axis=1)
        self.broken_input_data = broken_test_data

    @pytest.fixture(autouse=True)
    def test_instance_1(self, viral_aligner_test_thresholds, viral_aligner_broken_test_data):
        self.test_broken_data_instance = virus.VirusClasPar(
            sample_id="ID-9999999",
            original_viral_aligner_df=self.broken_input_data,
            virus_thresholds_dict=self.thresholds,  # type: ignore
            server="mscape",
        )

    def test_test_class(self):
        assert (
            s := self.test_broken_data_instance.sample_id
        ) == "ID-9999999", f"Sanity checking the test class setup failed - expected ID-9999999, got {s}"

    def test_nothing_breaks_if_missing_column(self):
        assert (
            e := self.test_broken_data_instance.exitcode
        ) == 1, f"Expected exitcode to be 1 for broken input data, got {e}"
        assert (r := self.test_broken_data_instance.results) == {}, f"Expected empty results dict, got {r}"
        self.test_broken_data_instance.get_virus_analysis_table()
        ec = self.test_broken_data_instance.exitcode
        assert ec == 1, f"Expected exitcode of 1, but got {ec}"
        print(
            "\n",
            "When testing with input data that is missing a column, the exitcode should be 1 before and after getting "
            "the analysis table - success.",
        )

    def test_save_outputs_to_csv(self, tmp_path):
        self.test_broken_data_instance.save_outputs_to_csv(tmp_path)
        print(f"Saving broken data csvs {tmp_path}")

    def test_write_to_json(self, tmp_path):
        self.test_broken_data_instance.get_virus_analysis_table()

        filename = tmp_path / f"{self.test_broken_data_instance.sample_id}_broken_viral_aligner_analysis_fields.json"
        self.test_broken_data_instance.analysis_table.write_analysis_to_json(filename)  # type: ignore
        print(f"Saving json to {filename}")
