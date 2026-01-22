"""
Create unit tests for modules in the tests/ folder. All functions in a repo should be unit tested and
tests should be run before and after any changes are made.
"""

import json

import pandas as pd
import pytest  # noqa: F401
from claspar import bacteria
from taxaplease import TaxaPlease

pd.set_option("display.max_colwidth", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)


tp = TaxaPlease()


@pytest.fixture(scope="module")
def metadata_json():
    with open("tests/test_metadata.json") as json_file:  # noqa: PTH123
        metadata = json.loads(json_file.read())
    assert metadata["sylph_results"] and metadata["classifier_calls"]
    return metadata


################
# Kraken Tests #
################


class TestKrakenBacteria:
    @pytest.fixture(autouse=True)
    def kraken_thresholds_dict(self):
        self.thresholds = {
            "READ_THRESHOLD": 10,
            "GENUS_RANK_THRESHOLD": 3,
            "GENUS_READ_PCT_THRESHOLD": 20,
        }

    @pytest.fixture(autouse=True)
    def test_kraken_data(self, metadata_json):
        self.test_input_data = pd.DataFrame(metadata_json["classifier_calls"])

    @pytest.fixture(autouse=True)
    def test_instance_1(self, kraken_thresholds_dict, test_kraken_data):
        kraken_assignments = bacteria.KrakenBacteria(
            sample_id="ID-12345678",
            original_classifier_df=self.test_input_data,
            kraken_bacteria_thresholds_dict=self.thresholds,  # type: ignore
            server="mscape",
        )
        self.kraken_class_instance = kraken_assignments

    def test_test_class(self):
        assert (s := self.kraken_class_instance.sample_id) == "ID-12345678", (
            f"Sanity checking the test class setup failed - expected ID-12345678, got {s}"
        )

    @pytest.mark.parametrize(
        "taxonid,isbacteria,parent_id,parent_name",
        [
            (139, True, 64895, "Borreliella"),
            (1410656, True, 859, "Fusobacterium necrophorum"),
            (2696357, False, 2788787, "unclassified Caudoviricetes"),
            (3052230, False, 11102, "Hepacivirus"),
        ],
    )
    def test__get_parent_taxonomy(self, taxonid, isbacteria, parent_id, parent_name):
        actual = self.kraken_class_instance._get_parent_taxonomy(taxonid)
        assert (a := actual[0]) == isbacteria, f"Expected {isbacteria} but got {a}"
        assert (b := actual[1]["taxid"]) == parent_id, f"Expected {parent_id} but got {b}"  # type: ignore
        assert (c := actual[2]) == parent_name, f"Expected {parent_name}, but got {c}"
        print(f"\nParent taxonomy dict for {taxonid} is {actual[1]}, all is as expected.")

    @pytest.mark.parametrize(
        "test,count_descendants,order_in_genus,pct_genus_reads,outcome",
        [
            ("Dominant taxa, high count", 100, 1, 80, "high"),
            ("thresholds", 10, 3, 20, "high"),
            ("low count, many species in genus", 5, 5, 5, "low"),
            ("high count, many species in genus", 1000, 1, 10, "low"),
        ],
    )
    def test__get_kraken_confidence_rating(self, test, count_descendants, order_in_genus, pct_genus_reads, outcome):
        actual = self.kraken_class_instance._get_kraken_confidence_rating(
            count_descendants=count_descendants, order_in_genus=order_in_genus, pct_genus_reads=pct_genus_reads
        )
        assert actual == outcome, f"Expected {outcome}, got {actual}"
        print(
            f"\nTesting '{test}' with {count_descendants} reads, {order_in_genus} rank in the genus, "
            f"{pct_genus_reads} percentage genus reads was given {outcome} confidence rating. All is as expected."
        )

    def test__process_kraken(self):
        actual_species_df, actual_genus_df = self.kraken_class_instance._process_kraken()
        total_species = len(actual_species_df)
        assert total_species == 80, f"Expected 80, got {total_species}"

        high_confidence_species = actual_species_df.loc[actual_species_df["kraken_confidence"] == "high"].shape[0]
        assert high_confidence_species == 6, f"Expected 6, got {high_confidence_species}"

        total_genera = len(actual_genus_df)
        assert total_genera == 56, f"Expected 56, got {total_genera}"

        genera_with_more_than_one_species = actual_genus_df.loc[actual_genus_df["total_species_identified"] > 1].shape[
            0
        ]
        assert genera_with_more_than_one_species == 14, f"Expected 14, got {genera_with_more_than_one_species}"

        print(
            f"\nProcessing kraken results for the test data revealed {total_species} total species, of which "
            f"{high_confidence_species} were high confidence. There were {total_genera} total genera, of which "
            f"{genera_with_more_than_one_species} had more than one species in it. (All as expected)"
        )

    def test_get_kraken_results(self):
        headline_result, main_result, species_df, genus_df = self.kraken_class_instance._get_kraken_results()
        assert "Sample ID-12345678 has 6 high confidence bacterial species classified by Kraken" in headline_result
        assert (x := len(main_result.keys())) == 6, f"Expected 6, got {x}"
        assert all(y := (len(row) == 3 for row in main_result.values())), f"Expected Each row to have 3, got {y}"
        assert (w := species_df.shape[0]) == 80, f"Expected 80, got {w}"
        assert (z := genus_df.shape[0]) == 56, (
            f"Expected 56, got {z}"
        )  # 56 rows in the genus dataframe using test data.
        print(f'Headline result is: "{headline_result}" - as expected.')

    def test_save_outputs_to_csv(self, tmp_path):
        self.kraken_class_instance.save_outputs_to_csv(tmp_path)
        print(f"Saving to {tmp_path}")

    def test_write_to_json(self, tmp_path):
        self.kraken_class_instance.get_kraken_bacteria_analysis_table()

        filename = tmp_path / f"{self.kraken_class_instance.sample_id}_kraken_bacteria_analysis_fields.json"
        self.kraken_class_instance.analysis_table.write_analysis_to_json(filename)  # type: ignore
        print(f"Saving json to {filename}")


class TestNoKrakenBacteria:
    @pytest.fixture(autouse=True)
    def kraken_thresholds_dict(self):
        self.thresholds = {
            "READ_THRESHOLD": 10,
            "GENUS_RANK_THRESHOLD": 3,
            "GENUS_READ_PCT_THRESHOLD": 20,
        }

    @pytest.fixture(autouse=True)
    def test_kraken_data(self):
        self.test_input_data = pd.DataFrame()

    @pytest.fixture(autouse=True)
    def test_instance_no_data(self, kraken_thresholds_dict, test_kraken_data):
        kraken_assignments = bacteria.KrakenBacteria(
            sample_id="ID-00000000",
            original_classifier_df=self.test_input_data,
            kraken_bacteria_thresholds_dict=self.thresholds,  # type: ignore
            server="mscape",
        )
        self.kraken_class_instance = kraken_assignments

    def test_test_class(self):
        assert (s := self.kraken_class_instance.sample_id) == "ID-00000000", (
            f"Sanity checking the test class setup failed - expected ID-0000000, got {s}"
        )

    def test__process_kraken(self):
        species_df, genus_df = self.kraken_class_instance._process_kraken()
        assert species_df.empty, f"Expected empty species df, got {species_df}"
        assert genus_df.empty, f"Expected empty genus df, got {genus_df}"

    def test__get_kraken_results(self):
        headline, result, species_df, genus_df = self.kraken_class_instance._get_kraken_results()
        assert "Sample ID-00000000 has 0 high confidence bacterial" in headline
        assert "out of 0 total assignments" in headline
        assert result == {}, f"Expected empty dict, got {result}"
        assert species_df.empty, f"Expected empty species df, got {species_df}"
        assert genus_df.empty, f"Expected empty genus df, got {genus_df}"

    def test_save_outputs_to_csv(self, tmp_path):
        self.kraken_class_instance.save_outputs_to_csv(tmp_path)
        print(f"Saving to {tmp_path}")

    def test_write_to_json(self, tmp_path):
        self.kraken_class_instance.get_kraken_bacteria_analysis_table()

        filename = tmp_path / f"{self.kraken_class_instance.sample_id}_kraken_bacteria_analysis_fields.json"
        self.kraken_class_instance.analysis_table.write_analysis_to_json(filename)  # type: ignore
        print(f"Saving json to {filename}")


###############
# Sylph Tests #
###############


class TestSylphBacteria:
    @pytest.fixture(autouse=True)
    def sylph_test_data(self, metadata_json):
        self.sylph_test_df = pd.DataFrame(metadata_json["sylph_results"])

    @pytest.fixture(autouse=True)
    def sylph_thresholds_dict(self):
        self.thresholds = {"CONTAINMENT_INDEX_THRESHOLD": 0.2, "EFFECTIVE_COVERAGE_THRESHOLD": 1.0}

    @pytest.fixture(autouse=True)
    def sylph_test_data_with_rank(self, sylph_test_data):
        sylph_test_data = self.sylph_test_df.copy()
        sylph_test_data["taxon_rank"] = sylph_test_data["taxon_id"].apply(lambda x: tp.get_record(x)["rank"])  # type: ignore
        assert all(sylph_test_data["taxon_rank"] == "species")
        self.sylph_df_with_rank = sylph_test_data

    @pytest.fixture(autouse=True)
    def test_instance_1(self, sylph_test_data, sylph_thresholds_dict):
        sylph_assignments = bacteria.SylphBacteria(
            sample_id="ID-12345678",
            original_sylph_df=self.sylph_test_df,
            sylph_bacteria_thresholds_dict=self.thresholds,
            server="mscape",
        )
        self.sylph_class_instance = sylph_assignments

    def test_test_class(self):
        assert (s := self.sylph_class_instance.sample_id) == "ID-12345678", (
            f"Sanity checking the test class setup failed - expected ID-12345678, got {s}"
        )

    def test__process_sylph_rank(self):
        expected_new_columns = pd.DataFrame(
            [
                [1869212, "Chitinophagaceae bacterium"],
                [520, "Bordetella pertussis"],
                [2104, "Mycoplasmoides pneumoniae"],
            ],
        )

        actual_new_columns = self.sylph_df_with_rank.apply(
            lambda x: self.sylph_class_instance._process_sylph_rank(x), axis=1, result_type="expand"
        )
        assert actual_new_columns.equals(expected_new_columns)
        print(f"\nChecking the sylph rank to get ID and species name - get {actual_new_columns}")

    def test__process_sylph_rank_strain(self):
        test_df = pd.DataFrame(
            {"taxon_id": [1121296], "human_readable": ["[Clostridium] aminophilum DSM 10710"], "taxon_rank": ["strain"]}
        )
        expected_new_columns = pd.DataFrame([[1526, "[Clostridium] aminophilum"]])

        actual_new_columns = test_df.apply(
            lambda x: self.sylph_class_instance._process_sylph_rank(x), axis=1, result_type="expand"
        )
        assert actual_new_columns.equals(expected_new_columns)
        print(
            f"\nTaxa with strain rank will be given the species level ID and human readable:"
            f"{test_df} would give {expected_new_columns}."
        )

    def test__process_sylph_rank_no_id(self):
        test_df = pd.DataFrame({"taxon_id": [""], "human_readable": [""], "taxon_rank": [""]})
        expected_new_columns = pd.DataFrame([[None, None]])

        actual_new_columns = test_df.apply(
            lambda x: self.sylph_class_instance._process_sylph_rank(x), axis=1, result_type="expand"
        )
        assert actual_new_columns.equals(expected_new_columns)
        print(
            f"\nTaxa without a rank will be given the species level ID and human readable:"
            f"{test_df} would give {expected_new_columns}."
        )

    @pytest.mark.parametrize(
        "test,containment_index,effective_coverage,outcome",
        [
            ("high containment, high coverage", 1.0, 10.0, "high"),
            ("thresholds", 0.2, 1.0, "high"),
            ("low containment, high coverage", 0.05, 0.5, "low"),
            ("low containment, low coverage", 0.1, 0.04, "low"),
        ],
    )
    def test__get_sylph_confidence_rating(self, test, containment_index, effective_coverage, outcome):
        actual_outcome = self.sylph_class_instance._get_sylph_confidence_rating(
            containment_index=containment_index, effective_coverage=effective_coverage
        )

        assert actual_outcome == outcome
        print(
            f"\nTest '{test}' with containment index {containment_index} and effective coverage {effective_coverage} "
            f"has confidence rating {outcome} (as expected)."
        )

    def test_process_sylph(self):
        actual_sylph_result = self.sylph_class_instance._process_sylph()
        species_count = actual_sylph_result.loc[actual_sylph_result["taxon_rank"] == "species"].shape[0]
        assert species_count == 3
        high_confidence_sylph = actual_sylph_result.loc[actual_sylph_result["sylph_confidence"] == "high"].shape[0]
        assert high_confidence_sylph == 3
        count_of_rows_with_genus_ids = actual_sylph_result.loc[actual_sylph_result["genus_id"].notnull()].shape[0]
        # Genus ID is null if there is no genus for the species, for example if it's unclassified.
        assert count_of_rows_with_genus_ids == 2
        print(
            f"\nProcessing sylph data - there were {len(actual_sylph_result)} sylph results, of which {species_count}"
            f" were species, and {high_confidence_sylph} had a high confidence. There were "
            f"{count_of_rows_with_genus_ids} classifications with a genus ID (all as expected)."
        )

    def test__get_sylph_results(self):
        headline, main_result, full_result_df = self.sylph_class_instance._get_sylph_results()
        assert (
            "Sample ID-12345678 has 3 high confidence bacterial (and archaeal) species classified by Sylph" in headline
        )
        assert len(main_result.keys()) == 3
        assert (len(row) == 3 for row in main_result.values())  # 3 rows of 3 columns
        assert full_result_df.shape == (3, 22)  # full dataframe to be published has 22 columns
        print(f"\nTest Sylph headline results: {headline} \nand the results:\n {main_result}")

    def test_get_sylph_analysis_table(self):
        analysis_table = self.sylph_class_instance.get_sylph_analysis_table()
        assert (p := analysis_table.pipeline_name) == "ClasPar", f'Expected pipeline name "ClasPar", got "{p}"'
        assert (n := analysis_table.name) == "bacteria-classifier-parser", (
            f'Expected name "bacteria-classifier-parser", got {n}'
        )
        assert "sylph" in (d := analysis_table.description), f'Expected "sylph" to be in the descritopn, got {d}'

    def test_save_outputs_to_csv(self, tmp_path):
        self.sylph_class_instance.save_outputs_to_csv(tmp_path)
        print(f"Saving to {tmp_path}")

    def test_write_to_json(self, tmp_path):
        self.sylph_class_instance.get_sylph_analysis_table()

        filename = tmp_path / f"{self.sylph_class_instance.sample_id}_sylph_analysis_fields.json"
        self.sylph_class_instance.analysis_table.write_analysis_to_json(filename)  # type: ignore
        print(f"Saving json to {filename}")


class TestNoSylphBacteria:
    @pytest.fixture(autouse=True)
    def sylph_test_data(self):
        self.sylph_test_df = pd.DataFrame()

    @pytest.fixture(autouse=True)
    def sylph_thresholds_dict(self):
        self.thresholds = {"CONTAINMENT_INDEX_THRESHOLD": 0.2, "EFFECTIVE_COVERAGE_THRESHOLD": 1.0}

    @pytest.fixture(autouse=True)
    def test_instance_1(self, sylph_test_data, sylph_thresholds_dict):
        sylph_assignments = bacteria.SylphBacteria(
            sample_id="ID-00000000",
            original_sylph_df=self.sylph_test_df,
            sylph_bacteria_thresholds_dict=self.thresholds,
            server="mscape",
        )
        self.sylph_class_instance = sylph_assignments

    def test_test_class(self):
        assert (s := self.sylph_class_instance.sample_id) == "ID-00000000", (
            f"Sanity checking the test class setup failed - expected ID-00000000, got {s}"
        )

    def test_process_sylph(self):
        actual_sylph_result = self.sylph_class_instance._process_sylph()
        assert actual_sylph_result.empty, f"Expected empty dataframe, got {actual_sylph_result}"

    def test_save_outputs_to_csv(self, tmp_path):
        self.sylph_class_instance.save_outputs_to_csv(tmp_path)
        print(f"Saving to {tmp_path}")

    def test_write_to_json(self, tmp_path):
        self.sylph_class_instance.get_sylph_analysis_table()

        filename = tmp_path / f"{self.sylph_class_instance.sample_id}_sylph_analysis_fields.json"
        self.sylph_class_instance.analysis_table.write_analysis_to_json(filename)  # type: ignore
        print(f"Saving json to {filename}")
