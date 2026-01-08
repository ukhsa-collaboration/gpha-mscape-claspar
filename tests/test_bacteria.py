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


@pytest.fixture
def metadata_json():
    with open("tests/test_metadata.json") as json_file:  # noqa: PTH123
        metadata = json.loads(json_file.read())
    assert metadata["sylph_results"] and metadata["classifier_calls"]
    return metadata


@pytest.fixture
def kraken_thresholds_dict():
    return {
        "READ_THRESHOLD": 10,
        "GENUS_RANK_THRESHOLD": 3,
        "GENUS_READ_PCT_THRESHOLD": 20,
    }


@pytest.fixture
def sylph_thresholds_dict():
    return {"CONTAINMENT_INDEX_THRESHOLD": 0.2, "EFFECTIVE_COVERAGE_THRESHOLD": 1.0}


@pytest.fixture
def test_kraken_data(metadata_json):
    return pd.DataFrame(metadata_json["classifier_calls"])


@pytest.fixture
def test_sylph_data(metadata_json):
    return pd.DataFrame(metadata_json["sylph_results"])


def test_process_collection_date():
    data = {
        "id": [1, 2, 3],
        "collection_date": ["2025-10-10", "", "1992-12-31"],
        "received_date": ["", "12-12-2025", ""],
    }
    df = pd.DataFrame(data)
    new_df = bacteria._process_collection_date(df)
    assert new_df.shape == (3, 6), "Actual dataframe does not have expected shape."
    assert (
        "collection_date_epi_week" in new_df.columns.values
    ), "Actual dataframe does not contain collection_date_epi_week column."
    print(f"Metadata with date handling looks like: {new_df}")


@pytest.mark.parametrize(
    "taxonid,isbacteria,parent_id,parent_name",
    [
        (139, True, 64895, "Borreliella"),
        (1410656, True, 859, "Fusobacterium necrophorum"),
        (2696357, False, 2788787, "unclassified Caudoviricetes"),
        (3052230, False, 11102, "Hepacivirus"),
    ],
)
def test__get_parent_taxonomy(taxonid, isbacteria, parent_id, parent_name):
    actual = bacteria._get_parent_taxonomy(taxonid, tp)
    assert actual[0] == isbacteria
    assert actual[1]["taxid"] == parent_id
    assert actual[2] == parent_name
    print(f"Parent taxonomy for {taxonid} is {actual}")


@pytest.mark.parametrize(
    "test,count_descendants,order_in_genus,pct_genus_reads,outcome",
    [
        ("Dominant taxa, high count", 100, 1, 80, "high"),
        ("thresholds", 10, 3, 20, "high"),
        ("low count, many species in genus", 5, 5, 5, "low"),
        ("high count, many species in genus", 1000, 1, 10, "low"),
    ],
)
def test__get_kraken_confidence_rating(
    test,
    count_descendants,
    order_in_genus,
    pct_genus_reads,
    outcome,
    kraken_thresholds_dict,
):
    actual = bacteria._get_kraken_confidence_rating(
        count_descendants=count_descendants,
        order_in_genus=order_in_genus,
        pct_genus_reads=pct_genus_reads,
        kraken_thresholds_dict=kraken_thresholds_dict,
    )
    assert actual == outcome
    print(
        f"\nTesting '{test}' with {count_descendants} reads, {order_in_genus} rank in the genus, "
        f"{pct_genus_reads} percentage genus reads was given {outcome} confidence rating."
    )


def test_process_kraken(test_kraken_data, kraken_thresholds_dict):
    actual_species_df, actual_genus_df = bacteria.process_kraken(
        test_kraken_data, kraken_thresholds_dict, tp
    )
    total_species = len(actual_species_df)
    assert total_species == 80

    high_confidence_species = actual_species_df.loc[
        actual_species_df["kraken_confidence"] == "high"
    ].shape[0]
    assert high_confidence_species == 6

    total_genera = len(actual_genus_df)
    assert total_genera == 56

    genera_with_more_than_one_species = actual_genus_df.loc[
        actual_genus_df["total_species_identified"] > 1
    ].shape[0]
    assert genera_with_more_than_one_species == 14

    print(
        f"\nProcessing kraken results for the test data revealed {total_species} total species, of which "
        f"{high_confidence_species} were high confidence. There were {total_genera} total genera, of which "
        f"{genera_with_more_than_one_species} had more than one species in it. (All as expected)"
    )


def test__process_sylph_rank(test_sylph_data):
    expected_new_columns = pd.DataFrame(
        [
            [1869212, "Chitinophagaceae bacterium"],
            [520, "Bordetella pertussis"],
            [2104, "Mycoplasmoides pneumoniae"],
        ],
    )
    print(expected_new_columns)
    test_sylph_data["taxon_rank"] = test_sylph_data["taxon_id"].apply(
        lambda x: tp.get_record(x)["rank"]
    )
    assert all(test_sylph_data["taxon_rank"] == "species")
    actual_new_columns = test_sylph_data.apply(
        lambda x: bacteria._process_sylph_rank(x), axis=1, result_type="expand"
    )
    assert actual_new_columns.equals(expected_new_columns)


@pytest.mark.parametrize(
    "test,containment_index,effective_coverage,outcome",
    [
        ("high containment, high coverage", 1.0, 10.0, "high"),
        ("thresholds", 0.2, 1.0, "high"),
        ("low containment, high coverage", 0.05, 0.5, "low"),
        ("low containment, low coverage", 0.1, 0.04, "low"),
    ],
)
def test__get_sylph_confidence_rating(
    test, containment_index, effective_coverage, outcome, sylph_thresholds_dict
):
    actual_outcome = bacteria._get_sylph_confidence_rating(
        containment_index=containment_index,
        effective_coverage=effective_coverage,
        sylph_thresholds_dict=sylph_thresholds_dict,
    )

    assert actual_outcome == outcome
    print(
        f"\nTest '{test}' with containment index {containment_index} and effective coverage {effective_coverage} "
        f"has confidence rating {outcome} (as expected)."
    )


def test_process_sylph(test_sylph_data, sylph_thresholds_dict):
    actual_sylph_result = bacteria.process_sylph(
        test_sylph_data, sylph_thresholds_dict, tp
    )
    species_count = actual_sylph_result.loc[
        actual_sylph_result["taxon_rank"] == "species"
    ].shape[0]
    assert species_count == 3
    high_confidence_sylph = actual_sylph_result.loc[
        actual_sylph_result["sylph_confidence"] == "high"
    ].shape[0]
    assert high_confidence_sylph == 3
    count_of_rows_with_genus_ids = actual_sylph_result.loc[
        actual_sylph_result["genus_id"].notnull()
    ].shape[0]
    # Genus ID is null if there is no genus for the species, for example if it's unclassified.
    assert count_of_rows_with_genus_ids == 2
    print(
        f"\nProcessing sylph data - there were {len(actual_sylph_result)} sylph results, of which {species_count} were "
        f"species, and {high_confidence_sylph} had a high confidence. There were {count_of_rows_with_genus_ids} "
        f"classifications with a genus ID (all as expected)."
    )
