import logging
from importlib import resources

import pandas as pd
import pytest

from claspar import setup

pd.set_option("display.max_colwidth", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)

LOGGER = logging.getLogger(__name__)


def test_read_config_file():
    with resources.as_file(resources.files("claspar.data").joinpath("filter_thresholds.yaml")) as config_file:
        thresholds, exitcodes = setup.read_config_file(config_file)
    assert exitcodes == [0, 0, 0]
    assert thresholds["sylph_filters"]
    print(f"\nConfig read in using function correctly - \n{thresholds}. Exit codes are {exitcodes}")


def test_check_filters():
    filters_to_check = ["up", "down"]
    thresholds_to_check = {"up": 10, "down": 10}
    check = setup.check_filters(filters_to_check, thresholds_to_check)
    assert check == 0  # If all the filters needed matches the config, all is well


def test_check_filters_extra_filter(caplog):
    filters_to_check = ["up", "down"]
    thresholds_to_check = {"up": 10, "down": 10, "left": 10}
    with caplog.at_level(logging.INFO):
        check = setup.check_filters(filters_to_check, thresholds_to_check)
    assert "Filter 'left' is not set to be used" in caplog.text
    assert check == 0  # exit code 0 if dicts match.
    print(f"If there is an extra filter in the config that is not used in the code, this is logged: {caplog.text}.")


def test_check_filters_missing_filter(caplog):
    filters_to_check = ["up", "down", "right"]
    thresholds_to_check = {"up": 10, "down": 10}
    with caplog.at_level(logging.INFO):
        check = setup.check_filters(filters_to_check, thresholds_to_check)
    assert "Filter 'right' has not been provided and is required. Exiting."
    assert check == 1


def test_parse_samplesheet():
    path_to_samplesheet = "tests/test_samplesheet.tsv"

    exit_code, dfs = setup.read_samplesheet(path_to_samplesheet)
    virus, sylph, kraken = dfs
    assert exit_code == 0
    assert (a := virus.shape[0]) == 35, f"Expected 35, got {a}"
    assert (b := sylph.shape[0]) == 3, f"Expected 3, got {b}"
    assert (c := kraken.shape[0]) == 491, f"Expected 491, got {c}"
    print(f"\nTesting the samplesheet parsing function - read in 3 dataframes, such as {sylph}")


def test_broken_samplesheet(tmp_path, caplog):
    samplesheet = pd.read_csv("tests/test_samplesheet.tsv", sep="\t")
    broken_samplesheet = samplesheet.rename(columns={"full_Onyx_json": "json"})
    broken_samplesheet_path = tmp_path / "broken_samplesheet_1.tsv"
    broken_samplesheet.to_csv(broken_samplesheet_path, index=False, sep="\t")

    with caplog.at_level(logging.DEBUG):
        exitcode, dfs = setup.read_samplesheet(broken_samplesheet_path)
        assert (s := "Expected column 'full_Onyx_json' not found in the samplesheet") in caplog.text, (
            f"Expected error message {s} but got {caplog.text}"
        )
        virus, sylph, kraken = dfs
        assert exitcode == 0, f"Expected exitcode 0, got {exitcode}"
        assert (a := virus.shape[0]) == 35, f"Expected 35, got {a}"
        assert (b := sylph.shape[0]) == 3, f"Expected 3, got {b}"
        assert (c := kraken.shape[0]) == 491, f"Expected 491, got {c}"
        print(f"Checking samplesheet parsing with broken samplesheet. Get log message {caplog.text} as expected.")


def test_one_column_samplesheet(tmp_path, caplog):
    samplesheet = pd.read_csv("tests/test_samplesheet.tsv", sep="\t")
    broken_samplesheet_2 = samplesheet.rename(columns={"full_Onyx_json": "json"})
    broken_samplesheet_2 = broken_samplesheet_2.drop("climb_id", axis=1)
    broken_samplesheet_2_path = tmp_path / "broken_samplesheet_1.tsv"
    broken_samplesheet_2.to_csv(broken_samplesheet_2_path, index=False, sep="\t")

    with caplog.at_level(logging.ERROR):
        exitcode, dfs = setup.read_samplesheet(broken_samplesheet_2_path)
        assert exitcode == 1, f"Expected exitcode 1, got {exitcode}"
        assert (
            lm := "Cannot parse the json from the sample sheet, dataframe has 1 row and 1 columns"
        ) in caplog.text, f"Expected log message {lm}, got {caplog.text}"
        assert all(df.empty for df in dfs), f"Expected all dataframes to be empty in list, got {dfs}"
        print(
            f"\nChecking samplesheet parsing with another broken samplesheet. Log message: \n{caplog.text} as expected."
        )
