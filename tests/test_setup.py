import logging
from importlib import resources

import pandas as pd

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
